import numpy as np

from . import kernels

from .. import kdtree

from typing import Optional, Union


class SPH2D(kdtree.KDTree2D):

    """
    Class for computing Smooth Particle Hydrodynamic (SPH) density and field estimation in 2D.
    """

    def __init__(self) -> None:
        """
        Initialises SPH in 2D.
        """
        kdtree.KDTree2D.__init__(self)
        self.kernel_type = 'cubic'


    def assign_mass(self, mass: Optional[np.ndarray] = None) -> None:
        """
        Assign particles mass, if none they will be assigned.
        """
        if mass is None:
            self.mass = np.ones(len(self.points))
        else:
            self.mass = mass


    def kernel(self, r: Union[float, np.ndarray], h: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Returns the SPH kernel value.

        Parameters
        ----------
        r : float or array
            Radius from a given neighbour.
        h : float or array
            Smoothing length.
        """
        if self.kernel_type == 'cubic':
            return kernels.cubic_kernel(r, h, dim=2)


    def setup(self, k: int = 20, mass: Optional[np.ndarray] = None) -> None:
        """
        Set basic settings for the SPH data.

        Parameters
        ----------
        k : int, optional
            Neighbours for nearest neighbour SPH calculation.
        mass : array, optional
            For assigning masses to the particles.
        """
        self.k = k
        self.assign_mass(mass=mass)


    def sph_estimate(
        self, x: Union[float, np.ndarray], y: Union[float, np.ndarray], f: Optional[np.ndarray] = None, dens: Optional[np.ndarray] = None
    ) -> Union[float, np.ndarray]:
        """
        Estimates a field based on SPH k neighbours.

        Parameters
        ----------
        x, y : float or array
            Cartesian coordinates.
        f : array
            Field values at KDTree points.
        dens : array
            Density estimation.

        Returns
        -------
        f_est : float or array
            Field estimation.
        """
        nind, ndist = self.nearest(x, y, k=self.k, return_dist=True)
        h = np.max(ndist, axis=1)
        if dens is None:
            w = np.array([self.mass[nind[i]]*self.kernel(ndist[i], h[i]) for i in range(0, len(h))])
            dens = np.sum(w, axis=1)
        if f is None:
            return dens
        else:
            # calculating density at each points is slow so instead calculate field and then divide by density.
            w = np.array([self.mass[nind[i]]*(self.f[nind[i]])*self.kernel(ndist[i], h[i]) for i in range(0, len(h))])
            f_est = np.sum(w, axis=1)
            f_est /= dens
            return f_est


    def get_density(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        """
        Calculates density based on SPH k neighbours.

        Parameters
        ----------
        x, y : array
            Cartesian coordinates.

        Returns
        -------
        dens : array
            Density estimation.
        """
        dens = self.sph_estimate(x, y, f=None)
        return dens


    def set_field(self, f: np.ndarray) -> None:
        """
        Sets the field values for SPH field estimation.

        Parameters
        ----------
        f : array
            Field values at KDTree points.
        """
        self.f = f


    def estimate(self, x: Union[float, np.ndarray], y: Union[float, np.ndarray], dens: np.ndarray = None) -> Union[float, np.ndarray]:
        """
        Estimates field at points given.

        Parameters
        ----------
        x, y : float or array
            Cartesian coordinates.
        dens : array
            Density at cartesian coordinates, if None this will be computed but will
            require double the time.

        Returns
        -------
        f_est : float or array
            Field estimation.
        """
        return self.sph_estimate(x, y, f=self.f, dens=dens)


    def clean(self) -> None:
        """
        Reinitialises the class.
        """
        self.__init__()
