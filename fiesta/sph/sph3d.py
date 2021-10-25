import numpy as np

from . import kernels

from .. import coords
from .. import kdtree


class SPH3D(kdtree.KDTree3D):

    """Class for computing Smooth Particle Hydrodynamic (SPH) density and field
    estimation in 3D.
    """

    def __init__(self):
        """Initialises SPH in 3D."""
        kdtree.KDTree3D.__init__(self)
        self.kernel_type = 'cubic'


    def assign_mass(self, mass=None):
        """Assign particles mass, if none they will be assigned"""
        if mass is None:
            self.mass = np.ones(len(self.points))
        else:
            self.mass = mass


    def kernel(self, r, h):
        """Returns the SPH kernel value.

        Parameters
        ----------
        r : float/array
            Radius from a given neighbour.
        h : float
            Smoothing length.
        """
        if self.kernel_type == 'cubic':
            return kernels.cubic_kernel(r, h, dim=3)


    def setup(self, k=50, mass=None):
        """Set basic settings for the SPH data.

        Parameters
        ----------
        k : int, optional
            Neighbours for nearest neighbour SPH calculation.
        mass : array, optional
            For assigning masses to the particles.
        """
        self.k = k
        self.assign_mass(mass=mass)


    def sph_estimate(self, x, y, z, f=None, dens=None):
        """Estimates a field based on SPH k neighbours.

        Parameters
        ----------
        x, y, z : array
            Cartesian coordinates.

        Returns
        -------
        f_est : array
            Field estimation.
        """
        nind, ndist = self.nearest(x, y, z, k=self.k, return_dist=True)
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


    def get_density(self, x, y, z):
        """Calculates density based on SPH k neighbours.

        Parameters
        ----------
        x, y, z : array
            Cartesian coordinates.

        Returns
        -------
        dens : array
            Density estimation.
        """
        dens = self.sph_estimate(x, y, z, f=None)
        return dens


    def set_field(self, f):
        """Sets the field values for SPH field estimation.

        Parameters
        ----------
        f : array
            Field values at KDTree points.
        """
        self.f = f


    def estimate(self, x, y, z, dens=None):
        """Estimates field at points given.

        Parameters
        ----------
        x, y, z : array
            Cartesian coordinates.
        dens : array
            Density at cartesian coordinates, if None this will be computed but will
            require double the time.

        Returns
        -------
        f_est : array
            Field estimation.
        """
        return self.sph_estimate(x, y, z, f=self.f, dens=dens)


    def clean(self):
        """Reinitialises the class."""
        self.__init__()
