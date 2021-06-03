import numpy as np
from scipy.spatial import KDTree as scKDTree


class KDTree3D:


    def __init__(self):
        """Initialises the 3D KDTree class"""
        self.points = None
        self.KD = None


    def _xyz2points(self, x, y, z):
        """Column stacks input coordinates.

        Parameters
        ----------
        x : array
            X coordinates.
        y : array
            Y coordinates.
        z : array
            Z coordinates.

        Return
        ------
        points : 2darray
            Column stacked array.
        """
        if np.isscalar(x) == True:
            points = np.array([[x, y, z]])
        else:
            points = np.column_stack((x, y, z))
        return points


    def build_tree(self, x, y, z):
        """Function for building the KDTree.

        Parameters
        ----------
        x : array
            X coordinates.
        y : array
            Y coordinates.
        z : array
            Z coordinates.
        """
        self.points = self._xyz2points(x, y, z)
        self.KD = scKDTree(self.points)


    def nearest(self, x, y, z, return_dist=False):
        """Returns the nearest index (and distance) of a point from the KDTree

        Parameters
        ----------
        x : array
            X coordinates.
        y : array
            Y coordinates.
        z : array
            Z coordinates.
        return_dist : bool, optional
            If True the distance to the nearest point will also be supplied.

        Return
        ------
        nind : int
            Index of nearest point.
        ndist : float, optional
            Distance to nearest point.
        """
        points = self._xyz2points(x, y, z)
        ndist, nind = self.KD.query(points)
        if return_dist == False:
            return nind
        else:
            return nind, ndist

    def clean(self):
        """Reinitilises the class."""
        self.__init__()
