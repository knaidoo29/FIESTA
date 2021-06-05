import numpy as np
from scipy.spatial import KDTree as scKDTree

from .. import coords

class KDTree2D:


    def __init__(self):
        """Initialises the 2D KDTree class"""
        self.points = None
        self.KD = None


    def build_tree(self, x, y):
        """Function for building the KDTree.

        Parameters
        ----------
        x : array
            X coordinates.
        y : array
            Y coordinates.
        """
        self.points = coords.xy2points(x, y)
        self.KD = scKDTree(self.points)


    def nearest(self, x, y, return_dist=False):
        """Returns the nearest index (and distance) of a point from the KDTree

        Parameters
        ----------
        x : array
            X coordinates.
        y : array
            Y coordinates.
        return_dist : bool, optional
            If True the distance to the nearest point will also be supplied.

        Return
        ------
        nind : int
            Index of nearest point.
        ndist : float, optional
            Distance to nearest point.
        """
        points = coords.xy2points(x, y)
        ndist, nind = self.KD.query(points)
        if return_dist == False:
            return nind
        else:
            return nind, ndist

    def clean(self):
        """Reinitilises the class."""
        self.__init__()
