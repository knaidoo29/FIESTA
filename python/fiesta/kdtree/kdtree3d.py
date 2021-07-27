import numpy as np
from scipy.spatial import KDTree as scKDTree

from .. import coords


class KDTree3D:


    def __init__(self):
        """Initialises the 3D KDTree class"""
        self.points = None
        self.KD = None
        self.usepara = None
        self.ncpu = None


    def build_tree(self, x, y, z, boxsize=None, usepara=False, ncpu=4):
        """Function for building the KDTree.

        Parameters
        ----------
        x : array
            X coordinates.
        y : array
            Y coordinates.
        z : array
            Z coordinates.
        boxsize : float, optional
            Periodic boundary boxsize.
        """
        self.usepara = usepara
        self.ncpu = ncpu
        self.points = coords.xyz2points(x, y, z)
        self.KD = scKDTree(self.points, boxsize=boxsize)


    def nearest(self, x, y, z, k=1, return_dist=False):
        """Returns the nearest index (and distance) of a point from the KDTree

        Parameters
        ----------
        x : array
            X coordinates.
        y : array
            Y coordinates.
        z : array
            Z coordinates.
        k : int
            Number of nearest points.
        return_dist : bool, optional
            If True the distance to the nearest point will also be supplied.

        Return
        ------
        nind : int
            Index of nearest point.
        ndist : float, optional
            Distance to nearest point.
        """
        points = coords.xyz2points(x, y, z)
        if self.usepara == False:
            ndist, nind = self.KD.query(points, k=k)
        else:
            ndist, nind = self.KD.query(points, k=k, workers=self.ncpu)
        if return_dist == False:
            return nind
        else:
            return nind, ndist


    def find_points_in_r(self, x, y, z, r):
        """Returns the nearest index (and distance) of a point from the KDTree.

        Parameters
        ----------
        x : array
            X coordinates.
        y : array
            Y coordinates.
        z : array
            Z coordinates.
        r : float
            Radial distance to find points in KDTree from input coordinates.

        Return
        ------
        ind : int
            Index of points in the KDTree that are within a distance r of the input
            coordinates..
        """
        points = coords.xyz2points(x, y, z)
        ind = self.KD.query_ball_point(points, r)
        return ind


    def clean(self):
        """Reinitilises the class."""
        self.__init__()
