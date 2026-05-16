import numpy as np
from scipy.spatial import KDTree as scKDTree

from .. import coords

from typing import Optional, Tuple


class KDTree2D:

    
    def __init__(self) -> None:
        """
        Initialises the 2D KDTree class.
        """
        self.points = None
        self.KD = None


    def build_tree(
        self,
        x: np.ndarray,
        y: np.ndarray,
        boxsize: Optional[float] = None,
    ) -> None:
        """
        Function for building the KDTree.

        Parameters
        ----------
        x : array
            X coordinates.
        y : array
            Y coordinates.
        boxsize : float, optional
            Periodic boundary boxsize.
        """
        self.points = coords.xy2points(x, y)
        self.KD = scKDTree(self.points, boxsize=boxsize)


    def nearest(
        self, x: np.ndarray, y: np.ndarray, k: int = 1, return_dist: bool = False
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Returns the nearest index (and distance) of a point from the KDTree.

        Parameters
        ----------
        x : array
            X coordinates.
        y : array
            Y coordinates.
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
        points = coords.xy2points(x, y)
        ndist, nind = self.KD.query(points, k=k)
        if return_dist == False:
            return nind
        else:
            return nind, ndist


    def find_points_in_r(self, x: np.ndarray, y: np.ndarray, r: float) -> np.ndarray:
        """
        Returns the nearest index (and distance) of a point from the KDTree.

        Parameters
        ----------
        x : array
            X coordinates.
        y : array
            Y coordinates.
        r : float
            Radial distance to find points in KDTree from input coordinates.

        Return
        ------
        ind : int
            Index of points in the KDTree that are within a distance r of the input
            coordinates..
        """
        points = coords.xy2points(x, y)
        ind = self.KD.query_ball_point(points, r)
        return ind

    
    def clean(self) -> None:  # pragma: no cover
        """
        Reinitilises the class.
        """
        self.__init__()
