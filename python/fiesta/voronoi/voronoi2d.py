import numpy as np
from scipy.spatial import Voronoi as scVoronoi

from .. import coords
from .. import src
from .. import utils


class Voronoi2D:


    def __init__(self):
        """Initialises Voronoi2D class"""
        self.points = None
        self.voronoi = None
        self.vertices = None
        self.cell = None
        self.regions = None
        self.ridge_vertices = None
        self.ridge_points = None


    def set_points(self, x, y):
        """Sets the points for voronoi cells.

        Parameters
        ----------
        x : array
            X-coordinates.
        y : array
            Y-coordinates.
        """
        self.points = coords.xy2points(x, y)

    # self.buffer periodic or random or mask.


    def construct(self):
        """Constructs the voronoi tesselation from the input points"""
        self.voronoi = scVoronoi(self.points)
        # voronoi vertices
        self.vertices = self.voronoi.vertices
        # voronoi cells which coorespond to the index of the input points
        self.cell = self.voronoi.point_region
        # voronoi regions/cells
        self.regions = self.voronoi.regions
        # voronoi ridge vertices
        self.ridge_vertices = self.voronoi.ridge_vertices
        # voronoi ridge points, i.e. points connecting each face
        self.ridge_points = self.voronoi.ridge_points


    def get_area(self, badval=np.nan):
        """Calculates the area of the voronoi cells.

        Parameters
        ----------
        badval : float, optional
            Bad values for the area are set to these values.

        Returns
        -------
        area : array
            Area for each voronoi cell.
        """
        # find ridge information
        ridge_length = np.array([len(self.ridge_vertices[i]) for i in range(0, len(self.ridge_vertices))])
        ridge_vertices = np.array(utils.flatten_list(self.ridge_vertices))
        ridge_end = np.cumsum(ridge_length) - 1
        ridge_start = ridge_end - ridge_length + 1

        # split points and vertices to x and y components
        xpoints = self.points[:, 0]
        ypoints = self.points[:, 1]
        xverts = self.vertices[:, 0]
        yverts = self.vertices[:, 1]

        # split ridge points, essentially the points that are connected by a ridge
        ridge_point1 = self.ridge_points[:, 0]
        ridge_point2 = self.ridge_points[:, 1]

        # calculate area
        area = src.voronoi_2d_area(xpoints, ypoints, xverts, yverts, ridge_point1, ridge_point2,
                                   ridge_vertices, ridge_start, ridge_end, len(xpoints),
                                   len(ridge_point1), len(xverts), len(ridge_vertices))

        # remove and change bad values.
        cond = np.where((area == -1.) | (area == 0.))[0]
        area[cond] = badval
        self.area = area
        return area


    def clean(self):
        """Reinitialises the class"""
        self.__init__()
