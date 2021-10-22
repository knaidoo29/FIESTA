import numpy as np
from scipy.spatial import Voronoi as scVoronoi

from .. import boundary
from .. import coords
from .. import src
from .. import utils


class Voronoi2D:


    def __init__(self):
        """Initialises Voronoi2D class"""
        self.points = None
        self.npart = None
        self.voronoi = None
        self.vertices = None
        self.cell = None
        self.regions = None
        self.ridge_vertices = None
        self.ridge_points = None
        self.extent = None
        self.boxsize = None
        self.buffer_length = None
        self.usebuffer = False
        self.nbuffer = None
        self.ispart = None
        self.useperiodic = False
        self.nperiodic = None


    def _extent(self):
        """"Calculates and stores minimum and maximum of
        input points in the variable extent."""
        if self.extent is None:
            xmin = np.min(self.points[:self.npart, 0])
            xmax = np.max(self.points[:self.npart, 0])
            ymin = np.min(self.points[:self.npart, 1])
            ymax = np.max(self.points[:self.npart, 1])
            self.extent = [xmin, xmax, ymin, ymax]


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
        self.npart = len(self.points)


    def set_buffer(self, boxsize, buffer_length):
        """Defines buffer particles to be placed around a given box.

        Parameters
        ----------
        boxsize : float
            Size of the box, assumed particles lie in the range [0., boxsize]
        buffer_length : float
            Length of the buffer region.
        """
        # check boxsize is consistent with particles.
        self._extent()
        assert self.extent[0] >= 0. and self.extent[1] <= boxsize, "X coordinates exceed the range of the box, check or redefine boxsize."
        assert self.extent[2] >= 0. and self.extent[3] <= boxsize, "Y coordinates exceed the range of the box, check or redefine boxsize."
        self.boxsize = boxsize
        self.buffer_length = buffer_length
        self.usebuffer = True
        self.useperiodic = False
        x_buffer, y_buffer = boundary.buffer_random_particles_2d(self.npart, self.boxsize, self.buffer_length)
        self.nbuffer = len(x_buffer)
        # redefine points to include buffer points and also define mask
        self.ispart = np.ones(self.npart + self.nbuffer)
        self.ispart[self.npart:] = 0.
        # concatenate points and buffer
        x_pnb = np.concatenate([self.points[:, 0], x_buffer])
        y_pnb = np.concatenate([self.points[:, 1], y_buffer])
        self.points = coords.xy2points(x_pnb, y_pnb)


    def set_periodic(self, boxsize, buffer_length):
        """Defines periodic particles to be placed around a given box.

        Parameters
        ----------
        boxsize : float
            Size of the box, assumed particles lie in the range [0., boxsize]
        buffer_length : float
            Length of the buffer region.
        """
        # check boxsize is consistent with particles.
        self._extent()
        assert self.extent[0] >= 0. and self.extent[1] <= boxsize, "X coordinates exceed the range of the box, check or redefine boxsize."
        assert self.extent[2] >= 0. and self.extent[3] <= boxsize, "Y coordinates exceed the range of the box, check or redefine boxsize."
        self.boxsize = boxsize
        self.buffer_length = buffer_length
        x_periodic, y_periodic = boundary.buffer_periodic_particles_2d(self.points[:, 0], self.points[:, 1], self.boxsize, self.buffer_length)
        self.nperiodic = len(x_periodic)
        self.usebuffer = False
        self.useperiodic = True
        # redefine points to include periodic points and also define mask
        self.ispart = np.ones(self.npart + self.nperiodic)
        self.ispart[self.npart:] = 0.
        # concatenate points and buffer
        x_pnb = np.concatenate([self.points[:, 0], x_periodic])
        y_pnb = np.concatenate([self.points[:, 1], y_periodic])
        self.points = coords.xy2points(x_pnb, y_pnb)


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
        split : int, optional
            Split calculation along one axis, default = 1, meaning no splitting.

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
        area = src.voronoi_2d_area(xpoints=xpoints, ypoints=ypoints, xverts=xverts, yverts=yverts,
                                   ridge_point1=ridge_point1, ridge_point2=ridge_point2,
                                   ridge_vertices=ridge_vertices, ridge_start=ridge_start,
                                   ridge_end=ridge_end, npoints=len(xpoints), nridge=len(ridge_point1),
                                   nvertices=len(xverts), nridge_vertices=len(ridge_vertices))

        # remove and change bad values.
        cond = np.where((area == -1.) | (area == 0.))[0]
        area[cond] = badval
        if self.ispart is None:
            cond = np.where(self.ispart == 1.)[0]
            area = area[cond]
        self.area = area
        return area


    def clean(self):
        """Reinitialises the class"""
        self.__init__()
