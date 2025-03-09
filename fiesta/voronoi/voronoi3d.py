import numpy as np
from scipy.spatial import Voronoi as scVoronoi

from .. import boundary
from .. import coords
from .. import src
from .. import utils


class Voronoi3D:


    def __init__(self):
        """Initialises Voronoi3D class"""
        self.points = None
        self.npart = None
        self.voronoi = None
        self.voronoi_volume = None
        self.voronoi_dens = None
        self.vertices = None
        self.cell = None
        self.regions = None
        self.ridge_vertices = None
        self.ridge_points = None
        self.extent = None
        self.boxsize = None
        self.buffer_length = None
        self.usebuffer = None
        self.nbuffer = None
        self.ispart = None
        self.useperiodic = None
        self.nperiodic = None


    def _extent(self):
        """"Calculates and stores minimum and maximum of
        input points in the variable extent."""
        if self.extent is None:
            xmin = np.min(self.points[:, 0])
            xmax = np.max(self.points[:, 0])
            ymin = np.min(self.points[:, 1])
            ymax = np.max(self.points[:, 1])
            zmin = np.min(self.points[:, 1])
            zmax = np.max(self.points[:, 1])
            self.extent = [xmin, xmax, ymin, ymax, zmin, zmax]

    def set_points(self, x, y, z):
        """Sets the points for voronoi cells.

        Parameters
        ----------
        x : array
            X-coordinates.
        y : array
            Y-coordinates.
        z : array
            Z-coordinates.
        """
        self.points = coords.xyz2points(x, y, z)
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
        assert self.extent[4] >= 0. and self.extent[5] <= boxsize, "Z coordinates exceed the range of the box, check or redefine boxsize."
        self.boxsize = boxsize
        self.buffer_length = buffer_length
        self.usebuffer = True
        self.useperiodic = False
        x_buffer, y_buffer, z_buffer = boundary.buffer_random_particles_3d(self.npart, self.boxsize, self.buffer_length)
        self.nbuffer = len(x_buffer)
        # redefine points to include buffer points and also define mask
        self.ispart = np.ones(self.npart + self.nbuffer)
        self.ispart[self.npart:] = 0.
        # concatenate points and buffer
        x_pnb = np.concatenate([self.points[:, 0], x_buffer])
        y_pnb = np.concatenate([self.points[:, 1], y_buffer])
        z_pnb = np.concatenate([self.points[:, 2], z_buffer])
        self.points = coords.xyz2points(x_pnb, y_pnb, z_pnb)


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
        assert self.extent[4] >= 0. and self.extent[5] <= boxsize, "Z coordinates exceed the range of the box, check or redefine boxsize."
        self.boxsize = boxsize
        self.buffer_length = buffer_length
        x_periodic, y_periodic, z_periodic = boundary.buffer_periodic_particles_3d(self.points[:, 0], self.points[:, 1], self.points[:, 2], self.boxsize, self.buffer_length)
        self.nperiodic = len(x_periodic)
        self.usebuffer = False
        self.useperiodic = True
        # redefine points to include periodic points and also define mask
        self.ispart = np.ones(self.npart + self.nperiodic)
        self.ispart[self.npart:] = 0.
        # concatenate points and buffer
        x_pnb = np.concatenate([self.points[:, 0], x_periodic])
        y_pnb = np.concatenate([self.points[:, 1], y_periodic])
        z_pnb = np.concatenate([self.points[:, 2], z_periodic])
        self.points = coords.xyz2points(x_pnb, y_pnb, z_pnb)


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


    def get_volume(self, badval=np.nan):
        """Calculates the volume of the voronoi cells.

        Parameters
        ----------
        badval : float, optional
            Bad values for the volume are set to these values.

        Returns
        -------
        volume : array
            Volume for each voronoi cell.
        """
        # find ridge information
        ridge_length = np.array([len(self.ridge_vertices[i]) for i in range(0, len(self.ridge_vertices))])
        ridge_vertices = np.array(utils.flatten_list(self.ridge_vertices))
        ridge_end = np.cumsum(ridge_length) - 1
        ridge_start = ridge_end - ridge_length + 1

        # split points and vertices to x, y and z components
        xpoints = self.points[:, 0]
        ypoints = self.points[:, 1]
        zpoints = self.points[:, 2]
        xverts = self.vertices[:, 0]
        yverts = self.vertices[:, 1]
        zverts = self.vertices[:, 2]

        # split ridge points, essentially the points that are connected by a ridge
        ridge_point1 = self.ridge_points[:, 0]
        ridge_point2 = self.ridge_points[:, 1]

        # calculate volume
        self.voronoi_volume = src.voronoi_3d_volume(xpoints, ypoints, zpoints, xverts, yverts, zverts,
            ridge_point1, ridge_point2, ridge_vertices, ridge_start, ridge_end)

        # remove and change bad values.
        cond = np.where((self.voronoi_volume == -1.) | (self.voronoi_volume == 0.))[0]
        self.voronoi_volume[cond] = badval


    def get_dens(self, mean_dens=np.nan):
        """Computes the density of voronoi cells.

        Parameters
        ----------
        mean_dens : float
            For points where a density cannot be calculated this is set to the
            mean_dens. However, if boxsize is set mean_dens is calculated automatically.
        """
        if self.voronoi_volume is None:
            self.get_volume()
        if self.boxsize is not None:
            mean_dens = self.npart/(self.boxsize**3.)
        self.voronoi_dens = np.ones(len(self.voronoi_volume))*mean_dens
        cond = np.where(np.isfinite(self.voronoi_volume) == True)[0]
        self.voronoi_dens[cond] = 1./self.voronoi_volume[cond]


    def clean(self):
        """Reinitialises the class"""
        self.__init__()
