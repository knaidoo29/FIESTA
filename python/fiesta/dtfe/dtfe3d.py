import numpy as np
from scipy.spatial import Delaunay as scDelaunay

from .. import boundary
from .. import coords
from .. import src


class Delaunay3D:


    def __init__(self):
        """Initialises Delaunay3D class"""
        self.points = None
        self.npart = None
        self.delaunay = None
        self.delaunay_simplices = None
        self.x0 = None
        self.y0 = None
        self.z0 = None
        self.f0 = None
        self.delf0 = None
        self.extent = None
        self.boxsize = None
        self.buffer_length = None
        self.usebuffer = False
        self.nbuffer = None
        self.ispart = None
        self.useperiodic = False
        self.nperiodic = None
        self.userboundary = False
        self.nuserboundary = None


    def _extent(self):
        """"Calculates and stores minimum and maximum of
        input points in the variable extent."""
        if self.extent is None:
            xmin = np.min(self.points[:self.npart, 0])
            xmax = np.max(self.points[:self.npart, 0])
            ymin = np.min(self.points[:self.npart, 1])
            ymax = np.max(self.points[:self.npart, 1])
            zmin = np.min(self.points[:self.npart, 2])
            zmax = np.max(self.points[:self.npart, 2])
            self.extent = [xmin, xmax, ymin, ymax, zmin, zmax]


    def set_points(self, x, y, z, ispart=None):
        """Sets the points for voronoi cells.

        Parameters
        ----------
        x : array
            X-coordinates.
        y : array
            Y-coordinates.
        z : array
            Z-coordinates.
        ispart : array
            Binary mask to indicate if input points already include boundary particles.
        """
        self.points = coords.xyz2points(x, y, z)
        if ispart is None:
            self.npart = len(self.points)
        else:
            self.npart = int(np.sum(ispart))
            self.ispart = ispart
            self.userboundary = True
            self.nuserboundary = len(self.points) - self.npart


    def set_buffer(self, boxsize, buffer_length):
        """Defines buffer particles to be placed around a given box.

        Parameters
        ----------
        boxsize : float
            Size of the box, assumed particles lie in the range [0., boxsize]
        buffer_length : float
            Length of the buffer region.
        """
        # check userboundary particles have not been predefined.
        assert self.userboundary == False, "User defined boundary particles already defined."
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
        z_pnb = np.concatenate([self.points[:, 2], y_buffer])
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
        # check userboundary particles have not been predefined.
        assert self.userboundary == False, "User defined boundary particles already defined."
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
        """Constructs Delaunay tesselation"""
        self.delaunay = scDelaunay(self.points)
        self.delaunay_simplices = self.delaunay.simplices


    def find_simplex(self, x, y, z):
        """Find the simplex the coordinates lie within."""
        points = coords.xyz2points(x, y, z)
        simplices = self.delaunay.find_simplex(points)
        return simplices


    def set_field(self, f, bufferval=0.):
        """Sets the field values of the input points.

        Parameters
        ----------
        f : array
            Field values.
        bufferval : float, optional
            Field values to assign boundary particles.
        """
        lenf = len(f)
        if self.userboundary == True:
            assert lenf == self.npart+self.nuserboundary, "f must be equal to input points."
        else:
            assert lenf == self.npart, "f must be equal to input points."
        x, y, z = self.points[:, 0], self.points[:, 1], self.points[:, 2]
        if self.usebuffer == True:
            f = np.concatenate([f, bufferval*np.ones(self.nbuffer)])
        elif self.useperiodic == True:
            f = np.concatenate([f, bufferval*np.ones(self.nperiodic)])
        del_vert0 = self.delaunay_simplices[:, 0]
        del_vert1 = self.delaunay_simplices[:, 1]
        del_vert2 = self.delaunay_simplices[:, 2]
        del_vert3 = self.delaunay_simplices[:, 3]
        self.x0 = x[del_vert0]
        self.y0 = y[del_vert0]
        self.z0 = z[del_vert0]
        self.f0 = f[del_vert0]
        self.delf0 = src.get_delf0_3d(x, y, z, f, del_vert0, del_vert1, del_vert2, del_vert3, len(x), len(del_vert0))


    def estimate(self, x, y, z):
        """Estimates a field from the Delaunay tesselation.

        Parameters
        ----------
        x : array
            X-coordinate for the field estimation.
        y : array
            Y-coordinate for the field estimation.
        z : array
            Z-coordinate for the field estimation.

        Returns
        -------
        f_est : array
            Estimates of the field
        """
        simplices = self.find_simplex(x, y, z)
        f_est = src.delaunay_estimate_3d(simplices, x, y, z, self.x0, self.y0, self.z0, self.f0, self.delf0, len(x), len(self.x0))
        return f_est


    def clean(self):
        self.__init__()
