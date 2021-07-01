import numpy as np

from .. import boundary
from .. import coords

from . import dtfe3d


class SeriesDelaunay3D:


    def __init__(self):
        """Initialises Delaunay3D class"""
        self.points = None
        self.npart = None
        self.f = None
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
            zmin = np.min(self.points[:self.npart, 2])
            zmax = np.max(self.points[:self.npart, 2])
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
        assert lenf == self.npart, "f must be equal to input points."
        if self.usebuffer == True:
            f = np.concatenate([f, bufferval*np.ones(self.nbuffer)])
        elif self.useperiodic == True:
            f = np.concatenate([f, bufferval*np.ones(self.nperiodic)])
        self.f = f


    def estimate(self, x_est, y_est, z_est, split=2):
        """Estimates a field from the Delaunay tesselation.

        Parameters
        ----------
        x_est : array
            X-coordinate for the field estimation.
        y_est : array
            Y-coordinate for the field estimation.
        z_est : array
            Z-coordinate for the field estimation.
        split : int
            Number of divisions to perform the DTFE.

        Returns
        -------
        f_est : array
            Estimates of the field
        """
        x_data = self.points[:, 0]
        y_data = self.points[:, 1]
        z_data = self.points[:, 2]
        split_edges = np.linspace(0., self.boxsize, split+1)
        split_centers = 0.5*(split_edges[:-1] + split_edges[1:])
        dsplit = split_edges[1] - split_edges[0]
        x_split = split_centers
        y_split = split_centers
        z_split = split_centers
        dx_est = 0.5*dsplit
        if self.buffer_length is None:
            dx_data = dx_est
        else:
            dx_data = dx_est + self.buffer_length
        f_est = np.zeros(len(x_est))
        D3D = dtfe3d.Delaunay3D()
        for i in range(0, len(x_split)):
            xcond_data = np.where((x_data >= x_split[i]-dx_data) & (x_data < x_split[i]+dx_data))[0]
            for j in range(0, len(y_split)):
                ycond_data = np.where((y_data[xcond_data] >= y_split[j]-dx_data) & (y_data[xcond_data] < y_split[j]+dx_data))[0]
                for k in range(0, len(z_split)):
                    zcond_data = np.where((z_data[xcond_data[ycond_data]] >= z_split[k]-dx_data) &
                                          (z_data[xcond_data[ycond_data]] < z_split[k]+dx_data))[0]
                    cond_data = xcond_data[ycond_data[zcond_data]]
                    x_d_lim = x_data[cond_data]
                    y_d_lim = y_data[cond_data]
                    z_d_lim = z_data[cond_data]
                    f_d_lim = self.f[cond_data]
                    cond_est = np.where((x_est >= x_split[i]-dx_est) & (x_est < x_split[i]+dx_est) &
                                        (y_est >= y_split[i]-dx_est) & (y_est < y_split[i]+dx_est) &
                                        (z_est >= z_split[i]-dx_est) & (z_est < z_split[i]+dx_est))[0]
                    if len(cond_est) > 0:
                        x_e_lim = x_est[cond_est]
                        y_e_lim = y_est[cond_est]
                        z_e_lim = z_est[cond_est]
                        D3D.set_points(x_d_lim, y_d_lim, z_d_lim)
                        D3D.construct()
                        D3D.set_field(f_d_lim)
                        f_e_lim = D3D.estimate(x_e_lim, y_e_lim, z_e_lim)
                        f_est[cond_est] = f_e_lim
                        D3D.clean()
        return f_est


    def clean(self):
        """Reinitialises the class"""
        self.__init__()
