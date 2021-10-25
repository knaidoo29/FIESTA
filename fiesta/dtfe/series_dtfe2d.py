import numpy as np

from .. import boundary
from .. import coords

from . import dtfe2d


class SeriesDelaunay2D:


    def __init__(self):
        """Initialises SeriesDelaunay2D class"""
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
            self.extent = [xmin, xmax, ymin, ymax]


    def set_points(self, x, y, ispart=None):
        """Sets the points for voronoi cells.

        Parameters
        ----------
        x : array
            X-coordinates.
        y : array
            Y-coordinates.
        ispart : array
            Binary mask to indicate if input points already include boundary particles.
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


    def estimate(self, x_est, y_est, split=2):
        """Estimates a field from the Delaunay tesselation in series, split^2 steps.

        Parameters
        ----------
        x_est : array
            X-coordinate for the field estimation.
        y_est : array
            Y-coordinate for the field estimation.
        split : int
            Number of divisions to perform the DTFE.

        Returns
        -------
        f_est : array
            Estimates of the field
        """
        x_data = self.points[:, 0]
        y_data = self.points[:, 1]
        split_edges = np.linspace(0., self.boxsize, split+1)
        split_centers = 0.5*(split_edges[:-1] + split_edges[1:])
        dsplit = split_edges[1] - split_edges[0]
        x_split = split_centers
        y_split = split_centers
        dx_est = 0.5*dsplit
        if self.buffer_length is None:
            dx_data = dx_est
        else:
            dx_data = dx_est + self.buffer_length
        f_est = np.zeros(len(x_est))
        D2D = dtfe2d.Delaunay2D()
        for i in range(0, len(x_split)):
            xcond_data = np.where((x_data >= x_split[i]-dx_data) & (x_data < x_split[i]+dx_data))[0]
            for j in range(0, len(y_split)):
                ycond_data = np.where((y_data[xcond_data] >= y_split[j]-dx_data) & (y_data[xcond_data] < y_split[j]+dx_data))[0]
                cond_data = xcond_data[ycond_data]
                x_d_lim = x_data[cond_data]
                y_d_lim = y_data[cond_data]
                f_d_lim = self.f[cond_data]
                cond_est = np.where((x_est >= x_split[i]-dx_est) & (x_est < x_split[i]+dx_est) &
                                    (y_est >= y_split[i]-dx_est) & (y_est < y_split[i]+dx_est))[0]
                if len(cond_est) > 0:
                    x_e_lim = x_est[cond_est]
                    y_e_lim = y_est[cond_est]
                    D2D.set_points(x_d_lim, y_d_lim)
                    D2D.construct()
                    D2D.set_field(f_d_lim)
                    f_e_lim = D2D.estimate(x_e_lim, y_e_lim)
                    f_est[cond_est] = f_e_lim
                    D2D.clean()
        return f_est


    def clean(self):
        """Reinitialises the class"""
        self.__init__()
