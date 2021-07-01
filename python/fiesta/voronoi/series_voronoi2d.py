import numpy as np

from .. import boundary
from .. import coords

from . import voronoi2d


class SeriesVoronoi2D:


    def __init__(self):
        """Initialises Voronoi2D class"""
        self.points = None
        self.npart = None
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


    def get_area(self, split=2, badval=np.nan):
        """Calculates the area of the voronoi cells.

        Parameters
        ----------
        split : int, optional
            Split calculation along one axis, default = 1, meaning no splitting.
        badval : float, optional
            Bad values for the area are set to these values.

        Returns
        -------
        area : array
            Area for each voronoi cell.
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
        area = np.zeros(len(x_data))
        V2D = voronoi2d.Voronoi2D()
        for i in range(0, len(x_split)):
            xcond_data = np.where((x_data >= x_split[i]-dx_data) & (x_data < x_split[i]+dx_data))[0]
            for j in range(0, len(y_split)):
                ycond_data = np.where((y_data[xcond_data] >= y_split[j]-dx_data) & (y_data[xcond_data] < y_split[j]+dx_data))[0]
                cond_data = xcond_data[ycond_data]
                x_d_lim = x_data[cond_data]
                y_d_lim = y_data[cond_data]
                cond_est = np.where((x_d_lim >= x_split[i]-dx_est) & (x_d_lim < x_split[i]+dx_est) &
                                    (y_d_lim >= y_split[j]-dx_est) & (y_d_lim < y_split[j]+dx_est))[0]
                V2D.set_points(x_d_lim, y_d_lim)
                V2D.construct()
                area_d_lim = V2D.get_area(badval=badval)
                area[cond_data[cond_est]] = area_d_lim[cond_est]
                V2D.clean()
        if self.ispart is not None:
            cond = np.where(self.ispart == 1.)[0]
            area = area[cond]
        self.area = area
        return area


    def clean(self):
        """Reinitialises the class"""
        self.__init__()
