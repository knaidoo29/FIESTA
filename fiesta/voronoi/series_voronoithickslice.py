import numpy as np

from .. import boundary
from .. import coords

from . import voronoithickslice

class SeriesVoronoiThickSlice:


    def __init__(self):
        """Initialises VoronoiThickSlice class"""
        self.points = None
        self.npart = None
        self.extent = None
        self.boxsize = None
        self.thickness = None
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
            self.extent = [xmin, xmax, ymin, ymax]


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


    def set_buffer(self, boxsize, thickness, buffer_length):
        """Defines buffer particles to be placed around a given box.

        Parameters
        ----------
        boxsize : float
            Size of the box, assumed particles lie in the range [0., boxsize]
        thickness : float
            Thickness of the slice.
        buffer_length : float
            Length of the buffer region.
        """
        # check boxsize is consistent with particles.
        self._extent()
        assert self.extent[0] >= 0. and self.extent[1] <= boxsize, "X coordinates exceed the range of the box, check or redefine boxsize."
        assert self.extent[2] >= 0. and self.extent[3] <= boxsize, "Y coordinates exceed the range of the box, check or redefine boxsize."
        self.boxsize = boxsize
        self.buffer_length = buffer_length
        self.thickness = thickness
        self.usebuffer = True
        self.useperiodic = False
        x_buffer, y_buffer, z_buffer = boundary.buffer_random_particles_3d_slice(self.npart, self.boxsize, self.thickness, self.buffer_length)
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
        self.boxsize = boxsize
        self.buffer_length = buffer_length
        x_periodic, y_periodic, z_periodic = boundary.buffer_periodic_particles_3d_slice(self.points[:, 0], self.points[:, 1], self.points[:, 2], self.boxsize, self.buffer_length)
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


    def get_volume(self, split=2, badval=np.nan):
        """Calculates the volume of the voronoi cells.

        Parameters
        ----------
        split : int, optional
            Split calculation along one axis, default = 1, meaning no splitting.
        badval : float, optional
            Bad values for the volume are set to these values.

        Returns
        -------
        volume : array
            Volume for each voronoi cell.
        """
        x_data = self.points[:, 0]
        y_data = self.points[:, 1]
        z_data = self.points[:, 2]
        split_edges = np.linspace(0., self.boxsize, split+1)
        split_centers = 0.5*(split_edges[:-1] + split_edges[1:])
        dsplit = split_edges[1] - split_edges[0]
        x_split = np.copy(split_centers)
        y_split = np.copy(split_centers)
        dx_est = 0.5*dsplit
        if self.buffer_length is None:
            dx_data = dx_est
        else:
            dx_data = dx_est + self.buffer_length
        volume = np.zeros(len(x_data))
        VTS = voronoithickslice.VoronoiThickSlice()
        for i in range(0, len(x_split)):
            xcond_data = np.where((x_data >= x_split[i]-dx_data) & (x_data < x_split[i]+dx_data))[0]
            for j in range(0, len(y_split)):
                ycond_data = np.where((y_data[xcond_data] >= y_split[j]-dx_data) & (y_data[xcond_data] < y_split[j]+dx_data))[0]
                cond_data = xcond_data[ycond_data]
                x_d_lim = x_data[cond_data]
                y_d_lim = y_data[cond_data]
                z_d_lim = z_data[cond_data]
                cond_est = np.where((x_d_lim >= x_split[i]-dx_est) & (x_d_lim < x_split[i]+dx_est) &
                                    (y_d_lim >= y_split[j]-dx_est) & (y_d_lim < y_split[j]+dx_est))[0]
                VTS.set_points(x_d_lim, y_d_lim, z_d_lim)
                VTS.construct()
                volume_d_lim = VTS.get_volume(badval=badval)
                volume[cond_data[cond_est]] = volume_d_lim[cond_est]
                VTS.clean()
        if self.ispart is not None:
            cond = np.where(self.ispart == 1.)[0]
            volume = volume[cond]
        self.volume = volume
        return volume


    def clean(self):
        """Reinitialises the class"""
        self.__init__()
