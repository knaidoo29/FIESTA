import numpy as np

import shift

from . import points


def split_limits_by_grid(boxsize, origin, ngrid, MPI):
    """Returns the ranges along one axis for data being split between nodes.

    Parameters
    ----------
    boxsize : float
        Box length.
    origin : float
        Origin of the box.
    ngrid : int
        Number of grid divisions along the axis.
    MPI : class object
        MPIutils MPI class object.

    Returns
    -------
    limits : float list
        Range.
    """
    edges, grid = shift.cart.grid1D(boxsize, ngrid, origin=origin)
    s1, s2 = MPI.split(len(grid))
    limits = [edges[s1[MPI.rank]], edges[s2[MPI.rank]]]
    return limits


def mpi_find_range_2D(fnames, freader, MPI):
    """Find ranges from each file and returns the ranges to node 0.

    Parameters
    ----------
    fnames : list str
        List of fname inputs for freader function.
    freader : func
        File reader function which returns ndarray with first two columns
        corresponding to x and y coordinate axis.
    MPI : class object
        MPIutils MPI class object.

    Returns
    -------
    ranges : ndarray
        Ranges with columns: xmin, xmax, ymin and ymax.
    """
    _fnames = MPI.split_array(fnames)
    xmins, xmaxs, ymins, ymaxs = [], [], [], []
    for i in range(0, len(_fnames)):
        data = freader(_fnames[i])
        xmins.append(np.min(data[:,0]))
        ymins.append(np.min(data[:,1]))
        xmaxs.append(np.max(data[:,0]))
        ymaxs.append(np.max(data[:,1]))
    xmins, xmaxs = np.array(xmins), np.array(xmaxs)
    ymins, ymaxs = np.array(ymins), np.array(ymaxs)
    ranges = points.coord2points([xmins, xmaxs, ymins, ymaxs])
    ranges = MPI.collect(ranges, outlist=True)
    if MPI.rank == 0:
        ranges = np.vstack(ranges)
    return ranges


def mpi_find_range_3D(fnames, freader, MPI):
    """Find ranges from each file and returns the ranges to node 0.

    Parameters
    ----------
    fnames : list str
        List of fname inputs for freader function.
    freader : func
        File reader function which returns ndarray with first two columns
        corresponding to x and y coordinate axis.
    MPI : class object
        MPIutils MPI class object.

    Returns
    -------
    ranges : ndarray
        Ranges with columns: xmin, xmax, ymin, ymax, zmin and zmax.
    """
    _fnames = MPI.split_array(fnames)
    xmins, xmaxs, ymins, ymaxs, zmins, zmaxs = [], [], [], [], [], []
    for i in range(0, len(_fnames)):
        data = freader(_fnames[i])
        xmins.append(np.min(data[:,0]))
        ymins.append(np.min(data[:,1]))
        zmins.append(np.min(data[:,2]))
        xmaxs.append(np.max(data[:,0]))
        ymaxs.append(np.max(data[:,1]))
        zmaxs.append(np.max(data[:,2]))
    xmins, xmaxs = np.array(xmins), np.array(xmaxs)
    ymins, ymaxs = np.array(ymins), np.array(ymaxs)
    zmins, zmaxs = np.array(zmins), np.array(zmaxs)
    ranges = points.coord2points([xmins, xmaxs, ymins, ymaxs, zmins, zmaxs])
    ranges = MPI.collect(ranges, outlist=True)
    if MPI.rank == 0:
        ranges = np.vstack(ranges)
    return ranges


def mpi_open_2D(fnames, freader, ranges, limits, MPI):
    """Find ranges from each file and returns the ranges to node 0.

    Parameters
    ----------
    fnames : list str
        List of fname inputs for freader function.
    freader : func
        File reader function which returns ndarray with first two columns
        corresponding to x and y coordinate axis.
    ranges : ndarray
        Ranges with columns: xmin, xmax, ymin and ymax.
    limits : list
        Ranges for each coordinate axis.
    MPI : class object
        MPIutils MPI class object.

    Returns
    -------
    ranges : ndarray
        Ranges with columns: xmin, xmax, ymin and ymax.
    """
    xmins, xmaxs, ymins, ymaxs = ranges[:,0], ranges[:, 1], ranges[:,2], ranges[:, 3]
    xmin, xmax, ymin, ymax = limits[0], limits[1], limits[2], limits[3]
    cond = np.where(((xmins >= xmin) & (ymins >= ymin) & (xmins < xmax) & (ymins < ymax)) |
                    ((xmaxs <= xmax) & (ymaxs <= ymax) & (xmaxs > xmin) & (ymaxs > ymin)))[0]
    datas = []
    for i in range(0, len(cond)):
        _data = freader(fnames[cond[i]])
        cond1 = np.where((_data[:,0] >= xmin) & (_data[:,0] < xmax) &
                         (_data[:,1] >= ymin) & (_data[:,1] < ymax))[0]
        datas.append(_data[cond1])
    datas = np.vstack(datas)
    return datas


def mpi_open_3D(fnames, freader, ranges, limits, MPI):
    """Find ranges from each file and returns the ranges to node 0.

    Parameters
    ----------
    fnames : list str
        List of fname inputs for freader function.
    freader : func
        File reader function which returns ndarray with first two columns
        corresponding to x and y coordinate axis.
    ranges : ndarray
        Ranges with columns: xmin, xmax, ymin and ymax.
    limits : list
        Ranges for each coordinate axis.
    MPI : class object
        MPIutils MPI class object.

    Returns
    -------
    ranges : ndarray
        Ranges with columns: xmin, xmax, ymin, ymax, zmin and zmax.
    """
    xmins, xmaxs, ymins, ymaxs, zmins, zmaxs = ranges[:,0], ranges[:,1], ranges[:,2], ranges[:,3], ranges[:,4], ranges[:,5]
    xmin, xmax, ymin, ymax, zmin, zmax = limits[0], limits[1], limits[2], limits[3], limits[4], limits[5]
    cond = np.where(((xmins >= xmin) & (ymins >= ymin) & (zmins >= zmin) &
                     (xmins < xmax) & (ymins < ymax) & (zmins < zmax)) |
                    ((xmaxs <= xmax) & (ymaxs <= ymax) & (zmaxs <= zmax) &
                     (xmaxs > xmin) & (ymaxs > ymin) & (zmaxs > zmin)))[0]
    datas = []
    for i in range(0, len(cond)):
        _data = freader(fnames[cond[i]])
        cond1 = np.where((_data[:,0] >= xmin) & (_data[:,0] < xmax) &
                         (_data[:,1] >= ymin) & (_data[:,1] < ymax) &
                         (_data[:,2] >= zmin) & (_data[:,2] < zmax))[0]
        datas.append(_data[cond1])
    datas = np.vstack(datas)
    return datas


def check_coords_at_MPI_0(x, MPI):
    """Check if coordinates are only inputed at the zeroth node.

    Parameter
    ---------
    x : array
        X-axis coordinates.
    MPI : class object
        MPIutils MPI class object.

    Returns
    -------
    checkatzero : bool
        True if points are only at the zeroth node.
    """
    if x is not None:
        nocoord = True
    else:
        nocoord = False
    nocoords = MPI.collect(nocoord, outlist=True)
    if MPI.rank == 0:
        if nocoords[0] == True:
            checkatzero = True
        for i in range(1, len(nocoords)):
            if nocoords[i] == True:
                checkatzero = False
        MPI.send(checkatzero, tag=11)
    else:
        checkatzero = MPI.recv(0, tag=11)
    return checkatzero



class MPI_SortByX:


    def __init__(self, MPI):
        self.MPI = MPI
        self.boxsize = None
        self.origin = None
        self.ngrid = None
        self.ngrid_rank = None
        self.data = None
        self.limits = None
        self.all_limits = None


    def settings(self, boxsize, ngrid, origin=0.):
        self.boxsize = boxsize
        self.ngrid = ngrid
        self.origin = origin


    def input(self, data):
        self.data = data


    def limits4grid(self):
        xedges, xgrid = shift.cart.mpi_grid1D(self.boxsize, self.ngrid, self.MPI, origin=self.origin)
        self.limits = [xedges[0], xedges[-1]]
        self.ngrid_rank = len(xgrid)
        all_limits = self.MPI.collect(self.limits, outlist=True)
        if self.MPI.rank == 0:
            all_limits = np.vstack(all_limits)
            self.MPI.send(all_limits, tag=11)
        else:
            all_limits = self.MPI.recv(0, tag=11)
        self.all_limits = all_limits


    def distribute(self):
        for i in range(0, self.MPI.size):
            if self.data is not None:
                cond = np.where((self.data[:,0] >= self.all_limits[i,0]) &
                                (self.data[:,0] < self.all_limits[i,1]))[0]
                _data = self.data[cond]
                _hasdata = True
            else:
                _data = None
                _hasdata = False
            _data = self.MPI.collect(_data, outlist=True)
            _hasdata = self.MPI.collect(_hasdata, outlist=True)
            if self.MPI.rank == 0:
                _hasdata = np.array(_hasdata).flatten()
                cond = np.where(_hasdata == True)[0]
                if len(cond) != 0:
                    _data = np.vstack([_data[c] for c in cond])
                else:
                    _data = None
                if i == 0:
                    _data4limit = np.copy(_data)
                else:
                    self.MPI.send(_data, to_rank=i, tag=10+i)
            else:
                if i != 0 and i == self.MPI.rank:
                    _data4limit = self.MPI.recv(0, tag=10+i)
            self.MPI.wait()
        return _data4limit


    def clean(self):
        self.__init__()
