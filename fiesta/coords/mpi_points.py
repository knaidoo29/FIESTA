import numpy as np

from . import points


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
