import numpy as np

import shift

from . import points

from typing import List, Union


def split_limits_by_grid(
    boxsize: float, origin: float, ngrid: int, MPI: object
) -> List[float]:
    """
    Returns the ranges along one axis for data being split between nodes.

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


def mpi_find_range_2D(fnames: str, freader: object, MPI: object) -> np.ndarray:
    """
    Find ranges from each file and returns the ranges to node 0.

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
        xmins.append(np.min(data[:, 0]))
        ymins.append(np.min(data[:, 1]))
        xmaxs.append(np.max(data[:, 0]))
        ymaxs.append(np.max(data[:, 1]))
    xmins, xmaxs = np.array(xmins), np.array(xmaxs)
    ymins, ymaxs = np.array(ymins), np.array(ymaxs)
    ranges = points.coord2points([xmins, xmaxs, ymins, ymaxs])
    ranges = MPI.gather(ranges)
    if MPI.rank == 0:
        ranges = np.vstack(ranges)
    return ranges


def mpi_find_range_3D(fnames: str, freader: object, MPI: object) -> np.ndarray:
    """
    Find ranges from each file and returns the ranges to node 0.

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
        xmins.append(np.min(data[:, 0]))
        ymins.append(np.min(data[:, 1]))
        zmins.append(np.min(data[:, 2]))
        xmaxs.append(np.max(data[:, 0]))
        ymaxs.append(np.max(data[:, 1]))
        zmaxs.append(np.max(data[:, 2]))
    xmins, xmaxs = np.array(xmins), np.array(xmaxs)
    ymins, ymaxs = np.array(ymins), np.array(ymaxs)
    zmins, zmaxs = np.array(zmins), np.array(zmaxs)
    ranges = points.coord2points([xmins, xmaxs, ymins, ymaxs, zmins, zmaxs])
    ranges = MPI.gather(ranges)
    if MPI.rank == 0:
        ranges = np.vstack(ranges)
    return ranges


def mpi_open_2D(
    fnames: str, freader: object, ranges: np.ndarray, limits: List[float]
) -> np.ndarray:
    """
    Find ranges from each file and returns the ranges to node 0.

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

    Returns
    -------
    datas : ndarray
        Data within range.
    """
    xmins, xmaxs, ymins, ymaxs = ranges[:, 0], ranges[:, 1], ranges[:, 2], ranges[:, 3]
    xmin, xmax, ymin, ymax = limits[0], limits[1], limits[2], limits[3]
    cond = np.where(
        ((xmins >= xmin) & (ymins >= ymin) & (xmins < xmax) & (ymins < ymax))
        | ((xmaxs <= xmax) & (ymaxs <= ymax) & (xmaxs > xmin) & (ymaxs > ymin))
    )[0]
    datas = []
    for i in range(0, len(cond)):
        _data = freader(fnames[cond[i]])
        cond1 = np.where(
            (_data[:, 0] >= xmin)
            & (_data[:, 0] < xmax)
            & (_data[:, 1] >= ymin)
            & (_data[:, 1] < ymax)
        )[0]
        datas.append(_data[cond1])
    datas = np.vstack(datas)
    return datas


def mpi_open_3D(
    fnames: str, freader: object, ranges: np.ndarray, limits: List[float]
) -> np.ndarray:
    """
    Find ranges from each file and returns the ranges to node 0.

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
    datas : ndarray
        Data within range.
    """
    xmins, xmaxs, ymins, ymaxs, zmins, zmaxs = (
        ranges[:, 0],
        ranges[:, 1],
        ranges[:, 2],
        ranges[:, 3],
        ranges[:, 4],
        ranges[:, 5],
    )
    xmin, xmax, ymin, ymax, zmin, zmax = (
        limits[0],
        limits[1],
        limits[2],
        limits[3],
        limits[4],
        limits[5],
    )
    cond = np.where(
        (
            (xmins >= xmin)
            & (ymins >= ymin)
            & (zmins >= zmin)
            & (xmins < xmax)
            & (ymins < ymax)
            & (zmins < zmax)
        )
        | (
            (xmaxs <= xmax)
            & (ymaxs <= ymax)
            & (zmaxs <= zmax)
            & (xmaxs > xmin)
            & (ymaxs > ymin)
            & (zmaxs > zmin)
        )
    )[0]
    datas = []
    for i in range(0, len(cond)):
        _data = freader(fnames[cond[i]])
        cond1 = np.where(
            (_data[:, 0] >= xmin)
            & (_data[:, 0] < xmax)
            & (_data[:, 1] >= ymin)
            & (_data[:, 1] < ymax)
            & (_data[:, 2] >= zmin)
            & (_data[:, 2] < zmax)
        )[0]
        datas.append(_data[cond1])
    datas = np.vstack(datas)
    return datas


def check_coords_at_MPI_0(x: np.ndarray, MPI: object) -> bool:
    """
    Check if coordinates are only inputed at the zeroth node.

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
    nocoords = MPI.gather(nocoord)
    if MPI.rank == 0:
        if nocoords[0] == True:
            checkatzero = True
        for i in range(1, len(nocoords)):
            if nocoords[i] == True:
                checkatzero = False
    else:
        checkatzero = False
    checkatzero = MPI.bcast(checkatzero)
    return checkatzero


class MPI_SortByX:


    def __init__(self, MPI: object) -> None:
        """
        Initialise MPI point class.

        Parameters
        ----------
        MPI : object
            shift.mpiutils mpi object.
        """
        self.MPI = MPI
        self.boxsize = None
        self.origin = None
        self.ngrid = None
        self.ngrid_rank = None
        self.data = None
        self.limits = None
        self.all_limits = None
        self.buffer_length = 0.0


    def settings(
        self,
        boxsize: Union[float, List[float]],
        ngrid: Union[int, List[int]],
        origin: Union[float, List[float]] = 0.0,
        buffer_length: float = 0.0,
    ) -> None:
        """
        Define MPI point setting and dimensions along the x axis.

        Parameters
        ----------
        boxsize : float or list
            Size of the x-axis range.
        ngrid : int
            Number of divisions along this grid being defined.
        origin : float or list, optional
            Origin for x-axis range.
        buffer_length : float, optional
            Size of the buffer length.
        """
        self.boxsize = boxsize
        self.ngrid = ngrid
        self.origin = origin
        assert buffer_length >= 0.0, "buffer_length cannot be < 0."
        self.buffer_length = buffer_length


    def input(self, data: np.ndarray) -> None:
        """
        Input coordinates.

        Parameters
        ----------
        daata : array
            A 2 dimensional array where the first column is the x-axis coordinate.
        """
        self.data = data


    def limits4grid(self) -> None:
        """
        Determine the limits along x-axis for the x-axis division into a grid.
        """
        if np.isscalar(self.boxsize):
            xboxsize = self.boxsize
        else:
            xboxsize = self.boxsize[0]
        if np.isscalar(self.origin):
            xorigin = self.origin
        else:
            xorigin = self.origin[0]
        if np.isscalar(self.ngrid):
            xngrid = self.ngrid
        else:
            xngrid = self.ngrid[0]
        xedges, xgrid = shift.cart.mpi_grid1D(
            xboxsize, xngrid, self.MPI, origin=xorigin
        )
        self.limits = [xedges[0], xedges[-1]]
        self.ngrid_rank = len(xgrid)
        all_limits = self.MPI.gather(self.limits)
        if self.MPI.rank == 0:
            all_limits = np.vstack(all_limits)
        all_limits = self.MPI.bcast(all_limits)
        self.all_limits = all_limits


    def _checkifdist(self) -> bool:
        """
        Checks whether points are distributed according to slab ranges.

        Returns
        -------
        ifdist : bool
            Determines whether the data has ben distributed along the x-axis.
        """
        if self.data is None:
            check = False
        else:
            xmin, xmax = (
                np.min(self.data[:, 0]) - self.buffer_length,
                np.max(self.data[:, 0]) + self.buffer_length,
            )
            if self.limits[0] <= xmin and self.limits[1] >= xmax:
                check = True
            else:
                check = False
        return bool(self.MPI.comm.allreduce(check, op=self.MPI.mpi.LAND))
    

    def distribute(self, include_internalbuffer: bool = False) -> np.ndarray:
        """
        Distribute points along the x-axis using an MPI slab decomposition.

        Optionally include particles within the internal buffer regions.
        Periodic boundary conditions are not applied.

        Parameters
        ----------
        include_internalbuffer : bool, optional
            Whether to include internal buffer particles.

        Returns
        -------
        np.ndarray
            Distributed particle data.
        """
        if self._checkifdist():
            return self.data

        MPI = self.MPI

        # Determine the array dtype and shape across ranks.
        # Some ranks may initially have self.data = None.
        metadata = MPI.comm.allgather(
            None if self.data is None
            else (self.data.dtype, self.data.shape[1:])
        )

        metadata = [m for m in metadata if m is not None]

        if not metadata:
            return None

        dtype, trailing_shape = metadata[0]

        if any(m != metadata[0] for m in metadata):
            raise ValueError(
                "All ranks must have the same data dtype and trailing shape."
            )

        if self.data is None:
            data = np.empty((0, *trailing_shape), dtype=dtype)
        else:
            data = self.data

        # Identify particles destined for each rank.
        x = data[:, 0]

        sendcounts = np.zeros(MPI.size, dtype=np.int64)
        indices = []

        for i in range(MPI.size):

            xmin, xmax = self.all_limits[i]

            if include_internalbuffer:
                xmin -= self.buffer_length
                xmax += self.buffer_length
            else:
                if i == 0:
                    xmin -= self.buffer_length
                if i == MPI.size - 1:
                    xmax += self.buffer_length

            # The final slab includes its upper boundary.
            if i == MPI.size - 1:
                cond = (x >= xmin) & (x <= xmax)
            else:
                cond = (x >= xmin) & (x < xmax)

            idx = np.flatnonzero(cond)

            indices.append(idx)
            sendcounts[i] = idx.size

        # Group outgoing particles by destination rank.
        indices = np.concatenate(indices)
        senddata = np.ascontiguousarray(data[indices])

        # Exchange the number of particles.
        recvcounts = np.empty(MPI.size, dtype=np.int64)

        MPI.Alltoall(sendcounts, recvcounts)

        # Convert particle counts into array-element counts.
        row_size = int(np.prod(data.shape[1:], dtype=np.int64))

        sendcounts *= row_size
        recvcounts *= row_size

        # Compute displacements.
        senddispls = np.zeros(MPI.size, dtype=np.int64)
        recvdispls = np.zeros(MPI.size, dtype=np.int64)

        senddispls[1:] = np.cumsum(sendcounts[:-1])
        recvdispls[1:] = np.cumsum(recvcounts[:-1])

        # Allocate the receiving array.
        recvdata = np.empty(
            (int(recvcounts.sum() // row_size), *data.shape[1:]),
            dtype=dtype,
        )

        # Redistribute particles.
        MPI.Alltoallv(
            senddata.reshape(-1),
            sendcounts,
            senddispls,
            recvdata.reshape(-1),
            recvcounts,
            recvdispls,
        )

        return recvdata

    
    def distribute_grid2D(
        self, x2d: np.ndarray, y2d: np.ndarray, f2d: np.ndarray
    ) -> np.ndarray:
        """
        Distributes a 2D grid data set via slab decomposition.

        Parameters
        ----------
        x2d : ndarray
            X-coordinates for the 2D grid.
        y2d : ndarray
            Y-coordinates for the 2D grid.
        f2d : ndarray
            Field values for the 2D grid.

        Returns
        -------
        f : ndarray
            Distributed 2D field grid.
        """
        if x2d is not None:
            data = np.column_stack(
                [x2d.flatten(), y2d.flatten(), f2d.flatten()]
            )
        else:
            data = None

        self.input(data)
        data = self.distribute()

        if np.isscalar(self.boxsize):
            ybox = self.boxsize
        else:
            ybox = self.boxsize[1]

        nygrid = self.ngrid

        if np.isscalar(self.origin):
            yorigin = self.origin
        else:
            yorigin = self.origin[1]

        dx = (
            self.limits[1] - self.limits[0]
        ) / float(self.ngrid_rank)

        dy = ybox / float(nygrid)

        xpixs = np.floor(
            (data[:, 0] - self.limits[0]) / dx
        ).astype(np.int32)

        ypixs = np.floor(
            (data[:, 1] - yorigin) / dy
        ).astype(np.int32)

        f = np.zeros(
            (self.ngrid_rank, nygrid)
        )

        f[xpixs, ypixs] = data[:, 2]

        return f


    def distribute_grid3D(
        self,
        x3d: np.ndarray,
        y3d: np.ndarray,
        z3d: np.ndarray,
        f3d: np.ndarray,
    ) -> np.ndarray:
        """
        Distributes a 3D grid data set via slab decomposition.

        Parameters
        ----------
        x3d : ndarray
            X-coordinates for the 3D grid.
        y3d : ndarray
            Y-coordinates for the 3D grid.
        z3d : ndarray
            Z-coordinates for the 3D grid.
        f3d : ndarray
            Field values for the 3D grid.

        Returns
        -------
        f : ndarray
            Distributed 3D field grid.
        """
        if x3d is not None:
            data = np.column_stack(
                [
                    x3d.flatten(),
                    y3d.flatten(),
                    z3d.flatten(),
                    f3d.flatten(),
                ]
            )
        else:
            data = None

        self.input(data)
        data = self.distribute()

        if np.isscalar(self.boxsize):
            ybox = self.boxsize
            zbox = self.boxsize
        else:
            ybox = self.boxsize[1]
            zbox = self.boxsize[2]

        nygrid = self.ngrid
        nzgrid = self.ngrid

        if np.isscalar(self.origin):
            yorigin = self.origin
            zorigin = self.origin
        else:
            yorigin = self.origin[1]
            zorigin = self.origin[2]

        dx = (
            self.limits[1] - self.limits[0]
        ) / float(self.ngrid_rank)

        dy = ybox / float(nygrid)
        dz = zbox / float(nzgrid)

        xpixs = np.floor(
            (data[:, 0] - self.limits[0]) / dx
        ).astype(np.int32)

        ypixs = np.floor(
            (data[:, 1] - yorigin) / dy
        ).astype(np.int32)

        zpixs = np.floor(
            (data[:, 2] - zorigin) / dz
        ).astype(np.int32)

        f = np.zeros(
            (self.ngrid_rank, nygrid, nzgrid)
        )

        f[xpixs, ypixs, zpixs] = data[:, 3]

        return f

    
    def clean(self) -> None:
        """
        Reinitialises the class.
        """
        self.__init__(self.MPI)


def distribute_points_by_x(
    data: np.ndarray, boxsize: float, ngrid: int, origin: float, MPI: object
) -> np.ndarray:
    """
    Distributes points along the x-axis based grid-based slab decomposition.

    Parameters
    ----------
    data : array
        A 2 dimensional array where the first column is the x-axis coordinate.
    boxsize : float
        Box length.
    ngrid : int
        Number of divisions along this grid being defined.
    origin : float
        Origin for x-axis range.
    MPI : class object
        MPIutils MPI class object.

    Returns
    -------
    Distributed data.
    """
    Sorter = MPI_SortByX(MPI)
    Sorter.settings(boxsize, ngrid, origin)
    Sorter.input(data)
    Sorter.limits4grid()
    data = Sorter.distribute()
    return data
