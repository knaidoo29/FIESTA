import numpy as np

import shift

from .. import coords, boundary
from . import sph4grid

from typing import Union, Tuple, List, Optional


def mpi_sph4grid2D(
    x: np.ndarray,
    y: np.ndarray,
    boxsize: Union[float, List[float]],
    ngrid: Union[int, List[int]],
    MPI: object,
    f: Optional[np.ndarray] = None,
    origin: Union[float, List[float]] = 0.0,
    mass: Optional[np.ndarray] = None,
    k: int = 20,
    periodic: bool = True,
    subsampling: int = 1,
    outputgrid: bool = False,
    buffer_length: float = 0.0,
) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """
    Estimates a field on a regular grid based on SPH k neighbours.

    Parameters
    ----------
    x, y : array
        Cartesian coordinates of the input points.
    boxsize : float or list of float
        Size of the box in each dimension.
    ngrid : int or list of int
        Number of grid points in each dimension.
    MPI : object
        shift.mpiutils MPI object.
    f : array, optional
        Field values at the input points.
    origin : float or list of float, optional
        Origin of the grid.
    mass: array, optional
        Mass of the input points. If not provided, it is assumed to be 1 for all points.
    k : int, optional
        Number of nearest neighbors to consider.
    periodic : bool or list of bool, optional
        Whether to use periodic boundary conditions.
    subsampling : int, optional
        Subsampling factor for the output grid.
    outputgrid : bool, optional
        Whether to output the grid coordinates.
    buffer_length : float, optional
        Buffer length, must be smaller than the slab x-axis length.

    Returns
    -------
    dgrid : array
        Density field on the regular grid, if f is None.
    fgrid : array
        Estimated field on the regular grid, if f is provided.
    x2d, y2d : array
        Grid coordinates, if outputgrid is True.
    """

    # define boxsize on each axis
    if np.isscalar(boxsize):
        xboxsize, yboxsize = boxsize, boxsize
    else:
        xboxsize, yboxsize = boxsize[0], boxsize[1]

    # define boxsize on each axis
    if np.isscalar(origin):
        xorigin, yorigin = origin, origin
    else:
        xorigin, yorigin = origin[0], origin[1]

    # define grid on each axis
    if np.isscalar(ngrid):
        nxgrid, nygrid = ngrid, ngrid
    else:
        nxgrid, nygrid = ngrid[0], ngrid[1]

    # collapse data
    if x is not None:
        if f is None:
            if mass is None:
                data = coords.coord2points([x, y, np.ones(len(x))])
            else:
                data = coords.coord2points([x, y, mass])
        else:
            if mass is None:
                data = coords.coord2points([x, y, f, np.ones(len(x))])
            else:
                data = coords.coord2points([x, y, f, mass])
    else:  # pragma: no cover
        data = None

    if f is not None:
        fieldexist = True
    else:
        fieldexist = False

    fieldexists = MPI.collect(fieldexist)

    if MPI.rank == 0:
        fieldexist = any(fieldexists)
    else:
        fieldexist = None

    fieldexist = MPI.broadcast(fieldexist)

    # sort coordinates and distribute by coordinate system
    MPI_SBX = coords.MPI_SortByX(MPI)
    MPI_SBX.settings(xboxsize, nxgrid, origin=xorigin)
    MPI_SBX.input(data)
    MPI_SBX.limits4grid()
    data = MPI_SBX.distribute()

    assert buffer_length <= xboxsize, "buffer_length must be smaller than slab length."

    if periodic:
        data = boundary.mpi_buffer_periodic_2D(
            data, boxsize, buffer_length, MPI, origin, include_internal=True
        )
    else:
        data = boundary.mpi_buffer_internal_2D(
            data, boxsize, buffer_length, MPI, origin
        )

    # coordinates only inside slab
    x, y = data[:, 0], data[:, 1]
    if fieldexist:
        f = data[:, 2]
        mass = data[:, 3]
    else:
        mass = data[:, 2]

    xedges, xgrid = shift.cart.mpi_grid1D(xboxsize, nxgrid, MPI, origin=xorigin)

    # local grid, we will use an underscore to differentiate from global
    _xorigin = xedges[0]
    _xboxsize = xedges[-1] - xedges[0]

    _nxgrid = len(xgrid)

    return sph4grid.sph4grid2D(
        x,
        y,
        [_xboxsize, yboxsize],
        [_nxgrid, nygrid],
        f=f,
        origin=[_xorigin, yorigin],
        mass=mass,
        k=k,
        periodic=False,
        subsampling=subsampling,
        outputgrid=outputgrid,
    )


def mpi_sph4grid3D(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    boxsize: Union[float, List[float]],
    ngrid: Union[int, List[int]],
    MPI: object,
    f: Optional[np.ndarray] = None,
    origin: Union[float, List[float]] = 0.0,
    mass: Optional[np.ndarray] = None,
    k: int = 20,
    periodic: bool = True,
    subsampling: int = 1,
    outputgrid: bool = False,
    buffer_length: float = 0.0,
) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """
    Estimates a field on a regular grid based on SPH k neighbours.

    Parameters
    ----------
    x, y, z : array
        Cartesian coordinates of the input points.
    boxsize : float or list of float
        Size of the box in each dimension.
    ngrid : int or list of int
        Number of grid points in each dimension.
    MPI : object
        shift.mpiutils MPI object.
    f : array, optional
        Field values at the input points.
    origin : float or list of float, optional
        Origin of the grid.
    mass: array, optional
        Mass of the input points. If not provided, it is assumed to be 1 for all points.
    k : int, optional
        Number of nearest neighbors to consider.
    periodic : bool or list of bool, optional
        Whether to use periodic boundary conditions.
    subsampling : int, optional
        Subsampling factor for the output grid.
    outputgrid : bool, optional
        Whether to output the grid coordinates.
    buffer_length : float, optional
        Buffer length, must be smaller than the slab x-axis length.

    Returns
    -------
    dgrid : array
        Density field on the regular grid, if f is None.
    fgrid : array
        Estimated field on the regular grid, if f is provided.
    x2d, y2d : array
        Grid coordinates, if outputgrid is True.
    """

    # define boxsize on each axis
    if np.isscalar(boxsize):
        xboxsize, yboxsize, zboxsize = boxsize, boxsize, boxsize
    else:
        xboxsize, yboxsize, zboxsize = boxsize[0], boxsize[1], boxsize[2]

    # define boxsize on each axis
    if np.isscalar(origin):
        xorigin, yorigin, zorigin = origin, origin, origin
    else:
        xorigin, yorigin, zorigin = origin[0], origin[1], origin[2]

    # define grid on each axis
    if np.isscalar(ngrid):
        nxgrid, nygrid, nzgrid = ngrid, ngrid, ngrid
    else:
        nxgrid, nygrid, nzgrid = ngrid[0], ngrid[1], ngrid[2]

    # collapse data
    if x is not None:
        if f is None:
            if mass is None:
                data = coords.coord2points([x, y, z, np.ones(len(x))])
            else:
                data = coords.coord2points([x, y, z, mass])
        else:
            if mass is None:
                data = coords.coord2points([x, y, z, f, np.ones(len(x))])
            else:
                data = coords.coord2points([x, y, z, f, mass])
    else:  # pragma: no cover
        data = None

    if f is not None:
        fieldexist = True
    else:
        fieldexist = False

    fieldexists = MPI.collect(fieldexist)

    if MPI.rank == 0:
        fieldexist = any(fieldexists)
    else:
        fieldexist = None

    fieldexist = MPI.broadcast(fieldexist)

    # sort coordinates and distribute by coordinate system
    MPI_SBX = coords.MPI_SortByX(MPI)
    MPI_SBX.settings(xboxsize, nxgrid, origin=xorigin)
    MPI_SBX.input(data)
    MPI_SBX.limits4grid()
    data = MPI_SBX.distribute()

    assert buffer_length <= xboxsize, "buffer_length must be smaller than slab length."

    if periodic:
        data = boundary.mpi_buffer_periodic_3D(
            data, boxsize, buffer_length, MPI, origin, include_internal=True
        )
    else:
        data = boundary.mpi_buffer_internal_3D(
            data, boxsize, buffer_length, MPI, origin
        )

    # coordinates only inside slab
    x, y, z = data[:, 0], data[:, 1], data[:, 2]
    if fieldexist:
        f = data[:, 3]
        mass = data[:, 4]
    else:
        mass = data[:, 3]

    xedges, xgrid = shift.cart.mpi_grid1D(xboxsize, nxgrid, MPI, origin=xorigin)

    # local grid, we will use an underscore to differentiate from global
    _xorigin = xedges[0]
    _xboxsize = xedges[-1] - xedges[0]

    _nxgrid = len(xgrid)

    return sph4grid.sph4grid3D(
        x,
        y,
        z,
        [_xboxsize, yboxsize, zboxsize],
        [_nxgrid, nygrid, nzgrid],
        f=f,
        origin=[_xorigin, yorigin, zorigin],
        mass=mass,
        k=k,
        periodic=False,
        subsampling=subsampling,
        outputgrid=outputgrid,
    )
