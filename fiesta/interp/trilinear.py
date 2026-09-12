import numpy as np
from numba import njit

import shift

from .. import coords

from typing import Union, List


@njit
def trilinear_periodic(
    fgrid: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    xbox: float,
    ybox: float,
    zbox: float,
    originx: float = 0.0,
    originy: float = 0.0,
    originz: float = 0.0,
    dtype: np.dtype = np.float64
) -> np.ndarray: # pragma: no cover
    """
    Trilinear interpolation on a periodic 3D cell-centred grid.

    Parameters
    ----------
    fgrid : ndarray
        Three-dimensional array containing the field values on the grid.
    x, y, z : ndarray
        Coordinates of the points at which the field is interpolated.
    xbox, ybox, zbox : float
        Physical size of the grid domain along each axis.
    originx, originy, originz : float, optional
        Physical coordinates of the lower boundary of the grid domain.
    dtype : np.dtype, optional
        Data type of the returned interpolation array.

    Returns
    -------
    f : ndarray
        Interpolated field values at the coordinates ``(x, y, z)``.

    Notes
    -----
    The field is assumed to be defined at cell centres. Periodic boundary
    conditions are applied independently along all three axes.
    """

    ngridx, ngridy, ngridz = fgrid.shape
    npart = len(x)

    f = np.empty(npart, dtype=dtype)

    idx = ngridx / xbox
    idy = ngridy / ybox
    idz = ngridz / zbox

    for i in range(npart):

        gx = (x[i] - originx) * idx - 0.5
        gy = (y[i] - originy) * idy - 0.5
        gz = (z[i] - originz) * idz - 0.5

        ix1_raw = int(np.floor(gx))
        iy1_raw = int(np.floor(gy))
        iz1_raw = int(np.floor(gz))

        tx = gx - ix1_raw
        ty = gy - iy1_raw
        tz = gz - iz1_raw

        ix1 = ix1_raw % ngridx
        iy1 = iy1_raw % ngridy
        iz1 = iz1_raw % ngridz

        ix2 = (ix1 + 1) % ngridx
        iy2 = (iy1 + 1) % ngridy
        iz2 = (iz1 + 1) % ngridz

        f000 = fgrid[ix1, iy1, iz1]
        f100 = fgrid[ix2, iy1, iz1]
        f010 = fgrid[ix1, iy2, iz1]
        f110 = fgrid[ix2, iy2, iz1]

        f001 = fgrid[ix1, iy1, iz2]
        f101 = fgrid[ix2, iy1, iz2]
        f011 = fgrid[ix1, iy2, iz2]
        f111 = fgrid[ix2, iy2, iz2]

        f00 = (1.0 - tx)*f000 + tx*f100
        f10 = (1.0 - tx)*f010 + tx*f110

        f01 = (1.0 - tx)*f001 + tx*f101
        f11 = (1.0 - tx)*f011 + tx*f111

        f0 = (1.0 - ty)*f00 + ty*f10
        f1 = (1.0 - ty)*f01 + ty*f11

        f[i] = (1.0 - tz)*f0 + tz*f1

    return f


@njit
def trilinear_nonperiodic(
    fgrid: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    xbox: float,
    ybox: float,
    zbox: float,
    originx: float = 0.0,
    originy: float = 0.0,
    originz: float = 0.0,
    dtype: np.dtype = np.float64
) -> np.ndarray: # pragma: no cover
    """
    Trilinear interpolation on a non-periodic 3D cell-centred grid.

    Parameters
    ----------
    fgrid : ndarray
        Three-dimensional array containing the field values on the grid.
    x, y, z : ndarray
        Coordinates of the points at which the field is interpolated.
    xbox, ybox, zbox : float
        Physical size of the grid domain along each axis.
    originx, originy, originz : float, optional
        Physical coordinates of the lower boundary of the grid domain.
    dtype : np.dtype, optional
        Data type of the returned interpolation array.

    Returns
    -------
    f : ndarray
        Interpolated field values at the coordinates ``(x, y, z)``.

    Notes
    -----
    The field is assumed to be defined at cell centres. At non-periodic
    boundaries, interpolation is clamped to the nearest grid value whenever
    the interpolation stencil would extend outside the domain.
    """

    ngridx, ngridy, ngridz = fgrid.shape
    npart = len(x)

    f = np.empty(npart, dtype=dtype)

    idx = ngridx / xbox
    idy = ngridy / ybox
    idz = ngridz / zbox

    for i in range(npart):

        gx = (x[i] - originx) * idx - 0.5
        gy = (y[i] - originy) * idy - 0.5
        gz = (z[i] - originz) * idz - 0.5

        ix1_raw = int(np.floor(gx))
        iy1_raw = int(np.floor(gy))
        iz1_raw = int(np.floor(gz))

        tx = gx - ix1_raw
        ty = gy - iy1_raw
        tz = gz - iz1_raw

        # x-axis
        if ix1_raw < 0:
            ix1 = 0
            ix2 = 0
        elif ix1_raw >= ngridx - 1:
            ix1 = ngridx - 1
            ix2 = ngridx - 1
        else:
            ix1 = ix1_raw
            ix2 = ix1 + 1

        # y-axis
        if iy1_raw < 0:
            iy1 = 0
            iy2 = 0
        elif iy1_raw >= ngridy - 1:
            iy1 = ngridy - 1
            iy2 = ngridy - 1
        else:
            iy1 = iy1_raw
            iy2 = iy1 + 1

        # z-axis
        if iz1_raw < 0:
            iz1 = 0
            iz2 = 0
        elif iz1_raw >= ngridz - 1:
            iz1 = ngridz - 1
            iz2 = ngridz - 1
        else:
            iz1 = iz1_raw
            iz2 = iz1 + 1

        f000 = fgrid[ix1, iy1, iz1]
        f100 = fgrid[ix2, iy1, iz1]
        f010 = fgrid[ix1, iy2, iz1]
        f110 = fgrid[ix2, iy2, iz1]

        f001 = fgrid[ix1, iy1, iz2]
        f101 = fgrid[ix2, iy1, iz2]
        f011 = fgrid[ix1, iy2, iz2]
        f111 = fgrid[ix2, iy2, iz2]

        f00 = (1.0 - tx)*f000 + tx*f100
        f10 = (1.0 - tx)*f010 + tx*f110

        f01 = (1.0 - tx)*f001 + tx*f101
        f11 = (1.0 - tx)*f011 + tx*f111

        f0 = (1.0 - ty)*f00 + ty*f10
        f1 = (1.0 - ty)*f01 + ty*f11

        f[i] = (1.0 - tz)*f0 + tz*f1

    return f


@njit
def trilinear_axisperiodic(
    fgrid: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    xbox: float,
    ybox: float,
    zbox: float,
    perix: int,
    periy: int,
    periz: int,
    originx: float = 0.0,
    originy: float = 0.0,
    originz: float = 0.0,
    dtype: np.dtype = np.float64
) -> np.ndarray: # pragma: no cover
    """
    Trilinear interpolation on a 3D cell-centred grid with configurable
    periodicity along each axis.

    Parameters
    ----------
    fgrid : ndarray
        Three-dimensional array containing the field values on the grid.
    x, y, z : ndarray
        Coordinates of the points at which the field is interpolated.
    xbox, ybox, zbox : float
        Physical size of the grid domain along each axis.
    perix, periy, periz : int
        Periodicity flags for the x-, y-, and z-directions respectively.
        A value of 1 applies periodic boundary conditions, while 0 applies
        non-periodic boundary conditions.
    originx, originy, originz : float, optional
        Physical coordinates of the lower boundary of the grid domain.
    dtype : np.dtype, optional
        Data type of the returned interpolation array.

    Returns
    -------
    f : ndarray
        Interpolated field values at the coordinates ``(x, y, z)``.

    Notes
    -----
    Periodic axes wrap across the corresponding domain boundary, while
    non-periodic axes are clamped to the nearest grid value whenever the
    interpolation stencil would extend outside the domain.
    """

    ngridx, ngridy, ngridz = fgrid.shape
    npart = len(x)

    f = np.empty(npart, dtype=dtype)

    idx = ngridx / xbox
    idy = ngridy / ybox
    idz = ngridz / zbox

    for i in range(npart):

        gx = (x[i] - originx) * idx - 0.5
        gy = (y[i] - originy) * idy - 0.5
        gz = (z[i] - originz) * idz - 0.5

        ix1_raw = int(np.floor(gx))
        iy1_raw = int(np.floor(gy))
        iz1_raw = int(np.floor(gz))

        tx = gx - ix1_raw
        ty = gy - iy1_raw
        tz = gz - iz1_raw

        # x-axis
        if perix == 1:
            ix1 = ix1_raw % ngridx
            ix2 = (ix1 + 1) % ngridx
        else:
            if ix1_raw < 0:
                ix1 = 0
                ix2 = 0
            elif ix1_raw >= ngridx - 1:
                ix1 = ngridx - 1
                ix2 = ngridx - 1
            else:
                ix1 = ix1_raw
                ix2 = ix1 + 1

        # y-axis
        if periy == 1:
            iy1 = iy1_raw % ngridy
            iy2 = (iy1 + 1) % ngridy
        else:
            if iy1_raw < 0:
                iy1 = 0
                iy2 = 0
            elif iy1_raw >= ngridy - 1:
                iy1 = ngridy - 1
                iy2 = ngridy - 1
            else:
                iy1 = iy1_raw
                iy2 = iy1 + 1

        # z-axis
        if periz == 1:
            iz1 = iz1_raw % ngridz
            iz2 = (iz1 + 1) % ngridz
        else:
            if iz1_raw < 0:
                iz1 = 0
                iz2 = 0
            elif iz1_raw >= ngridz - 1:
                iz1 = ngridz - 1
                iz2 = ngridz - 1
            else:
                iz1 = iz1_raw
                iz2 = iz1 + 1

        f000 = fgrid[ix1, iy1, iz1]
        f100 = fgrid[ix2, iy1, iz1]
        f010 = fgrid[ix1, iy2, iz1]
        f110 = fgrid[ix2, iy2, iz1]

        f001 = fgrid[ix1, iy1, iz2]
        f101 = fgrid[ix2, iy1, iz2]
        f011 = fgrid[ix1, iy2, iz2]
        f111 = fgrid[ix2, iy2, iz2]

        f00 = (1.0 - tx)*f000 + tx*f100
        f10 = (1.0 - tx)*f010 + tx*f110

        f01 = (1.0 - tx)*f001 + tx*f101
        f11 = (1.0 - tx)*f011 + tx*f111

        f0 = (1.0 - ty)*f00 + ty*f10
        f1 = (1.0 - ty)*f01 + ty*f11

        f[i] = (1.0 - tz)*f0 + tz*f1

    return f


def trilinear(
    fgrid: np.ndarray,
    boxsize: Union[float, List[float]],
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    origin: Union[float, List[float]] = 0.0,
    fill_value: float = np.nan,
    periodic: bool = True,
    dtype: np.dtype = np.float64
) -> np.ndarray:
    """Trilinear interpolation from a 3D grid defined in box of [0., boxsize].

    Parameter
    ---------
    fgrid : array
        Field values on a 3D grid.
    boxsize : float or list
        Box size in one or all axes.
    x : array
        x coordinate values.
    y : array
        y coordinate values.
    z : array
        z coordinate values.
    origin : float or list, optional
        Origin for the axes.
    fill_value : float, optional
        Fill outside boundary values.
    periodic : bool, optional
        Determines whether to interpolate on a periodic grid.
    dtype : np.dtype, optional
        Data type of the returned interpolation array.

    Returns
    -------
    f : array
        Field interpolation values.
    """
    if np.isscalar(boxsize):
        xbox = boxsize
        ybox = boxsize
        zbox = boxsize
    else:
        xbox, ybox, zbox = boxsize[0], boxsize[1], boxsize[2]
    if np.isscalar(origin):
        originx = origin
        originy = origin
        originz = origin
    else:
        originx = origin[0]
        originy = origin[1]
        originz = origin[2]
    # check if particles are inside the box
    inside = (
        (x >= originx)
        & (x < originx + xbox)
        & (y >= originy)
        & (y < originy + ybox)
        & (z >= originz)
        & (z < originz + zbox)
    )
    if np.all(inside):
        # All particles are within the boundaries so no boundary management is necessary.
        if np.isscalar(periodic):
            if periodic == True:
                f = trilinear_periodic(fgrid, x, y, z, xbox, ybox, zbox, originx, originy, originz, dtype=dtype)
            else:
                f = trilinear_nonperiodic(fgrid, x, y, z, xbox, ybox, zbox, originx, originy, originz, dtype=dtype)
        else:
            if periodic[0] is True:
                perix = 1
            else:
                perix = 0
            if periodic[1] is True:
                periy = 1
            else:
                periy = 0
            if periodic[2] is True:
                periz = 1
            else:
                periz = 0
            f = trilinear_axisperiodic(fgrid, x, y, z, xbox, ybox, zbox, perix, periy, periz, originx, originy, originz, dtype=dtype)
    else:
        f = np.full(len(x), fill_value, dtype=dtype)
        if np.isscalar(periodic):
            if periodic == True:
                f[inside] = trilinear_periodic(fgrid, x[inside], y[inside], z[inside], xbox, ybox, zbox, originx, originy, originz, dtype=dtype)
            else:
                f[inside] = trilinear_nonperiodic(fgrid, x[inside], y[inside], z[inside], xbox, ybox, zbox, originx, originy, originz, dtype=dtype)
        else:
            if periodic[0] is True:
                perix = 1
            else:
                perix = 0
            if periodic[1] is True:
                periy = 1
            else:
                periy = 0
            if periodic[2] is True:
                periz = 1
            else:
                periz = 0
            f[inside] = trilinear_axisperiodic(fgrid, x[inside], y[inside], z[inside], xbox, ybox, zbox, perix, periy, periz, originx, originy, originz, dtype=dtype)
    return f


def mpi_trilinear(
    fgrid: np.ndarray,
    ngrid: Union[int, List[int]],
    boxsize: Union[float, List[float]],
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    MPI: object,
    origin: Union[float, List[float]] = 0.0,
    fill_value: float = np.nan,
    periodic: bool = True,
    check_distributed: bool = True,
    dtype: np.dtype = np.float64
) -> np.ndarray:
    """
    Trilinear interpolation from a 3D grid defined in box of [0., boxsize].

    Parameter
    ---------
    fgrid : array
        Field values on a 3D grid
    ngrid : int or int list
        Grid dimensions.
    boxsize : float or list
        Box size.
    x : array
        x coordinate values.
    y : array
        y coordinate values.
    z : array
        z coordinate values.
    MPI : class object
        shift.mpiutils MPI object.
    origin : float or list, optional
        Origin for x and y coordinates.
    fill_value : float, optional
        Fill outside boundary values.
    periodic : bool, optional
        Determines whether to interpolate on a periodic grid.
    check_distributed : bool, optional
        If True, checks if the coordinates are distributed correctly across MPI ranks.
    dtype : np.dtype
        Data type for the output array.

    Returns
    -------
    x : array
        x coordinate values, redistributed for specific slab.
    y : array
        y coordinate values, redistributed for specific slab.
    z : array
        z coordinate values, redistributed for specific slab.
    f : array
        Field interpolation values.
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

    # define grid on each axis
    if np.isscalar(periodic):
        xperiodic, yperiodic, zperiodic = periodic, periodic, periodic
    else:
        xperiodic, yperiodic, zperiodic = periodic[0], periodic[1], periodic[2]
    
    # define pixel length across each axis
    dx = xboxsize / nxgrid

    xedges, xgrid = shift.cart.mpi_grid1D(xboxsize, nxgrid, MPI, origin=xorigin)

    _xorigin = xedges[0]
    _xboxsize = xedges[-1] - xedges[0]

    if np.shape(fgrid)[0] == len(xgrid):
        correct_shape = True
    else:
        correct_shape = False

    correct_shapes = MPI.collect([correct_shape])

    if MPI.rank == 0:
        if all(correct_shapes):
            correct_shape = True
        else:
            correct_shape = False
    
    correct_shape = MPI.broadcast(correct_shape)
    
    if correct_shape:

        fgrid_sendup = np.array([MPI.send_up(fgrid[-1])])
        fgrid_senddown = np.array([MPI.send_down(fgrid[0])])

        if MPI.rank == 0:
            if xperiodic:
                fgrid = np.concatenate([fgrid_sendup, fgrid, fgrid_senddown])
                _xorigin -= dx
                _xboxsize += 2*dx
            else:
                fgrid = np.concatenate([fgrid, fgrid_senddown])
                _xboxsize += dx
        elif MPI.rank == MPI.size - 1:
            if xperiodic:
                fgrid = np.concatenate([fgrid_sendup, fgrid, fgrid_senddown])
                _xorigin -= dx
                _xboxsize += 2*dx
            else:
                fgrid = np.concatenate([fgrid_sendup, fgrid])
                _xorigin -= dx
                _xboxsize += dx
        else:
            fgrid = np.concatenate([fgrid_sendup, fgrid, fgrid_senddown])
            _xorigin -= dx
            _xboxsize += 2*dx
        
        if check_distributed:
            if x is not None:
                data = coords.xyz2points(x, y, z)
            else: 
                data = None
            data = coords.distribute_points_by_x(data, boxsize, ngrid, origin, MPI)
            x, y, z = coords.points2xyz(data)
        f = trilinear(
            fgrid, [_xboxsize, yboxsize, zboxsize], x, y, z, [_xorigin, yorigin, zorigin], 
            fill_value=fill_value, periodic=[False, yperiodic, zperiodic], dtype=dtype
        )
        return x, y, z, f
    else:
        MPI.mpi_print_zero("ERROR: Shape of fgrid does not match expectation for distributed array")
        return np.nan, np.nan, np.nan, np.nan
        

        


