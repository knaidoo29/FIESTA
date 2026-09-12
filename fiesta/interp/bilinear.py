import numpy as np
from numba import njit

import shift

from .. import coords

from typing import Union, List


@njit
def bilinear_periodic(
    fgrid: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    xbox: float,
    ybox: float,
    originx: float = 0.0,
    originy: float = 0.0,
    dtype: np.dtype = np.float64
) -> np.ndarray: # pragma: no cover
    """
    Bilinear interpolation on a periodic 2D cell-centred grid.

    Parameters
    ----------
    fgrid : ndarray
        Two-dimensional array containing the field values on the grid.
    x, y : ndarray
        Coordinates of the points at which the field is interpolated.
    xbox, ybox : float
        Physical size of the grid domain along the x- and y-directions.
    originx, originy : float, optional
        Physical coordinates of the lower boundary of the grid domain.
    dtype : np.dtype, optional
        Data type of the returned interpolation array.

    Returns
    -------
    f : ndarray
        Interpolated field values at the coordinates ``(x, y)``.

    Notes
    -----
    The field is assumed to be defined at cell centres. Periodic boundary
    conditions are applied independently along both axes. Coordinates are
    converted to dimensionless grid coordinates before identifying the
    neighbouring cells and interpolation weights.
    """

    ngridx, ngridy = fgrid.shape
    npart = len(x)

    f = np.empty(npart, dtype=dtype)

    idx = ngridx / xbox
    idy = ngridy / ybox

    for i in range(npart):

        gx = (x[i] - originx) * idx - 0.5
        gy = (y[i] - originy) * idy - 0.5

        ix1_raw = int(np.floor(gx))
        iy1_raw = int(np.floor(gy))

        tx = gx - ix1_raw
        ty = gy - iy1_raw

        ix1 = ix1_raw % ngridx
        iy1 = iy1_raw % ngridy

        ix2 = (ix1 + 1) % ngridx
        iy2 = (iy1 + 1) % ngridy

        f11 = fgrid[ix1, iy1]
        f12 = fgrid[ix2, iy1]
        f21 = fgrid[ix1, iy2]
        f22 = fgrid[ix2, iy2]

        f1 = (1.0 - tx)*f11 + tx*f12
        f2 = (1.0 - tx)*f21 + tx*f22

        f[i] = (1.0 - ty)*f1 + ty*f2

    return f


@njit
def bilinear_nonperiodic(
    fgrid: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    xbox: float,
    ybox: float,
    originx: float = 0.0,
    originy: float = 0.0,
    dtype: np.dtype = np.float64
) -> np.ndarray: # pragma: no cover
    """
    Bilinear interpolation on a non-periodic 2D cell-centred grid.

    Parameters
    ----------
    fgrid : ndarray
        Two-dimensional array containing the field values on the grid.
    x, y : ndarray
        Coordinates of the points at which the field is interpolated.
    xbox, ybox : float
        Physical size of the grid domain along the x- and y-directions.
    originx, originy : float, optional
        Physical coordinates of the lower boundary of the grid domain.
    dtype : np.dtype, optional
        Data type of the returned interpolation array.

    Returns
    -------
    f : ndarray
        Interpolated field values at the coordinates ``(x, y)``.

    Notes
    -----
    The field is assumed to be defined at cell centres. At non-periodic
    boundaries, interpolation is clamped to the nearest grid value whenever
    a neighbouring interpolation cell would lie outside the domain.
    Coordinates are converted to dimensionless grid coordinates before
    identifying the neighbouring cells and interpolation weights.
    """
    ngridx, ngridy = fgrid.shape
    npart = len(x)

    f = np.empty(npart, dtype=dtype)

    idx = ngridx / xbox
    idy = ngridy / ybox

    for i in range(npart):

        gx = (x[i] - originx) * idx - 0.5
        gy = (y[i] - originy) * idy - 0.5

        # x direction
        ix1_raw = int(np.floor(gx))
        tx = gx - ix1_raw

        if ix1_raw < 0:
            ix1 = 0
            ix2 = 0

        elif ix1_raw >= ngridx - 1:
            ix1 = ngridx - 1
            ix2 = ngridx - 1

        else:
            ix1 = ix1_raw
            ix2 = ix1 + 1

        # y direction
        iy1_raw = int(np.floor(gy))
        ty = gy - iy1_raw

        if iy1_raw < 0:
            iy1 = 0
            iy2 = 0

        elif iy1_raw >= ngridy - 1:
            iy1 = ngridy - 1
            iy2 = ngridy - 1

        else:
            iy1 = iy1_raw
            iy2 = iy1 + 1

        f11 = fgrid[ix1, iy1]
        f12 = fgrid[ix2, iy1]
        f21 = fgrid[ix1, iy2]
        f22 = fgrid[ix2, iy2]

        f1 = (1.0 - tx)*f11 + tx*f12
        f2 = (1.0 - tx)*f21 + tx*f22

        f[i] = (1.0 - ty)*f1 + ty*f2

    return f


@njit
def bilinear_axisperiodic(
    fgrid: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    xbox: float,
    ybox: float,
    perix: int,
    periy: int,
    originx: float = 0.0,
    originy: float = 0.0,
    dtype: np.dtype = np.float64
) -> np.ndarray: # pragma: no cover
    """
    Bilinear interpolation on a 2D cell-centred grid with configurable
    periodicity along each axis.

    Parameters
    ----------
    fgrid : ndarray
        Two-dimensional array containing the field values on the grid.
    x, y : ndarray
        Coordinates of the points at which the field is interpolated.
    xbox, ybox : float
        Physical size of the grid domain along the x- and y-directions.
    perix, periy : int
        Periodicity flags for the x- and y-directions respectively.
        A value of 1 applies periodic boundary conditions, while 0 applies
        non-periodic boundary conditions.
    originx, originy : float, optional
        Physical coordinates of the lower boundary of the grid domain.
    dtype : np.dtype, optional
        Data type of the returned interpolation array.

    Returns
    -------
    f : ndarray
        Interpolated field values at the coordinates ``(x, y)``.

    Notes
    -----
    The field is assumed to be defined at cell centres. Periodic axes wrap
    across the corresponding domain boundary, while non-periodic axes are
    clamped to the nearest grid value whenever the interpolation stencil
    would extend outside the domain. Coordinates are converted to
    dimensionless grid coordinates before identifying neighbouring cells
    and interpolation weights.
    """

    ngridx, ngridy = fgrid.shape
    npart = len(x)

    f = np.empty(npart, dtype=dtype)

    idx = ngridx / xbox
    idy = ngridy / ybox

    for i in range(npart):

        gx = (x[i] - originx) * idx - 0.5
        gy = (y[i] - originy) * idy - 0.5

        # x
        ix1_raw = int(np.floor(gx))
        tx = gx - ix1_raw

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

        # y
        iy1_raw = int(np.floor(gy))
        ty = gy - iy1_raw

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

        f11 = fgrid[ix1, iy1]
        f12 = fgrid[ix2, iy1]
        f21 = fgrid[ix1, iy2]
        f22 = fgrid[ix2, iy2]

        f1 = (1.0 - tx)*f11 + tx*f12
        f2 = (1.0 - tx)*f21 + tx*f22

        f[i] = (1.0 - ty)*f1 + ty*f2

    return f


def bilinear(
    fgrid: np.ndarray,
    boxsize: Union[float, List[float]],
    x: np.ndarray,
    y: np.ndarray,
    origin: Union[float, List[float]] = 0.0,
    fill_value: float = np.nan,
    periodic: bool = True,
    dtype: np.dtype = np.float64
) -> np.ndarray:
    """
    Bilinear interpolation from a 2D grid defined in box of [0., boxsize].

    Parameter
    ---------
    fgrid : array
        Field values on a 2D grid.
    boxsize : float or list
        Box size.
    x : array
        x coordinate values.
    y : array
        y coordinate values.
    origin : float or list, optional
        Origin for x and y coordinates.
    fill_value : float, optional
        Fill outside boundary values.
    periodic : bool, optional
        Determines whether to interpolate on a periodic grid.
    dtype : np.dtype
        Data type for the output array.
    
    Returns
    -------
    f : array
        Field interpolation values.
    """
    if np.isscalar(boxsize):
        xbox = boxsize
        ybox = boxsize
    else:
        xbox, ybox = boxsize[0], boxsize[1]
    if np.isscalar(origin):
        originx = origin
        originy = origin
    else:
        originx = origin[0]
        originy = origin[1]
    # check if particles are inside the box
    inside = (
        (x >= originx)
        & (x < originx + xbox)
        & (y >= originy)
        & (y < originy + ybox)
    )
    if np.all(inside):
        # All particles are within the boundaries so no boundary management is necessary.
        if np.isscalar(periodic):
            if periodic == True:
                f = bilinear_periodic(fgrid, x, y, xbox, ybox, originx, originy, dtype=dtype)
            else:
                f = bilinear_nonperiodic(fgrid, x, y, xbox, ybox, originx, originy, dtype=dtype)
        else:
            if periodic[0] is True:
                perix = 1
            else:
                perix = 0
            if periodic[1] is True:
                periy = 1
            else:
                periy = 0
            f = bilinear_axisperiodic(fgrid, x, y, xbox, ybox, perix, periy, originx, originy, dtype=dtype)
    else:
        # Some particles are outside the boundary.
        f = np.full(len(x), fill_value, dtype=dtype)
        if np.isscalar(periodic):
            if periodic == True:
                f[inside] = bilinear_periodic(fgrid, x[inside], y[inside], xbox, ybox, originx, originy, dtype=dtype)
            else:
                f[inside] = bilinear_nonperiodic(fgrid, x[inside], y[inside], xbox, ybox, originx, originy, dtype=dtype)
        else:
            if periodic[0] is True:
                perix = 1
            else:
                perix = 0
            if periodic[1] is True:
                periy = 1
            else:
                periy = 0
            f[inside] = bilinear_axisperiodic(fgrid, x[inside], y[inside], xbox, ybox, perix, periy, originx, originy, dtype=dtype)
    return f


def mpi_bilinear(
    fgrid: np.ndarray,
    ngrid: Union[int, List[int]],
    boxsize: Union[float, List[float]],
    x: np.ndarray,
    y: np.ndarray,
    MPI: object,
    origin: Union[float, List[float]] = 0.0,
    fill_value: float = np.nan,
    periodic: bool = True,
    check_distributed: bool = True,
    dtype: np.dtype = np.float64
) -> np.ndarray:
    """
    Bilinear interpolation from a 2D grid defined in box of [0., boxsize].

    Parameter
    ---------
    fgrid : array
        Field values on a 2D grid
    ngrid : int or int list
        Grid dimensions.
    boxsize : float or list
        Box size.
    x : array
        x coordinate values.
    y : array
        y coordinate values.
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
    f : array
        Field interpolation values.
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

    # define grid on each axis
    if np.isscalar(periodic):
        xperiodic, yperiodic = periodic, periodic
    else:
        xperiodic, yperiodic = periodic[0], periodic[1]
    
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
                data = coords.xy2points(x, y)
            else:
                data = None
            data = coords.distribute_points_by_x(data, boxsize, ngrid, origin, MPI)
            x, y = coords.points2xy(data)
        f = bilinear(
            fgrid, [_xboxsize, yboxsize], x, y, [_xorigin, yorigin], fill_value=fill_value, 
            periodic=[False, yperiodic], dtype=dtype
        )
        return x, y, f
    else:
        MPI.mpi_print_zero("ERROR: Shape of fgrid does not match expectation for distributed array")
        return np.nan, np.nan, np.nan
        

        


