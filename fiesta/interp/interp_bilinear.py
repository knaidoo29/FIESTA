import numpy as np
from numba import njit

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
