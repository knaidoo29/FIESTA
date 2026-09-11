import numpy as np

from numba import njit
# from .. import src

from typing import Union, List

@njit
def part2grid_ngp_2d(
    x: np.ndarray,
    y: np.ndarray,
    f: np.ndarray,
    xlength: float,
    ylength: float,
    xmin: float,
    ymin: float,
    nxgrid: int,
    nygrid: int,
    dtype: np.dtype = np.float64
) -> np.ndarray:
    """
    Nearest-grid-point assignment in 2D.

    Parameters
    ----------
    x, y : array
        Cartesian coordinates.
    f : array
        Field values at the particle coordinates.
    xlength, ylength : float
        Length of the box along the x and y coordinates.
    xmin, ymin : float
        Minimum values along the x and y axes.
    nxgrid, nygrid : int
        Number of grid cells along x and y.
    dtype : np.dtype, optional
        Data type of the returned interpolation array.

    Returns
    -------
    fgrid : ndarray
        NGP field assignment.
    """
    npart = len(x)

    idx = nxgrid / xlength
    idy = nygrid / ylength

    wngp = idx * idy

    fgrid = np.zeros((nxgrid, nygrid), dtype=dtype)

    for i in range(npart):

        ix = int(np.floor((x[i] - xmin) * idx))
        iy = int(np.floor((y[i] - ymin) * idy))

        if (
            (ix >= 0)
            and (ix < nxgrid)
            and (iy >= 0)
            and (iy < nygrid)
        ):
            fgrid[ix, iy] += f[i] * wngp

    return fgrid



@njit
def part2grid_ngp_2d_unit(
    x: np.ndarray,
    y: np.ndarray,
    xlength: float,
    ylength: float,
    xmin: float,
    ymin: float,
    nxgrid: int,
    nygrid: int,
    dtype: np.dtype = np.float64
) -> np.ndarray:
    """
    Nearest-grid-point assignment in 2D.

    Parameters
    ----------
    x, y : array
        Cartesian coordinates.
    xlength, ylength : float
        Length of the box along the x and y coordinates.
    xmin, ymin : float
        Minimum values along the x and y axes.
    nxgrid, nygrid : int
        Number of grid cells along x and y.
    dtype : np.dtype, optional
        Data type of the returned interpolation array.
        
    Returns
    -------
    fgrid : ndarray
        NGP field assignment.
    """
    npart = len(x)

    idx = nxgrid / xlength
    idy = nygrid / ylength

    wngp = idx * idy

    fgrid = np.zeros((nxgrid, nygrid), dtype=dtype)

    for i in range(npart):

        ix = int(np.floor((x[i] - xmin) * idx))
        iy = int(np.floor((y[i] - ymin) * idy))

        if (
            (ix >= 0)
            and (ix < nxgrid)
            and (iy >= 0)
            and (iy < nygrid)
        ):
            fgrid[ix, iy] += wngp

    return fgrid


@njit
def part2grid_cic_2d(
    x: np.ndarray,
    y: np.ndarray,
    f: np.ndarray,
    xlength: float,
    ylength: float,
    xmin: float,
    ymin: float,
    nxgrid: int,
    nygrid: int,
    periodx: bool,
    periody: bool,
    dtype: np.dtype = np.float64
) -> np.ndarray:
    """
    Assign particle values to a two-dimensional grid using the
    cloud-in-cell (CIC) assignment scheme.

    Each particle contributes to the two nearest cell centres along each
    coordinate axis, giving a total stencil of four neighbouring grid
    cells. The contribution to each cell is weighted linearly according
    to the particle position relative to the neighbouring cell centres.

    The deposited field is normalized by the grid-cell area. When the
    complete stencil lies within the domain, or periodic boundaries are
    applied, the weights associated with each particle sum to

        1 / (dx * dy).

    Parameters
    ----------
    x, y : ndarray
        One-dimensional arrays containing the Cartesian coordinates of
        the particles. Both arrays must have the same length.
    f : ndarray
        One-dimensional array containing the value associated with each
        particle. Must have the same length as ``x`` and ``y``.
    xlength, ylength : float
        Physical lengths of the grid domain along the x and y axes.
    xmin, ymin : float
        Minimum coordinates of the grid domain along the x and y axes.
    nxgrid, nygrid : int
        Number of grid cells along the x and y axes.
    periodx, periody : bool
        Whether periodic boundary conditions are applied along the x
        and y axes, respectively.
    dtype : numpy dtype, optional
        Data type used for the output grid. Default is ``np.float64``.

    Returns
    -------
    fgrid : ndarray
        Two-dimensional array of shape ``(nxgrid, nygrid)`` containing
        the CIC-assigned field, with data type given by ``dtype``.

    Notes
    -----
    The grid is treated as cell-centred. The dimensionless particle
    coordinate relative to the grid-cell centres is

        gx = (x - xmin) / dx - 0.5
        gy = (y - ymin) / dy - 0.5,

    where

        dx = xlength / nxgrid
        dy = ylength / nygrid.

    For a given axis, if

        i0 = floor(g)
        t = g - i0,

    then the particle contributes to cells ``i0`` and ``i0 + 1`` with
    dimensionless weights ``1 - t`` and ``t``, respectively.

    For non-periodic axes, stencil contributions falling outside the
    grid domain are discarded. For periodic axes, stencil indices are
    wrapped onto the opposite side of the grid.
    """
    npart = len(x)

    idx = nxgrid / xlength
    idy = nygrid / ylength

    norm = idx * idy

    fgrid = np.zeros((nxgrid, nygrid), dtype=dtype)

    for i in range(npart):

        gx = (x[i] - xmin) * idx - 0.5
        gy = (y[i] - ymin) * idy - 0.5

        ix0 = int(np.floor(gx))
        iy0 = int(np.floor(gy))

        tx = gx - ix0
        ty = gy - iy0

        ix1 = ix0 + 1
        iy1 = iy0 + 1

        wx0 = 1.0 - tx
        wx1 = tx

        wy0 = 1.0 - ty
        wy1 = ty

        if periodx:
            if ix0 < 0:
                ix0 += nxgrid
            elif ix0 >= nxgrid:
                ix0 -= nxgrid

            if ix1 < 0:
                ix1 += nxgrid
            elif ix1 >= nxgrid:
                ix1 -= nxgrid

        if periody:
            if iy0 < 0:
                iy0 += nygrid
            elif iy0 >= nygrid:
                iy0 -= nygrid

            if iy1 < 0:
                iy1 += nygrid
            elif iy1 >= nygrid:
                iy1 -= nygrid

        fp = f[i] * norm

        if (0 <= ix0 < nxgrid) and (0 <= iy0 < nygrid):
            fgrid[ix0, iy0] += fp * wx0 * wy0

        if (0 <= ix0 < nxgrid) and (0 <= iy1 < nygrid):
            fgrid[ix0, iy1] += fp * wx0 * wy1

        if (0 <= ix1 < nxgrid) and (0 <= iy0 < nygrid):
            fgrid[ix1, iy0] += fp * wx1 * wy0

        if (0 <= ix1 < nxgrid) and (0 <= iy1 < nygrid):
            fgrid[ix1, iy1] += fp * wx1 * wy1

    return fgrid



@njit
def part2grid_cic_2d_unit(
    x: np.ndarray,
    y: np.ndarray,
    xlength: float,
    ylength: float,
    xmin: float,
    ymin: float,
    nxgrid: int,
    nygrid: int,
    periodx: bool,
    periody: bool,
    dtype: np.dtype = np.float64
) -> np.ndarray:
    """
    Assign particle values to a two-dimensional grid using the
    cloud-in-cell (CIC) assignment scheme.

    Each particle contributes to the two nearest cell centres along each
    coordinate axis, giving a total stencil of four neighbouring grid
    cells. The contribution to each cell is weighted linearly according
    to the particle position relative to the neighbouring cell centres.

    The deposited field is normalized by the grid-cell area. When the
    complete stencil lies within the domain, or periodic boundaries are
    applied, the weights associated with each particle sum to

        1 / (dx * dy).

    Parameters
    ----------
    x, y : ndarray
        One-dimensional arrays containing the Cartesian coordinates of
        the particles. Both arrays must have the same length.
    xlength, ylength : float
        Physical lengths of the grid domain along the x and y axes.
    xmin, ymin : float
        Minimum coordinates of the grid domain along the x and y axes.
    nxgrid, nygrid : int
        Number of grid cells along the x and y axes.
    periodx, periody : bool
        Whether periodic boundary conditions are applied along the x
        and y axes, respectively.
    dtype : numpy dtype, optional
        Data type used for the output grid. Default is ``np.float64``.

    Returns
    -------
    fgrid : ndarray
        Two-dimensional array of shape ``(nxgrid, nygrid)`` containing
        the CIC-assigned field, with data type given by ``dtype``.

    Notes
    -----
    The grid is treated as cell-centred. The dimensionless particle
    coordinate relative to the grid-cell centres is

        gx = (x - xmin) / dx - 0.5
        gy = (y - ymin) / dy - 0.5,

    where

        dx = xlength / nxgrid
        dy = ylength / nygrid.

    For a given axis, if

        i0 = floor(g)
        t = g - i0,

    then the particle contributes to cells ``i0`` and ``i0 + 1`` with
    dimensionless weights ``1 - t`` and ``t``, respectively.

    For non-periodic axes, stencil contributions falling outside the
    grid domain are discarded. For periodic axes, stencil indices are
    wrapped onto the opposite side of the grid.
    """
    npart = len(x)

    idx = nxgrid / xlength
    idy = nygrid / ylength

    norm = idx * idy

    fgrid = np.zeros((nxgrid, nygrid), dtype=dtype)

    for i in range(npart):

        gx = (x[i] - xmin) * idx - 0.5
        gy = (y[i] - ymin) * idy - 0.5

        ix0 = int(np.floor(gx))
        iy0 = int(np.floor(gy))

        tx = gx - ix0
        ty = gy - iy0

        ix1 = ix0 + 1
        iy1 = iy0 + 1

        wx0 = 1.0 - tx
        wx1 = tx

        wy0 = 1.0 - ty
        wy1 = ty

        if periodx:
            if ix0 < 0:
                ix0 += nxgrid
            elif ix0 >= nxgrid:
                ix0 -= nxgrid

            if ix1 < 0:
                ix1 += nxgrid
            elif ix1 >= nxgrid:
                ix1 -= nxgrid

        if periody:
            if iy0 < 0:
                iy0 += nygrid
            elif iy0 >= nygrid:
                iy0 -= nygrid

            if iy1 < 0:
                iy1 += nygrid
            elif iy1 >= nygrid:
                iy1 -= nygrid

        fp = norm

        if (0 <= ix0 < nxgrid) and (0 <= iy0 < nygrid):
            fgrid[ix0, iy0] += fp * wx0 * wy0

        if (0 <= ix0 < nxgrid) and (0 <= iy1 < nygrid):
            fgrid[ix0, iy1] += fp * wx0 * wy1

        if (0 <= ix1 < nxgrid) and (0 <= iy0 < nygrid):
            fgrid[ix1, iy0] += fp * wx1 * wy0

        if (0 <= ix1 < nxgrid) and (0 <= iy1 < nygrid):
            fgrid[ix1, iy1] += fp * wx1 * wy1

    return fgrid


@njit
def part2grid_tsc_2d(
    x: np.ndarray,
    y: np.ndarray,
    f: np.ndarray,
    xlength: float,
    ylength: float,
    xmin: float,
    ymin: float,
    nxgrid: int,
    nygrid: int,
    periodx: bool,
    periody: bool,
    dtype: np.dtype = np.float64
) -> np.ndarray:
    """
    Assign particle values to a two-dimensional grid using the
    triangular-shaped-cloud (TSC) assignment scheme.

    Each particle contributes to three neighbouring cell centres along
    each coordinate axis, giving a total stencil of nine grid cells. The
    contribution to each cell is determined by the quadratic TSC kernel
    evaluated from the particle displacement relative to the centre of
    its containing grid cell.

    The deposited field is normalized by the grid-cell area. When the
    complete stencil lies within the domain, or periodic boundaries are
    applied, the weights associated with each particle sum to

        1 / (dx * dy).

    Parameters
    ----------
    x, y : ndarray
        One-dimensional arrays containing the Cartesian coordinates of
        the particles. Both arrays must have the same length.
    f : ndarray
        One-dimensional array containing the value associated with each
        particle. Must have the same length as ``x`` and ``y``.
    xlength, ylength : float
        Physical lengths of the grid domain along the x and y axes.
    xmin, ymin : float
        Minimum coordinates of the grid domain along the x and y axes.
    nxgrid, nygrid : int
        Number of grid cells along the x and y axes.
    periodx, periody : bool
        Whether periodic boundary conditions are applied along the x
        and y axes, respectively.
    dtype : numpy dtype, optional
        Data type used for the output grid. Default is ``np.float64``.

    Returns
    -------
    fgrid : ndarray
        Two-dimensional array of shape ``(nxgrid, nygrid)`` containing
        the TSC-assigned field, with data type given by ``dtype``.

    Notes
    -----
    For each axis, the particle position is expressed relative to the
    centre of its containing grid cell. For the x direction,

        gx = (x - xmin) / dx
        ix = floor(gx)
        sx = gx - ix - 0.5,

    where

        dx = xlength / nxgrid,

    and ``sx`` lies in the interval [-0.5, 0.5).

    The particle contributes to cells ``ix - 1``, ``ix`` and ``ix + 1``
    with dimensionless weights

        w_- = 0.5 * (0.5 - sx)**2
        w_0 = 0.75 - sx**2
        w_+ = 0.5 * (0.5 + sx)**2.

    The same construction is applied independently along the y axis,
    and the two-dimensional assignment weight is the product of the
    corresponding one-dimensional weights.

    For non-periodic axes, stencil contributions falling outside the
    grid domain are discarded. For periodic axes, stencil indices are
    wrapped onto the opposite side of the grid.
    """
    npart = len(x)

    idx = nxgrid / xlength
    idy = nygrid / ylength

    norm = idx * idy

    fgrid = np.zeros((nxgrid, nygrid), dtype=dtype)

    for i in range(npart):

        # Particle position in grid-cell coordinates.
        gx = (x[i] - xmin) * idx
        gy = (y[i] - ymin) * idy

        ix0 = int(np.floor(gx))
        iy0 = int(np.floor(gy))

        # Offset from the centre of the containing cell.
        sx = gx - ix0 - 0.5
        sy = gy - iy0 - 0.5

        # Three neighbouring grid cells.
        ixm = ix0 - 1
        ixp = ix0 + 1

        iym = iy0 - 1
        iyp = iy0 + 1

        # TSC weights along x.
        wxm = 0.5 * (0.5 - sx) * (0.5 - sx)
        wx0 = 0.75 - sx * sx
        wxp = 0.5 * (0.5 + sx) * (0.5 + sx)

        # TSC weights along y.
        wym = 0.5 * (0.5 - sy) * (0.5 - sy)
        wy0 = 0.75 - sy * sy
        wyp = 0.5 * (0.5 + sy) * (0.5 + sy)

        # Periodic wrapping.
        if periodx:
            if ixm < 0:
                ixm += nxgrid
            elif ixm >= nxgrid:
                ixm -= nxgrid

            if ix0 < 0:
                ix0 += nxgrid
            elif ix0 >= nxgrid:
                ix0 -= nxgrid

            if ixp < 0:
                ixp += nxgrid
            elif ixp >= nxgrid:
                ixp -= nxgrid

        if periody:
            if iym < 0:
                iym += nygrid
            elif iym >= nygrid:
                iym -= nygrid

            if iy0 < 0:
                iy0 += nygrid
            elif iy0 >= nygrid:
                iy0 -= nygrid

            if iyp < 0:
                iyp += nygrid
            elif iyp >= nygrid:
                iyp -= nygrid

        fp = f[i] * norm

        # x = ixm
        if (0 <= ixm < nxgrid) and (0 <= iym < nygrid):
            fgrid[ixm, iym] += fp * wxm * wym

        if (0 <= ixm < nxgrid) and (0 <= iy0 < nygrid):
            fgrid[ixm, iy0] += fp * wxm * wy0

        if (0 <= ixm < nxgrid) and (0 <= iyp < nygrid):
            fgrid[ixm, iyp] += fp * wxm * wyp

        # x = ix0
        if (0 <= ix0 < nxgrid) and (0 <= iym < nygrid):
            fgrid[ix0, iym] += fp * wx0 * wym

        if (0 <= ix0 < nxgrid) and (0 <= iy0 < nygrid):
            fgrid[ix0, iy0] += fp * wx0 * wy0

        if (0 <= ix0 < nxgrid) and (0 <= iyp < nygrid):
            fgrid[ix0, iyp] += fp * wx0 * wyp

        # x = ixp
        if (0 <= ixp < nxgrid) and (0 <= iym < nygrid):
            fgrid[ixp, iym] += fp * wxp * wym

        if (0 <= ixp < nxgrid) and (0 <= iy0 < nygrid):
            fgrid[ixp, iy0] += fp * wxp * wy0

        if (0 <= ixp < nxgrid) and (0 <= iyp < nygrid):
            fgrid[ixp, iyp] += fp * wxp * wyp

    return fgrid



@njit
def part2grid_tsc_2d_unit(
    x: np.ndarray,
    y: np.ndarray,
    xlength: float,
    ylength: float,
    xmin: float,
    ymin: float,
    nxgrid: int,
    nygrid: int,
    periodx: bool,
    periody: bool,
    dtype: np.dtype = np.float64
) -> np.ndarray:
    """
    Assign particle values to a two-dimensional grid using the
    triangular-shaped-cloud (TSC) assignment scheme.

    Each particle contributes to three neighbouring cell centres along
    each coordinate axis, giving a total stencil of nine grid cells. The
    contribution to each cell is determined by the quadratic TSC kernel
    evaluated from the particle displacement relative to the centre of
    its containing grid cell.

    The deposited field is normalized by the grid-cell area. When the
    complete stencil lies within the domain, or periodic boundaries are
    applied, the weights associated with each particle sum to

        1 / (dx * dy).

    Parameters
    ----------
    x, y : ndarray
        One-dimensional arrays containing the Cartesian coordinates of
        the particles. Both arrays must have the same length.
    xlength, ylength : float
        Physical lengths of the grid domain along the x and y axes.
    xmin, ymin : float
        Minimum coordinates of the grid domain along the x and y axes.
    nxgrid, nygrid : int
        Number of grid cells along the x and y axes.
    periodx, periody : bool
        Whether periodic boundary conditions are applied along the x
        and y axes, respectively.
    dtype : numpy dtype, optional
        Data type used for the output grid. Default is ``np.float64``.

    Returns
    -------
    fgrid : ndarray
        Two-dimensional array of shape ``(nxgrid, nygrid)`` containing
        the TSC-assigned field, with data type given by ``dtype``.

    Notes
    -----
    For each axis, the particle position is expressed relative to the
    centre of its containing grid cell. For the x direction,

        gx = (x - xmin) / dx
        ix = floor(gx)
        sx = gx - ix - 0.5,

    where

        dx = xlength / nxgrid,

    and ``sx`` lies in the interval [-0.5, 0.5).

    The particle contributes to cells ``ix - 1``, ``ix`` and ``ix + 1``
    with dimensionless weights

        w_- = 0.5 * (0.5 - sx)**2
        w_0 = 0.75 - sx**2
        w_+ = 0.5 * (0.5 + sx)**2.

    The same construction is applied independently along the y axis,
    and the two-dimensional assignment weight is the product of the
    corresponding one-dimensional weights.

    For non-periodic axes, stencil contributions falling outside the
    grid domain are discarded. For periodic axes, stencil indices are
    wrapped onto the opposite side of the grid.
    """
    npart = len(x)

    idx = nxgrid / xlength
    idy = nygrid / ylength

    norm = idx * idy

    fgrid = np.zeros((nxgrid, nygrid), dtype=dtype)

    for i in range(npart):

        # Particle position in grid-cell coordinates.
        gx = (x[i] - xmin) * idx
        gy = (y[i] - ymin) * idy

        ix0 = int(np.floor(gx))
        iy0 = int(np.floor(gy))

        # Offset from the centre of the containing cell.
        sx = gx - ix0 - 0.5
        sy = gy - iy0 - 0.5

        # Three neighbouring grid cells.
        ixm = ix0 - 1
        ixp = ix0 + 1

        iym = iy0 - 1
        iyp = iy0 + 1

        # TSC weights along x.
        wxm = 0.5 * (0.5 - sx) * (0.5 - sx)
        wx0 = 0.75 - sx * sx
        wxp = 0.5 * (0.5 + sx) * (0.5 + sx)

        # TSC weights along y.
        wym = 0.5 * (0.5 - sy) * (0.5 - sy)
        wy0 = 0.75 - sy * sy
        wyp = 0.5 * (0.5 + sy) * (0.5 + sy)

        # Periodic wrapping.
        if periodx:
            if ixm < 0:
                ixm += nxgrid
            elif ixm >= nxgrid:
                ixm -= nxgrid

            if ix0 < 0:
                ix0 += nxgrid
            elif ix0 >= nxgrid:
                ix0 -= nxgrid

            if ixp < 0:
                ixp += nxgrid
            elif ixp >= nxgrid:
                ixp -= nxgrid

        if periody:
            if iym < 0:
                iym += nygrid
            elif iym >= nygrid:
                iym -= nygrid

            if iy0 < 0:
                iy0 += nygrid
            elif iy0 >= nygrid:
                iy0 -= nygrid

            if iyp < 0:
                iyp += nygrid
            elif iyp >= nygrid:
                iyp -= nygrid

        fp = norm

        # x = ixm
        if (0 <= ixm < nxgrid) and (0 <= iym < nygrid):
            fgrid[ixm, iym] += fp * wxm * wym

        if (0 <= ixm < nxgrid) and (0 <= iy0 < nygrid):
            fgrid[ixm, iy0] += fp * wxm * wy0

        if (0 <= ixm < nxgrid) and (0 <= iyp < nygrid):
            fgrid[ixm, iyp] += fp * wxm * wyp

        # x = ix0
        if (0 <= ix0 < nxgrid) and (0 <= iym < nygrid):
            fgrid[ix0, iym] += fp * wx0 * wym

        if (0 <= ix0 < nxgrid) and (0 <= iy0 < nygrid):
            fgrid[ix0, iy0] += fp * wx0 * wy0

        if (0 <= ix0 < nxgrid) and (0 <= iyp < nygrid):
            fgrid[ix0, iyp] += fp * wx0 * wyp

        # x = ixp
        if (0 <= ixp < nxgrid) and (0 <= iym < nygrid):
            fgrid[ixp, iym] += fp * wxp * wym

        if (0 <= ixp < nxgrid) and (0 <= iy0 < nygrid):
            fgrid[ixp, iy0] += fp * wxp * wy0

        if (0 <= ixp < nxgrid) and (0 <= iyp < nygrid):
            fgrid[ixp, iyp] += fp * wxp * wyp

    return fgrid


@njit
def part2grid_pcs_2d(
    x: np.ndarray,
    y: np.ndarray,
    f: np.ndarray,
    xlength: float,
    ylength: float,
    xmin: float,
    ymin: float,
    nxgrid: int,
    nygrid: int,
    periodx: bool,
    periody: bool,
    dtype: np.dtype = np.float64
) -> np.ndarray:
    """
    Assign particle values to a two-dimensional grid using the
    piecewise-cubic-spline (PCS) assignment scheme.

    Each particle contributes to four neighbouring cell centres along
    each coordinate axis, giving a total stencil of 16 grid cells. The
    contribution to each cell is determined by a cubic B-spline kernel
    evaluated from the particle position relative to the neighbouring
    grid-cell centres.

    The deposited field is normalized by the grid-cell area. When the
    complete stencil lies within the domain, or periodic boundaries are
    applied, the weights associated with each particle sum to

        1 / (dx * dy).

    Parameters
    ----------
    x, y : ndarray
        One-dimensional arrays containing the Cartesian coordinates of
        the particles. Both arrays must have the same length.
    f : ndarray
        One-dimensional array containing the value associated with each
        particle. Must have the same length as ``x`` and ``y``.
    xlength, ylength : float
        Physical lengths of the grid domain along the x and y axes.
    xmin, ymin : float
        Minimum coordinates of the grid domain along the x and y axes.
    nxgrid, nygrid : int
        Number of grid cells along the x and y axes.
    periodx, periody : bool
        Whether periodic boundary conditions are applied along the x
        and y axes, respectively.
    dtype : numpy dtype, optional
        Data type used for the output grid. Default is ``np.float64``.

    Returns
    -------
    fgrid : ndarray
        Two-dimensional array of shape ``(nxgrid, nygrid)`` containing
        the PCS-assigned field, with data type given by ``dtype``.

    Notes
    -----
    The grid is treated as cell-centred. For the x direction, define

        gx = (x - xmin) / dx - 0.5
        ix = floor(gx)
        tx = gx - ix,

    where

        dx = xlength / nxgrid,

    and ``tx`` lies in the interval [0, 1).

    The particle contributes to cells ``ix - 1``, ``ix``, ``ix + 1``
    and ``ix + 2`` with dimensionless cubic B-spline weights

        w0 = (1 - tx)**3 / 6

        w1 = (4 - 6*tx**2 + 3*tx**3) / 6

        w2 = (1 + 3*tx + 3*tx**2 - 3*tx**3) / 6

        w3 = tx**3 / 6.

    The same construction is applied independently along the y axis,
    and the two-dimensional assignment weight is the product of the
    corresponding one-dimensional weights.

    For non-periodic axes, stencil contributions falling outside the
    grid domain are discarded. For periodic axes, stencil indices are
    wrapped onto the opposite side of the grid.
    """
    npart = len(x)

    idx = nxgrid / xlength
    idy = nygrid / ylength

    norm = idx * idy

    fgrid = np.zeros((nxgrid, nygrid), dtype=dtype)

    for i in range(npart):

        gx = (x[i] - xmin) * idx - 0.5
        gy = (y[i] - ymin) * idy - 0.5

        ix1 = int(np.floor(gx))
        iy1 = int(np.floor(gy))

        tx = gx - ix1
        ty = gy - iy1

        ix0 = ix1 - 1
        ix2 = ix1 + 1
        ix3 = ix1 + 2

        iy0 = iy1 - 1
        iy2 = iy1 + 1
        iy3 = iy1 + 2

        tx2 = tx * tx
        tx3 = tx2 * tx

        ty2 = ty * ty
        ty3 = ty2 * ty

        omt_x = 1.0 - tx
        omt_y = 1.0 - ty

        wx0 = omt_x * omt_x * omt_x / 6.0
        wx1 = (4.0 - 6.0 * tx2 + 3.0 * tx3) / 6.0
        wx2 = (1.0 + 3.0 * tx + 3.0 * tx2 - 3.0 * tx3) / 6.0
        wx3 = tx3 / 6.0

        wy0 = omt_y * omt_y * omt_y / 6.0
        wy1 = (4.0 - 6.0 * ty2 + 3.0 * ty3) / 6.0
        wy2 = (1.0 + 3.0 * ty + 3.0 * ty2 - 3.0 * ty3) / 6.0
        wy3 = ty3 / 6.0

        if periodx:
            if ix0 < 0:
                ix0 += nxgrid
            elif ix0 >= nxgrid:
                ix0 -= nxgrid

            if ix1 < 0:
                ix1 += nxgrid
            elif ix1 >= nxgrid:
                ix1 -= nxgrid

            if ix2 < 0:
                ix2 += nxgrid
            elif ix2 >= nxgrid:
                ix2 -= nxgrid

            if ix3 < 0:
                ix3 += nxgrid
            elif ix3 >= nxgrid:
                ix3 -= nxgrid

        if periody:
            if iy0 < 0:
                iy0 += nygrid
            elif iy0 >= nygrid:
                iy0 -= nygrid

            if iy1 < 0:
                iy1 += nygrid
            elif iy1 >= nygrid:
                iy1 -= nygrid

            if iy2 < 0:
                iy2 += nygrid
            elif iy2 >= nygrid:
                iy2 -= nygrid

            if iy3 < 0:
                iy3 += nygrid
            elif iy3 >= nygrid:
                iy3 -= nygrid

        fp = f[i] * norm

        # ix0
        if (0 <= ix0 < nxgrid) and (0 <= iy0 < nygrid):
            fgrid[ix0, iy0] += fp * wx0 * wy0
        if (0 <= ix0 < nxgrid) and (0 <= iy1 < nygrid):
            fgrid[ix0, iy1] += fp * wx0 * wy1
        if (0 <= ix0 < nxgrid) and (0 <= iy2 < nygrid):
            fgrid[ix0, iy2] += fp * wx0 * wy2
        if (0 <= ix0 < nxgrid) and (0 <= iy3 < nygrid):
            fgrid[ix0, iy3] += fp * wx0 * wy3

        # ix1
        if (0 <= ix1 < nxgrid) and (0 <= iy0 < nygrid):
            fgrid[ix1, iy0] += fp * wx1 * wy0
        if (0 <= ix1 < nxgrid) and (0 <= iy1 < nygrid):
            fgrid[ix1, iy1] += fp * wx1 * wy1
        if (0 <= ix1 < nxgrid) and (0 <= iy2 < nygrid):
            fgrid[ix1, iy2] += fp * wx1 * wy2
        if (0 <= ix1 < nxgrid) and (0 <= iy3 < nygrid):
            fgrid[ix1, iy3] += fp * wx1 * wy3

        # ix2
        if (0 <= ix2 < nxgrid) and (0 <= iy0 < nygrid):
            fgrid[ix2, iy0] += fp * wx2 * wy0
        if (0 <= ix2 < nxgrid) and (0 <= iy1 < nygrid):
            fgrid[ix2, iy1] += fp * wx2 * wy1
        if (0 <= ix2 < nxgrid) and (0 <= iy2 < nygrid):
            fgrid[ix2, iy2] += fp * wx2 * wy2
        if (0 <= ix2 < nxgrid) and (0 <= iy3 < nygrid):
            fgrid[ix2, iy3] += fp * wx2 * wy3

        # ix3
        if (0 <= ix3 < nxgrid) and (0 <= iy0 < nygrid):
            fgrid[ix3, iy0] += fp * wx3 * wy0
        if (0 <= ix3 < nxgrid) and (0 <= iy1 < nygrid):
            fgrid[ix3, iy1] += fp * wx3 * wy1
        if (0 <= ix3 < nxgrid) and (0 <= iy2 < nygrid):
            fgrid[ix3, iy2] += fp * wx3 * wy2
        if (0 <= ix3 < nxgrid) and (0 <= iy3 < nygrid):
            fgrid[ix3, iy3] += fp * wx3 * wy3

    return fgrid


@njit
def part2grid_pcs_2d_unit(
    x: np.ndarray,
    y: np.ndarray,
    xlength: float,
    ylength: float,
    xmin: float,
    ymin: float,
    nxgrid: int,
    nygrid: int,
    periodx: bool,
    periody: bool,
    dtype: np.dtype = np.float64
) -> np.ndarray:
    """
    Assign particle values to a two-dimensional grid using the
    piecewise-cubic-spline (PCS) assignment scheme.

    Each particle contributes to four neighbouring cell centres along
    each coordinate axis, giving a total stencil of 16 grid cells. The
    contribution to each cell is determined by a cubic B-spline kernel
    evaluated from the particle position relative to the neighbouring
    grid-cell centres.

    The deposited field is normalized by the grid-cell area. When the
    complete stencil lies within the domain, or periodic boundaries are
    applied, the weights associated with each particle sum to

        1 / (dx * dy).

    Parameters
    ----------
    x, y : ndarray
        One-dimensional arrays containing the Cartesian coordinates of
        the particles. Both arrays must have the same length.
    xlength, ylength : float
        Physical lengths of the grid domain along the x and y axes.
    xmin, ymin : float
        Minimum coordinates of the grid domain along the x and y axes.
    nxgrid, nygrid : int
        Number of grid cells along the x and y axes.
    periodx, periody : bool
        Whether periodic boundary conditions are applied along the x
        and y axes, respectively.
    dtype : numpy dtype, optional
        Data type used for the output grid. Default is ``np.float64``.

    Returns
    -------
    fgrid : ndarray
        Two-dimensional array of shape ``(nxgrid, nygrid)`` containing
        the PCS-assigned field, with data type given by ``dtype``.

    Notes
    -----
    The grid is treated as cell-centred. For the x direction, define

        gx = (x - xmin) / dx - 0.5
        ix = floor(gx)
        tx = gx - ix,

    where

        dx = xlength / nxgrid,

    and ``tx`` lies in the interval [0, 1).

    The particle contributes to cells ``ix - 1``, ``ix``, ``ix + 1``
    and ``ix + 2`` with dimensionless cubic B-spline weights

        w0 = (1 - tx)**3 / 6

        w1 = (4 - 6*tx**2 + 3*tx**3) / 6

        w2 = (1 + 3*tx + 3*tx**2 - 3*tx**3) / 6

        w3 = tx**3 / 6.

    The same construction is applied independently along the y axis,
    and the two-dimensional assignment weight is the product of the
    corresponding one-dimensional weights.

    For non-periodic axes, stencil contributions falling outside the
    grid domain are discarded. For periodic axes, stencil indices are
    wrapped onto the opposite side of the grid.
    """
    npart = len(x)

    idx = nxgrid / xlength
    idy = nygrid / ylength

    norm = idx * idy

    fgrid = np.zeros((nxgrid, nygrid), dtype=dtype)

    for i in range(npart):

        gx = (x[i] - xmin) * idx - 0.5
        gy = (y[i] - ymin) * idy - 0.5

        ix1 = int(np.floor(gx))
        iy1 = int(np.floor(gy))

        tx = gx - ix1
        ty = gy - iy1

        ix0 = ix1 - 1
        ix2 = ix1 + 1
        ix3 = ix1 + 2

        iy0 = iy1 - 1
        iy2 = iy1 + 1
        iy3 = iy1 + 2

        tx2 = tx * tx
        tx3 = tx2 * tx

        ty2 = ty * ty
        ty3 = ty2 * ty

        omt_x = 1.0 - tx
        omt_y = 1.0 - ty

        wx0 = omt_x * omt_x * omt_x / 6.0
        wx1 = (4.0 - 6.0 * tx2 + 3.0 * tx3) / 6.0
        wx2 = (1.0 + 3.0 * tx + 3.0 * tx2 - 3.0 * tx3) / 6.0
        wx3 = tx3 / 6.0

        wy0 = omt_y * omt_y * omt_y / 6.0
        wy1 = (4.0 - 6.0 * ty2 + 3.0 * ty3) / 6.0
        wy2 = (1.0 + 3.0 * ty + 3.0 * ty2 - 3.0 * ty3) / 6.0
        wy3 = ty3 / 6.0

        if periodx:
            if ix0 < 0:
                ix0 += nxgrid
            elif ix0 >= nxgrid:
                ix0 -= nxgrid

            if ix1 < 0:
                ix1 += nxgrid
            elif ix1 >= nxgrid:
                ix1 -= nxgrid

            if ix2 < 0:
                ix2 += nxgrid
            elif ix2 >= nxgrid:
                ix2 -= nxgrid

            if ix3 < 0:
                ix3 += nxgrid
            elif ix3 >= nxgrid:
                ix3 -= nxgrid

        if periody:
            if iy0 < 0:
                iy0 += nygrid
            elif iy0 >= nygrid:
                iy0 -= nygrid

            if iy1 < 0:
                iy1 += nygrid
            elif iy1 >= nygrid:
                iy1 -= nygrid

            if iy2 < 0:
                iy2 += nygrid
            elif iy2 >= nygrid:
                iy2 -= nygrid

            if iy3 < 0:
                iy3 += nygrid
            elif iy3 >= nygrid:
                iy3 -= nygrid

        fp = norm

        # ix0
        if (0 <= ix0 < nxgrid) and (0 <= iy0 < nygrid):
            fgrid[ix0, iy0] += fp * wx0 * wy0
        if (0 <= ix0 < nxgrid) and (0 <= iy1 < nygrid):
            fgrid[ix0, iy1] += fp * wx0 * wy1
        if (0 <= ix0 < nxgrid) and (0 <= iy2 < nygrid):
            fgrid[ix0, iy2] += fp * wx0 * wy2
        if (0 <= ix0 < nxgrid) and (0 <= iy3 < nygrid):
            fgrid[ix0, iy3] += fp * wx0 * wy3

        # ix1
        if (0 <= ix1 < nxgrid) and (0 <= iy0 < nygrid):
            fgrid[ix1, iy0] += fp * wx1 * wy0
        if (0 <= ix1 < nxgrid) and (0 <= iy1 < nygrid):
            fgrid[ix1, iy1] += fp * wx1 * wy1
        if (0 <= ix1 < nxgrid) and (0 <= iy2 < nygrid):
            fgrid[ix1, iy2] += fp * wx1 * wy2
        if (0 <= ix1 < nxgrid) and (0 <= iy3 < nygrid):
            fgrid[ix1, iy3] += fp * wx1 * wy3

        # ix2
        if (0 <= ix2 < nxgrid) and (0 <= iy0 < nygrid):
            fgrid[ix2, iy0] += fp * wx2 * wy0
        if (0 <= ix2 < nxgrid) and (0 <= iy1 < nygrid):
            fgrid[ix2, iy1] += fp * wx2 * wy1
        if (0 <= ix2 < nxgrid) and (0 <= iy2 < nygrid):
            fgrid[ix2, iy2] += fp * wx2 * wy2
        if (0 <= ix2 < nxgrid) and (0 <= iy3 < nygrid):
            fgrid[ix2, iy3] += fp * wx2 * wy3

        # ix3
        if (0 <= ix3 < nxgrid) and (0 <= iy0 < nygrid):
            fgrid[ix3, iy0] += fp * wx3 * wy0
        if (0 <= ix3 < nxgrid) and (0 <= iy1 < nygrid):
            fgrid[ix3, iy1] += fp * wx3 * wy1
        if (0 <= ix3 < nxgrid) and (0 <= iy2 < nygrid):
            fgrid[ix3, iy2] += fp * wx3 * wy2
        if (0 <= ix3 < nxgrid) and (0 <= iy3 < nygrid):
            fgrid[ix3, iy3] += fp * wx3 * wy3

    return fgrid


def part2grid2D(
    x: np.ndarray,
    y: np.ndarray,
    boxsize: Union[float, List[float]],
    ngrid: Union[int, List[int]],
    f: np.ndarray = None,
    method: str = "TSC",
    periodic: Union[bool, List[bool]] = True,
    origin: Union[float, List[float]] = 0.0,
    dtype: np.dtype = np.float64
) -> np.ndarray:
    """
    Returns the density contrast for the nearest grid point grid assignment.

    Parameters
    ----------
    x : array
        X coordinates of the particle.
    y : array
        Y coordinates of the particle.
    boxsize : float or list
        Box size.
    ngrid : int
        Grid divisions across one axis.
    f : array
        Value of each particle to be assigned to the grid.
    method : str, optional
        Grid assignment scheme, either 'NGP', 'CIC', 'TSC' or 'PCS'.
    periodic : bool or list, optional
        Assign particles with periodic boundaries.
    origin : float or list, optional
        Origin.
    dtype : numpy dtype, optional
        Data type used for the output grid. Default is ``np.float64``.
    
    Returns
    -------
    fgrid : array
        Grid assigned values.
    """
    if np.isscalar(boxsize):
        xlength, ylength = boxsize, boxsize
    else:
        xlength, ylength = boxsize[0], boxsize[1]
    if np.isscalar(origin):
        xmin = origin
        ymin = origin
    else:
        xmin, ymin = origin[0], origin[1]
    if np.isscalar(ngrid):
        nxgrid, nygrid = int(ngrid), int(ngrid)
    else:
        nxgrid, nygrid = int(ngrid[0]), int(ngrid[1])
    if np.isscalar(periodic):
        periodx = periodic
        periody = periodic
    else:
        periodx, periody = periodic[0], periodic[1]
    if method == "NGP":
        if f is None:
            fgrid = part2grid_ngp_2d_unit(
                x, y, xlength, ylength, xmin, ymin, nxgrid, nygrid, dtype=dtype
            )
        else:
            fgrid = part2grid_ngp_2d(
                x, y, f, xlength, ylength, xmin, ymin, nxgrid, nygrid, dtype=dtype
            )
    elif method == "CIC":
        if f is None:
            fgrid = part2grid_cic_2d_unit(
                x, y, xlength, ylength, xmin, ymin, nxgrid, nygrid, periodx, periody, dtype=dtype
            )
        else:
            fgrid = part2grid_cic_2d(
                x, y, f, xlength, ylength, xmin, ymin, nxgrid, nygrid, periodx, periody, dtype=dtype
            )
    elif method == "TSC":
        if f is None:
            fgrid = part2grid_tsc_2d_unit(
                x, y, xlength, ylength, xmin, ymin, nxgrid, nygrid, periodx, periody, dtype=dtype
            )
        else:
            fgrid = part2grid_tsc_2d(
                x, y, f, xlength, ylength, xmin, ymin, nxgrid, nygrid, periodx, periody, dtype=dtype
            )
    elif method == "PCS":
        if f is None:
            fgrid = part2grid_pcs_2d_unit(
                x, y, xlength, ylength, xmin, ymin, nxgrid, nygrid, periodx, periody, dtype=dtype
            )
        else:
            fgrid = part2grid_pcs_2d(
                x, y, f, xlength, ylength, xmin, ymin, nxgrid, nygrid, periodx, periody, dtype=dtype
            )
    return fgrid


@njit
def part2grid_ngp_3d(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    f: np.ndarray,
    xlength: float,
    ylength: float,
    zlength: float,
    xmin: float,
    ymin: float,
    zmin: float,
    nxgrid: int,
    nygrid: int,
    nzgrid: int,
    dtype: np.dtype = np.float64
) -> np.ndarray:
    """
    Assign particle values to a three-dimensional grid using the
    nearest-grid-point (NGP) assignment scheme.

    Each particle is assigned to the grid cell containing its position.
    The deposited value is normalized by the cell volume, such that a
    particle value ``f[i]`` contributes

        f[i] / (dx * dy * dz)

    to its corresponding grid cell.

    Particles lying outside the specified grid domain are ignored.

    Parameters
    ----------
    x, y, z : ndarray
        One-dimensional arrays containing the Cartesian coordinates of
        the particles. All arrays must have the same length.
    f : ndarray
        One-dimensional array containing the value associated with each
        particle. Must have the same length as ``x``, ``y`` and ``z``.
    xlength, ylength, zlength : float
        Physical lengths of the grid domain along the x, y and z axes.
    xmin, ymin, zmin : float
        Minimum coordinates of the grid domain along the x, y and z axes.
    nxgrid, nygrid, nzgrid : int
        Number of grid cells along the x, y and z axes.
    dtype : numpy dtype, optional
        Data type used for the output grid. Default is ``np.float64``.

    Returns
    -------
    fgrid : ndarray
        Three-dimensional array of shape
        ``(nxgrid, nygrid, nzgrid)`` containing the NGP-assigned field.

    Notes
    -----
    The grid spacings are

        dx = xlength / nxgrid
        dy = ylength / nygrid
        dz = zlength / nzgrid

    and particle positions are mapped to cell indices according to

        ix = floor((x - xmin) / dx)
        iy = floor((y - ymin) / dy)
        iz = floor((z - zmin) / dz)

    No periodic wrapping is applied in this routine.
    """
    npart = len(x)

    idx = nxgrid / xlength
    idy = nygrid / ylength
    idz = nzgrid / zlength

    wngp = idx * idy * idz

    fgrid = np.zeros((nxgrid, nygrid, nzgrid), dtype=dtype)

    for i in range(npart):

        ix = int(np.floor((x[i] - xmin) * idx))
        iy = int(np.floor((y[i] - ymin) * idy))
        iz = int(np.floor((z[i] - zmin) * idz))

        if (
            (ix >= 0)
            and (ix < nxgrid)
            and (iy >= 0)
            and (iy < nygrid)
            and (iz >= 0)
            and (iz < nzgrid)
        ):
            fgrid[ix, iy, iz] += f[i] * wngp

    return fgrid



@njit
def part2grid_ngp_3d_unit(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    xlength: float,
    ylength: float,
    zlength: float,
    xmin: float,
    ymin: float,
    zmin: float,
    nxgrid: int,
    nygrid: int,
    nzgrid: int,
    dtype: np.dtype = np.float64
) -> np.ndarray:
    """
    Assign particle values to a three-dimensional grid using the
    nearest-grid-point (NGP) assignment scheme.

    Each particle is assigned to the grid cell containing its position.
    The deposited value is normalized by the cell volume, such that a
    particle value ``f[i]`` contributes

        f[i] / (dx * dy * dz)

    to its corresponding grid cell.

    Particles lying outside the specified grid domain are ignored.

    Parameters
    ----------
    x, y, z : ndarray
        One-dimensional arrays containing the Cartesian coordinates of
        the particles. All arrays must have the same length.
    xlength, ylength, zlength : float
        Physical lengths of the grid domain along the x, y and z axes.
    xmin, ymin, zmin : float
        Minimum coordinates of the grid domain along the x, y and z axes.
    nxgrid, nygrid, nzgrid : int
        Number of grid cells along the x, y and z axes.
    dtype : numpy dtype, optional
        Data type used for the output grid. Default is ``np.float64``.
        
    Returns
    -------
    fgrid : ndarray
        Three-dimensional array of shape
        ``(nxgrid, nygrid, nzgrid)`` containing the NGP-assigned field.

    Notes
    -----
    The grid spacings are

        dx = xlength / nxgrid
        dy = ylength / nygrid
        dz = zlength / nzgrid

    and particle positions are mapped to cell indices according to

        ix = floor((x - xmin) / dx)
        iy = floor((y - ymin) / dy)
        iz = floor((z - zmin) / dz)

    No periodic wrapping is applied in this routine.
    """
    npart = len(x)

    idx = nxgrid / xlength
    idy = nygrid / ylength
    idz = nzgrid / zlength

    wngp = idx * idy * idz

    fgrid = np.zeros((nxgrid, nygrid, nzgrid), dtype=dtype)

    for i in range(npart):

        ix = int(np.floor((x[i] - xmin) * idx))
        iy = int(np.floor((y[i] - ymin) * idy))
        iz = int(np.floor((z[i] - zmin) * idz))

        if (
            (ix >= 0)
            and (ix < nxgrid)
            and (iy >= 0)
            and (iy < nygrid)
            and (iz >= 0)
            and (iz < nzgrid)
        ):
            fgrid[ix, iy, iz] += wngp

    return fgrid


@njit
def part2grid_cic_3d(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    f: np.ndarray,
    xlength: float,
    ylength: float,
    zlength: float,
    xmin: float,
    ymin: float,
    zmin: float,
    nxgrid: int,
    nygrid: int,
    nzgrid: int,
    periodx: bool,
    periody: bool,
    periodz: bool,
    dtype: np.dtype = np.float64
) -> np.ndarray:
    """
    Assign particle values to a three-dimensional grid using the
    cloud-in-cell (CIC) assignment scheme.

    Each particle contributes to the two nearest cell centres along each
    coordinate axis, giving a total stencil of eight neighbouring grid
    cells in three dimensions. The contribution to each cell is weighted
    linearly according to the particle's position relative to the two
    neighbouring cell centres.

    The deposited field is normalized by the cell volume, such that the
    sum of the CIC weights associated with a particle is

        1 / (dx * dy * dz)

    when the full stencil lies within the grid or periodic boundaries are
    applied.

    Parameters
    ----------
    x, y, z : ndarray
        One-dimensional arrays containing the Cartesian coordinates of
        the particles. All arrays must have the same length.
    f : ndarray
        One-dimensional array containing the value associated with each
        particle. Must have the same length as ``x``, ``y`` and ``z``.
    xlength, ylength, zlength : float
        Physical lengths of the grid domain along the x, y and z axes.
    xmin, ymin, zmin : float
        Minimum coordinates of the grid domain along the x, y and z axes.
    nxgrid, nygrid, nzgrid : int
        Number of grid cells along the x, y and z axes.
    periodx, periody, periodz : bool
        Whether periodic boundary conditions are applied along the
        corresponding axis.
    dtype : numpy dtype, optional
        Data type used for the output grid. Default is ``np.float64``.
    
    Returns
    -------
    fgrid : ndarray
        Three-dimensional array of shape
        ``(nxgrid, nygrid, nzgrid)`` containing the CIC-assigned field.

    Notes
    -----
    The grid is treated as cell-centred. The dimensionless coordinate
    relative to the cell centres is therefore

        gx = (x - xmin) / dx - 0.5

    and similarly for y and z.

    For each axis, if

        i0 = floor(gx)
        t  = gx - i0

    then the particle contributes to grid cells ``i0`` and ``i0 + 1``
    with weights ``1 - t`` and ``t``, respectively.

    For non-periodic axes, contributions falling outside the grid domain
    are discarded. For periodic axes, stencil indices are wrapped onto
    the opposite side of the grid.
    """
    npart = len(x)

    idx = nxgrid / xlength
    idy = nygrid / ylength
    idz = nzgrid / zlength

    norm = idx * idy * idz

    fgrid = np.zeros((nxgrid, nygrid, nzgrid), dtype=dtype)

    for i in range(npart):

        gx = (x[i] - xmin) * idx - 0.5
        gy = (y[i] - ymin) * idy - 0.5
        gz = (z[i] - zmin) * idz - 0.5

        ix0 = int(np.floor(gx))
        iy0 = int(np.floor(gy))
        iz0 = int(np.floor(gz))

        tx = gx - ix0
        ty = gy - iy0
        tz = gz - iz0

        ix1 = ix0 + 1
        iy1 = iy0 + 1
        iz1 = iz0 + 1

        wx0 = 1.0 - tx
        wx1 = tx

        wy0 = 1.0 - ty
        wy1 = ty

        wz0 = 1.0 - tz
        wz1 = tz

        if periodx:
            if ix0 < 0:
                ix0 += nxgrid
            elif ix0 >= nxgrid:
                ix0 -= nxgrid

            if ix1 < 0:
                ix1 += nxgrid
            elif ix1 >= nxgrid:
                ix1 -= nxgrid

        if periody:
            if iy0 < 0:
                iy0 += nygrid
            elif iy0 >= nygrid:
                iy0 -= nygrid

            if iy1 < 0:
                iy1 += nygrid
            elif iy1 >= nygrid:
                iy1 -= nygrid

        if periodz:
            if iz0 < 0:
                iz0 += nzgrid
            elif iz0 >= nzgrid:
                iz0 -= nzgrid

            if iz1 < 0:
                iz1 += nzgrid
            elif iz1 >= nzgrid:
                iz1 -= nzgrid

        fp = f[i] * norm

        if (
            (0 <= ix0 < nxgrid)
            and (0 <= iy0 < nygrid)
            and (0 <= iz0 < nzgrid)
        ):
            fgrid[ix0, iy0, iz0] += fp * wx0 * wy0 * wz0

        if (
            (0 <= ix0 < nxgrid)
            and (0 <= iy0 < nygrid)
            and (0 <= iz1 < nzgrid)
        ):
            fgrid[ix0, iy0, iz1] += fp * wx0 * wy0 * wz1

        if (
            (0 <= ix0 < nxgrid)
            and (0 <= iy1 < nygrid)
            and (0 <= iz0 < nzgrid)
        ):
            fgrid[ix0, iy1, iz0] += fp * wx0 * wy1 * wz0

        if (
            (0 <= ix0 < nxgrid)
            and (0 <= iy1 < nygrid)
            and (0 <= iz1 < nzgrid)
        ):
            fgrid[ix0, iy1, iz1] += fp * wx0 * wy1 * wz1

        if (
            (0 <= ix1 < nxgrid)
            and (0 <= iy0 < nygrid)
            and (0 <= iz0 < nzgrid)
        ):
            fgrid[ix1, iy0, iz0] += fp * wx1 * wy0 * wz0

        if (
            (0 <= ix1 < nxgrid)
            and (0 <= iy0 < nygrid)
            and (0 <= iz1 < nzgrid)
        ):
            fgrid[ix1, iy0, iz1] += fp * wx1 * wy0 * wz1

        if (
            (0 <= ix1 < nxgrid)
            and (0 <= iy1 < nygrid)
            and (0 <= iz0 < nzgrid)
        ):
            fgrid[ix1, iy1, iz0] += fp * wx1 * wy1 * wz0

        if (
            (0 <= ix1 < nxgrid)
            and (0 <= iy1 < nygrid)
            and (0 <= iz1 < nzgrid)
        ):
            fgrid[ix1, iy1, iz1] += fp * wx1 * wy1 * wz1

    return fgrid


@njit
def part2grid_cic_3d_unit(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    xlength: float,
    ylength: float,
    zlength: float,
    xmin: float,
    ymin: float,
    zmin: float,
    nxgrid: int,
    nygrid: int,
    nzgrid: int,
    periodx: bool,
    periody: bool,
    periodz: bool,
    dtype: np.dtype = np.float64
) -> np.ndarray:
    """
    Assign particle values to a three-dimensional grid using the
    cloud-in-cell (CIC) assignment scheme.

    Each particle contributes to the two nearest cell centres along each
    coordinate axis, giving a total stencil of eight neighbouring grid
    cells in three dimensions. The contribution to each cell is weighted
    linearly according to the particle's position relative to the two
    neighbouring cell centres.

    The deposited field is normalized by the cell volume, such that the
    sum of the CIC weights associated with a particle is

        1 / (dx * dy * dz)

    when the full stencil lies within the grid or periodic boundaries are
    applied.

    Parameters
    ----------
    x, y, z : ndarray
        One-dimensional arrays containing the Cartesian coordinates of
        the particles. All arrays must have the same length.
    f : ndarray
        One-dimensional array containing the value associated with each
        particle. Must have the same length as ``x``, ``y`` and ``z``.
    xlength, ylength, zlength : float
        Physical lengths of the grid domain along the x, y and z axes.
    xmin, ymin, zmin : float
        Minimum coordinates of the grid domain along the x, y and z axes.
    nxgrid, nygrid, nzgrid : int
        Number of grid cells along the x, y and z axes.
    periodx, periody, periodz : bool
        Whether periodic boundary conditions are applied along the
        corresponding axis.
    dtype : numpy dtype, optional
        Data type used for the output grid. Default is ``np.float64``.
            
    Returns
    -------
    fgrid : ndarray
        Three-dimensional array of shape
        ``(nxgrid, nygrid, nzgrid)`` containing the CIC-assigned field.

    Notes
    -----
    The grid is treated as cell-centred. The dimensionless coordinate
    relative to the cell centres is therefore

        gx = (x - xmin) / dx - 0.5

    and similarly for y and z.

    For each axis, if

        i0 = floor(gx)
        t  = gx - i0

    then the particle contributes to grid cells ``i0`` and ``i0 + 1``
    with weights ``1 - t`` and ``t``, respectively.

    For non-periodic axes, contributions falling outside the grid domain
    are discarded. For periodic axes, stencil indices are wrapped onto
    the opposite side of the grid.
    """
    npart = len(x)

    idx = nxgrid / xlength
    idy = nygrid / ylength
    idz = nzgrid / zlength

    norm = idx * idy * idz

    fgrid = np.zeros((nxgrid, nygrid, nzgrid), dtype=dtype)

    for i in range(npart):

        gx = (x[i] - xmin) * idx - 0.5
        gy = (y[i] - ymin) * idy - 0.5
        gz = (z[i] - zmin) * idz - 0.5

        ix0 = int(np.floor(gx))
        iy0 = int(np.floor(gy))
        iz0 = int(np.floor(gz))

        tx = gx - ix0
        ty = gy - iy0
        tz = gz - iz0

        ix1 = ix0 + 1
        iy1 = iy0 + 1
        iz1 = iz0 + 1

        wx0 = 1.0 - tx
        wx1 = tx

        wy0 = 1.0 - ty
        wy1 = ty

        wz0 = 1.0 - tz
        wz1 = tz

        if periodx:
            if ix0 < 0:
                ix0 += nxgrid
            elif ix0 >= nxgrid:
                ix0 -= nxgrid

            if ix1 < 0:
                ix1 += nxgrid
            elif ix1 >= nxgrid:
                ix1 -= nxgrid

        if periody:
            if iy0 < 0:
                iy0 += nygrid
            elif iy0 >= nygrid:
                iy0 -= nygrid

            if iy1 < 0:
                iy1 += nygrid
            elif iy1 >= nygrid:
                iy1 -= nygrid

        if periodz:
            if iz0 < 0:
                iz0 += nzgrid
            elif iz0 >= nzgrid:
                iz0 -= nzgrid

            if iz1 < 0:
                iz1 += nzgrid
            elif iz1 >= nzgrid:
                iz1 -= nzgrid

        fp = norm

        if (
            (0 <= ix0 < nxgrid)
            and (0 <= iy0 < nygrid)
            and (0 <= iz0 < nzgrid)
        ):
            fgrid[ix0, iy0, iz0] += fp * wx0 * wy0 * wz0

        if (
            (0 <= ix0 < nxgrid)
            and (0 <= iy0 < nygrid)
            and (0 <= iz1 < nzgrid)
        ):
            fgrid[ix0, iy0, iz1] += fp * wx0 * wy0 * wz1

        if (
            (0 <= ix0 < nxgrid)
            and (0 <= iy1 < nygrid)
            and (0 <= iz0 < nzgrid)
        ):
            fgrid[ix0, iy1, iz0] += fp * wx0 * wy1 * wz0

        if (
            (0 <= ix0 < nxgrid)
            and (0 <= iy1 < nygrid)
            and (0 <= iz1 < nzgrid)
        ):
            fgrid[ix0, iy1, iz1] += fp * wx0 * wy1 * wz1

        if (
            (0 <= ix1 < nxgrid)
            and (0 <= iy0 < nygrid)
            and (0 <= iz0 < nzgrid)
        ):
            fgrid[ix1, iy0, iz0] += fp * wx1 * wy0 * wz0

        if (
            (0 <= ix1 < nxgrid)
            and (0 <= iy0 < nygrid)
            and (0 <= iz1 < nzgrid)
        ):
            fgrid[ix1, iy0, iz1] += fp * wx1 * wy0 * wz1

        if (
            (0 <= ix1 < nxgrid)
            and (0 <= iy1 < nygrid)
            and (0 <= iz0 < nzgrid)
        ):
            fgrid[ix1, iy1, iz0] += fp * wx1 * wy1 * wz0

        if (
            (0 <= ix1 < nxgrid)
            and (0 <= iy1 < nygrid)
            and (0 <= iz1 < nzgrid)
        ):
            fgrid[ix1, iy1, iz1] += fp * wx1 * wy1 * wz1

    return fgrid


@njit
def part2grid_tsc_3d(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    f: np.ndarray,
    xlength: float,
    ylength: float,
    zlength: float,
    xmin: float,
    ymin: float,
    zmin: float,
    nxgrid: int,
    nygrid: int,
    nzgrid: int,
    periodx: bool,
    periody: bool,
    periodz: bool,
    dtype: np.dtype = np.float64
) -> np.ndarray:
    """
    Assign particle values to a three-dimensional grid using the
    triangular-shaped-cloud (TSC) assignment scheme.

    Each particle contributes to three neighbouring cell centres along
    each coordinate axis, giving a total stencil of 27 grid cells in
    three dimensions. The one-dimensional TSC weights are calculated
    from the particle displacement relative to the centre of its
    containing grid cell.

    The deposited field is normalized by the grid-cell volume. When the
    full stencil lies within the domain, or periodic boundaries are
    applied, the weights associated with each particle sum to

        1 / (dx * dy * dz).

    Parameters
    ----------
    x, y, z : ndarray
        One-dimensional arrays containing the Cartesian coordinates of
        the particles. All arrays must have the same length.
    f : ndarray
        One-dimensional array containing the value associated with each
        particle. Must have the same length as ``x``, ``y`` and ``z``.
    xlength, ylength, zlength : float
        Physical lengths of the grid domain along the x, y and z axes.
    xmin, ymin, zmin : float
        Minimum coordinates of the grid domain along the x, y and z axes.
    nxgrid, nygrid, nzgrid : int
        Number of grid cells along the x, y and z axes.
    periodx, periody, periodz : bool
        Whether periodic boundary conditions are applied along the
        corresponding axis.
    dtype : numpy dtype, optional
        Data type used for the output grid. Default is ``np.float64``.

    Returns
    -------
    fgrid : ndarray
        Three-dimensional array of shape
        ``(nxgrid, nygrid, nzgrid)`` containing the TSC-assigned field.

    Notes
    -----
    For each axis, the particle position is expressed relative to the
    centre of its containing cell. If

        s = (x - xmin) / dx - floor((x - xmin) / dx) - 0.5,

    then ``s`` lies in the interval [-0.5, 0.5), and the three
    dimensionless TSC weights are

        w_- = 0.5 * (0.5 - s)**2
        w_0 = 0.75 - s**2
        w_+ = 0.5 * (0.5 + s)**2.

    For non-periodic axes, stencil contributions falling outside the
    grid domain are discarded.
    """
    npart = len(x)

    idx = nxgrid / xlength
    idy = nygrid / ylength
    idz = nzgrid / zlength

    norm = idx * idy * idz

    fgrid = np.zeros((nxgrid, nygrid, nzgrid), dtype=dtype)

    for i in range(npart):

        gx = (x[i] - xmin) * idx
        gy = (y[i] - ymin) * idy
        gz = (z[i] - zmin) * idz

        ixc = int(np.floor(gx))
        iyc = int(np.floor(gy))
        izc = int(np.floor(gz))

        sx = gx - ixc - 0.5
        sy = gy - iyc - 0.5
        sz = gz - izc - 0.5

        ixm = ixc - 1
        ixp = ixc + 1

        iym = iyc - 1
        iyp = iyc + 1

        izm = izc - 1
        izp = izc + 1

        wxm = 0.5 * (0.5 - sx) * (0.5 - sx)
        wxc = 0.75 - sx * sx
        wxp = 0.5 * (0.5 + sx) * (0.5 + sx)

        wym = 0.5 * (0.5 - sy) * (0.5 - sy)
        wyc = 0.75 - sy * sy
        wyp = 0.5 * (0.5 + sy) * (0.5 + sy)

        wzm = 0.5 * (0.5 - sz) * (0.5 - sz)
        wzc = 0.75 - sz * sz
        wzp = 0.5 * (0.5 + sz) * (0.5 + sz)

        if periodx:
            if ixm < 0:
                ixm += nxgrid
            elif ixm >= nxgrid:
                ixm -= nxgrid

            if ixc < 0:
                ixc += nxgrid
            elif ixc >= nxgrid:
                ixc -= nxgrid

            if ixp < 0:
                ixp += nxgrid
            elif ixp >= nxgrid:
                ixp -= nxgrid

        if periody:
            if iym < 0:
                iym += nygrid
            elif iym >= nygrid:
                iym -= nygrid

            if iyc < 0:
                iyc += nygrid
            elif iyc >= nygrid:
                iyc -= nygrid

            if iyp < 0:
                iyp += nygrid
            elif iyp >= nygrid:
                iyp -= nygrid

        if periodz:
            if izm < 0:
                izm += nzgrid
            elif izm >= nzgrid:
                izm -= nzgrid

            if izc < 0:
                izc += nzgrid
            elif izc >= nzgrid:
                izc -= nzgrid

            if izp < 0:
                izp += nzgrid
            elif izp >= nzgrid:
                izp -= nzgrid

        fp = f[i] * norm

        ix = (ixm, ixc, ixp)
        iy = (iym, iyc, iyp)
        iz = (izm, izc, izp)

        wx = (wxm, wxc, wxp)
        wy = (wym, wyc, wyp)
        wz = (wzm, wzc, wzp)

        for jx in range(3):
            if (ix[jx] < 0) or (ix[jx] >= nxgrid):
                continue

            for jy in range(3):
                if (iy[jy] < 0) or (iy[jy] >= nygrid):
                    continue

                wxy = fp * wx[jx] * wy[jy]

                for jz in range(3):
                    if (iz[jz] < 0) or (iz[jz] >= nzgrid):
                        continue

                    fgrid[ix[jx], iy[jy], iz[jz]] += wxy * wz[jz]

    return fgrid


@njit
def part2grid_tsc_3d_unit(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    xlength: float,
    ylength: float,
    zlength: float,
    xmin: float,
    ymin: float,
    zmin: float,
    nxgrid: int,
    nygrid: int,
    nzgrid: int,
    periodx: bool,
    periody: bool,
    periodz: bool,
    dtype: np.dtype = np.float64
) -> np.ndarray:
    """
    Assign particle values to a three-dimensional grid using the
    triangular-shaped-cloud (TSC) assignment scheme.

    Each particle contributes to three neighbouring cell centres along
    each coordinate axis, giving a total stencil of 27 grid cells in
    three dimensions. The one-dimensional TSC weights are calculated
    from the particle displacement relative to the centre of its
    containing grid cell.

    The deposited field is normalized by the grid-cell volume. When the
    full stencil lies within the domain, or periodic boundaries are
    applied, the weights associated with each particle sum to

        1 / (dx * dy * dz).

    Parameters
    ----------
    x, y, z : ndarray
        One-dimensional arrays containing the Cartesian coordinates of
        the particles. All arrays must have the same length.
    xlength, ylength, zlength : float
        Physical lengths of the grid domain along the x, y and z axes.
    xmin, ymin, zmin : float
        Minimum coordinates of the grid domain along the x, y and z axes.
    nxgrid, nygrid, nzgrid : int
        Number of grid cells along the x, y and z axes.
    periodx, periody, periodz : bool
        Whether periodic boundary conditions are applied along the
        corresponding axis.
    dtype : numpy dtype, optional
        Data type used for the output grid. Default is ``np.float64``.

    Returns
    -------
    fgrid : ndarray
        Three-dimensional array of shape
        ``(nxgrid, nygrid, nzgrid)`` containing the TSC-assigned field.

    Notes
    -----
    For each axis, the particle position is expressed relative to the
    centre of its containing cell. If

        s = (x - xmin) / dx - floor((x - xmin) / dx) - 0.5,

    then ``s`` lies in the interval [-0.5, 0.5), and the three
    dimensionless TSC weights are

        w_- = 0.5 * (0.5 - s)**2
        w_0 = 0.75 - s**2
        w_+ = 0.5 * (0.5 + s)**2.

    For non-periodic axes, stencil contributions falling outside the
    grid domain are discarded.
    """
    npart = len(x)

    idx = nxgrid / xlength
    idy = nygrid / ylength
    idz = nzgrid / zlength

    norm = idx * idy * idz

    fgrid = np.zeros((nxgrid, nygrid, nzgrid), dtype=dtype)

    for i in range(npart):

        gx = (x[i] - xmin) * idx
        gy = (y[i] - ymin) * idy
        gz = (z[i] - zmin) * idz

        ixc = int(np.floor(gx))
        iyc = int(np.floor(gy))
        izc = int(np.floor(gz))

        sx = gx - ixc - 0.5
        sy = gy - iyc - 0.5
        sz = gz - izc - 0.5

        ixm = ixc - 1
        ixp = ixc + 1

        iym = iyc - 1
        iyp = iyc + 1

        izm = izc - 1
        izp = izc + 1

        wxm = 0.5 * (0.5 - sx) * (0.5 - sx)
        wxc = 0.75 - sx * sx
        wxp = 0.5 * (0.5 + sx) * (0.5 + sx)

        wym = 0.5 * (0.5 - sy) * (0.5 - sy)
        wyc = 0.75 - sy * sy
        wyp = 0.5 * (0.5 + sy) * (0.5 + sy)

        wzm = 0.5 * (0.5 - sz) * (0.5 - sz)
        wzc = 0.75 - sz * sz
        wzp = 0.5 * (0.5 + sz) * (0.5 + sz)

        if periodx:
            if ixm < 0:
                ixm += nxgrid
            elif ixm >= nxgrid:
                ixm -= nxgrid

            if ixc < 0:
                ixc += nxgrid
            elif ixc >= nxgrid:
                ixc -= nxgrid

            if ixp < 0:
                ixp += nxgrid
            elif ixp >= nxgrid:
                ixp -= nxgrid

        if periody:
            if iym < 0:
                iym += nygrid
            elif iym >= nygrid:
                iym -= nygrid

            if iyc < 0:
                iyc += nygrid
            elif iyc >= nygrid:
                iyc -= nygrid

            if iyp < 0:
                iyp += nygrid
            elif iyp >= nygrid:
                iyp -= nygrid

        if periodz:
            if izm < 0:
                izm += nzgrid
            elif izm >= nzgrid:
                izm -= nzgrid

            if izc < 0:
                izc += nzgrid
            elif izc >= nzgrid:
                izc -= nzgrid

            if izp < 0:
                izp += nzgrid
            elif izp >= nzgrid:
                izp -= nzgrid

        fp = norm

        ix = (ixm, ixc, ixp)
        iy = (iym, iyc, iyp)
        iz = (izm, izc, izp)

        wx = (wxm, wxc, wxp)
        wy = (wym, wyc, wyp)
        wz = (wzm, wzc, wzp)

        for jx in range(3):
            if (ix[jx] < 0) or (ix[jx] >= nxgrid):
                continue

            for jy in range(3):
                if (iy[jy] < 0) or (iy[jy] >= nygrid):
                    continue

                wxy = fp * wx[jx] * wy[jy]

                for jz in range(3):
                    if (iz[jz] < 0) or (iz[jz] >= nzgrid):
                        continue

                    fgrid[ix[jx], iy[jy], iz[jz]] += wxy * wz[jz]

    return fgrid


@njit
def part2grid_pcs_3d(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    f: np.ndarray,
    xlength: float,
    ylength: float,
    zlength: float,
    xmin: float,
    ymin: float,
    zmin: float,
    nxgrid: int,
    nygrid: int,
    nzgrid: int,
    periodx: bool,
    periody: bool,
    periodz: bool,
    dtype: np.dtype = np.float64
) -> np.ndarray:
    """
    Assign particle values to a three-dimensional grid using the
    piecewise-cubic-spline (PCS) assignment scheme.

    Each particle contributes to four neighbouring cell centres along
    each coordinate axis, giving a total stencil of 64 grid cells in
    three dimensions. The one-dimensional PCS weights correspond to a
    cubic B-spline assignment kernel.

    The deposited field is normalized by the grid-cell volume. When the
    complete stencil lies within the domain, or periodic boundaries are
    applied, the weights associated with each particle sum to

        1 / (dx * dy * dz).

    Parameters
    ----------
    x, y, z : ndarray
        One-dimensional arrays containing the Cartesian coordinates of
        the particles. All arrays must have the same length.
    f : ndarray
        One-dimensional array containing the value associated with each
        particle. Must have the same length as ``x``, ``y`` and ``z``.
    xlength, ylength, zlength : float
        Physical lengths of the grid domain along the x, y and z axes.
    xmin, ymin, zmin : float
        Minimum coordinates of the grid domain along the x, y and z axes.
    nxgrid, nygrid, nzgrid : int
        Number of grid cells along the x, y and z axes.
    periodx, periody, periodz : bool
        Whether periodic boundary conditions are applied along the
        corresponding axis.
    dtype : np.dtype, optional
        Data type of the returned interpolation array.
    
    Returns
    -------
    fgrid : ndarray
        Three-dimensional array of shape
        ``(nxgrid, nygrid, nzgrid)`` containing the PCS-assigned field.

    Notes
    -----
    The grid is treated as cell-centred. For each axis, define

        g = (x - xmin) / dx - 0.5
        i = floor(g)
        t = g - i,

    where ``t`` lies in [0, 1). The four neighbouring grid cells are
    ``i - 1``, ``i``, ``i + 1`` and ``i + 2``, with dimensionless
    cubic B-spline weights

        w0 = (1 - t)**3 / 6
        w1 = (4 - 6*t**2 + 3*t**3) / 6
        w2 = (1 + 3*t + 3*t**2 - 3*t**3) / 6
        w3 = t**3 / 6.

    For non-periodic axes, stencil contributions falling outside the
    grid domain are discarded.
    """
    npart = len(x)

    idx = nxgrid / xlength
    idy = nygrid / ylength
    idz = nzgrid / zlength

    norm = idx * idy * idz

    fgrid = np.zeros((nxgrid, nygrid, nzgrid), dtype=dtype)

    for i in range(npart):

        gx = (x[i] - xmin) * idx - 0.5
        gy = (y[i] - ymin) * idy - 0.5
        gz = (z[i] - zmin) * idz - 0.5

        ix1 = int(np.floor(gx))
        iy1 = int(np.floor(gy))
        iz1 = int(np.floor(gz))

        tx = gx - ix1
        ty = gy - iy1
        tz = gz - iz1

        ix0 = ix1 - 1
        ix2 = ix1 + 1
        ix3 = ix1 + 2

        iy0 = iy1 - 1
        iy2 = iy1 + 1
        iy3 = iy1 + 2

        iz0 = iz1 - 1
        iz2 = iz1 + 1
        iz3 = iz1 + 2

        tx2 = tx * tx
        ty2 = ty * ty
        tz2 = tz * tz

        tx3 = tx2 * tx
        ty3 = ty2 * ty
        tz3 = tz2 * tz

        omtx = 1.0 - tx
        omty = 1.0 - ty
        omtz = 1.0 - tz

        wx0 = omtx * omtx * omtx / 6.0
        wx1 = (4.0 - 6.0 * tx2 + 3.0 * tx3) / 6.0
        wx2 = (1.0 + 3.0 * tx + 3.0 * tx2 - 3.0 * tx3) / 6.0
        wx3 = tx3 / 6.0

        wy0 = omty * omty * omty / 6.0
        wy1 = (4.0 - 6.0 * ty2 + 3.0 * ty3) / 6.0
        wy2 = (1.0 + 3.0 * ty + 3.0 * ty2 - 3.0 * ty3) / 6.0
        wy3 = ty3 / 6.0

        wz0 = omtz * omtz * omtz / 6.0
        wz1 = (4.0 - 6.0 * tz2 + 3.0 * tz3) / 6.0
        wz2 = (1.0 + 3.0 * tz + 3.0 * tz2 - 3.0 * tz3) / 6.0
        wz3 = tz3 / 6.0

        if periodx:
            if ix0 < 0:
                ix0 += nxgrid
            elif ix0 >= nxgrid:
                ix0 -= nxgrid

            if ix1 < 0:
                ix1 += nxgrid
            elif ix1 >= nxgrid:
                ix1 -= nxgrid

            if ix2 < 0:
                ix2 += nxgrid
            elif ix2 >= nxgrid:
                ix2 -= nxgrid

            if ix3 < 0:
                ix3 += nxgrid
            elif ix3 >= nxgrid:
                ix3 -= nxgrid

        if periody:
            if iy0 < 0:
                iy0 += nygrid
            elif iy0 >= nygrid:
                iy0 -= nygrid

            if iy1 < 0:
                iy1 += nygrid
            elif iy1 >= nygrid:
                iy1 -= nygrid

            if iy2 < 0:
                iy2 += nygrid
            elif iy2 >= nygrid:
                iy2 -= nygrid

            if iy3 < 0:
                iy3 += nygrid
            elif iy3 >= nygrid:
                iy3 -= nygrid

        if periodz:
            if iz0 < 0:
                iz0 += nzgrid
            elif iz0 >= nzgrid:
                iz0 -= nzgrid

            if iz1 < 0:
                iz1 += nzgrid
            elif iz1 >= nzgrid:
                iz1 -= nzgrid

            if iz2 < 0:
                iz2 += nzgrid
            elif iz2 >= nzgrid:
                iz2 -= nzgrid

            if iz3 < 0:
                iz3 += nzgrid
            elif iz3 >= nzgrid:
                iz3 -= nzgrid

        fp = f[i] * norm

        ix = (ix0, ix1, ix2, ix3)
        iy = (iy0, iy1, iy2, iy3)
        iz = (iz0, iz1, iz2, iz3)

        wx = (wx0, wx1, wx2, wx3)
        wy = (wy0, wy1, wy2, wy3)
        wz = (wz0, wz1, wz2, wz3)

        for jx in range(4):
            if (ix[jx] < 0) or (ix[jx] >= nxgrid):
                continue

            for jy in range(4):
                if (iy[jy] < 0) or (iy[jy] >= nygrid):
                    continue

                wxy = fp * wx[jx] * wy[jy]

                for jz in range(4):
                    if (iz[jz] < 0) or (iz[jz] >= nzgrid):
                        continue

                    fgrid[ix[jx], iy[jy], iz[jz]] += wxy * wz[jz]

    return fgrid


@njit
def part2grid_pcs_3d_unit(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    xlength: float,
    ylength: float,
    zlength: float,
    xmin: float,
    ymin: float,
    zmin: float,
    nxgrid: int,
    nygrid: int,
    nzgrid: int,
    periodx: bool,
    periody: bool,
    periodz: bool,
    dtype: np.dtype = np.float64
) -> np.ndarray:
    """
    Assign particle values to a three-dimensional grid using the
    piecewise-cubic-spline (PCS) assignment scheme.

    Each particle contributes to four neighbouring cell centres along
    each coordinate axis, giving a total stencil of 64 grid cells in
    three dimensions. The one-dimensional PCS weights correspond to a
    cubic B-spline assignment kernel.

    The deposited field is normalized by the grid-cell volume. When the
    complete stencil lies within the domain, or periodic boundaries are
    applied, the weights associated with each particle sum to

        1 / (dx * dy * dz).

    Parameters
    ----------
    x, y, z : ndarray
        One-dimensional arrays containing the Cartesian coordinates of
        the particles. All arrays must have the same length.
    xlength, ylength, zlength : float
        Physical lengths of the grid domain along the x, y and z axes.
    xmin, ymin, zmin : float
        Minimum coordinates of the grid domain along the x, y and z axes.
    nxgrid, nygrid, nzgrid : int
        Number of grid cells along the x, y and z axes.
    periodx, periody, periodz : bool
        Whether periodic boundary conditions are applied along the
        corresponding axis.
    dtype : np.dtype, optional
        Data type of the returned interpolation array.

    Returns
    -------
    fgrid : ndarray
        Three-dimensional array of shape
        ``(nxgrid, nygrid, nzgrid)`` containing the PCS-assigned field.

    Notes
    -----
    The grid is treated as cell-centred. For each axis, define

        g = (x - xmin) / dx - 0.5
        i = floor(g)
        t = g - i,

    where ``t`` lies in [0, 1). The four neighbouring grid cells are
    ``i - 1``, ``i``, ``i + 1`` and ``i + 2``, with dimensionless
    cubic B-spline weights

        w0 = (1 - t)**3 / 6
        w1 = (4 - 6*t**2 + 3*t**3) / 6
        w2 = (1 + 3*t + 3*t**2 - 3*t**3) / 6
        w3 = t**3 / 6.

    For non-periodic axes, stencil contributions falling outside the
    grid domain are discarded.
    """
    npart = len(x)

    idx = nxgrid / xlength
    idy = nygrid / ylength
    idz = nzgrid / zlength

    norm = idx * idy * idz

    fgrid = np.zeros((nxgrid, nygrid, nzgrid), dtype=dtype)

    for i in range(npart):

        gx = (x[i] - xmin) * idx - 0.5
        gy = (y[i] - ymin) * idy - 0.5
        gz = (z[i] - zmin) * idz - 0.5

        ix1 = int(np.floor(gx))
        iy1 = int(np.floor(gy))
        iz1 = int(np.floor(gz))

        tx = gx - ix1
        ty = gy - iy1
        tz = gz - iz1

        ix0 = ix1 - 1
        ix2 = ix1 + 1
        ix3 = ix1 + 2

        iy0 = iy1 - 1
        iy2 = iy1 + 1
        iy3 = iy1 + 2

        iz0 = iz1 - 1
        iz2 = iz1 + 1
        iz3 = iz1 + 2

        tx2 = tx * tx
        ty2 = ty * ty
        tz2 = tz * tz

        tx3 = tx2 * tx
        ty3 = ty2 * ty
        tz3 = tz2 * tz

        omtx = 1.0 - tx
        omty = 1.0 - ty
        omtz = 1.0 - tz

        wx0 = omtx * omtx * omtx / 6.0
        wx1 = (4.0 - 6.0 * tx2 + 3.0 * tx3) / 6.0
        wx2 = (1.0 + 3.0 * tx + 3.0 * tx2 - 3.0 * tx3) / 6.0
        wx3 = tx3 / 6.0

        wy0 = omty * omty * omty / 6.0
        wy1 = (4.0 - 6.0 * ty2 + 3.0 * ty3) / 6.0
        wy2 = (1.0 + 3.0 * ty + 3.0 * ty2 - 3.0 * ty3) / 6.0
        wy3 = ty3 / 6.0

        wz0 = omtz * omtz * omtz / 6.0
        wz1 = (4.0 - 6.0 * tz2 + 3.0 * tz3) / 6.0
        wz2 = (1.0 + 3.0 * tz + 3.0 * tz2 - 3.0 * tz3) / 6.0
        wz3 = tz3 / 6.0

        if periodx:
            if ix0 < 0:
                ix0 += nxgrid
            elif ix0 >= nxgrid:
                ix0 -= nxgrid

            if ix1 < 0:
                ix1 += nxgrid
            elif ix1 >= nxgrid:
                ix1 -= nxgrid

            if ix2 < 0:
                ix2 += nxgrid
            elif ix2 >= nxgrid:
                ix2 -= nxgrid

            if ix3 < 0:
                ix3 += nxgrid
            elif ix3 >= nxgrid:
                ix3 -= nxgrid

        if periody:
            if iy0 < 0:
                iy0 += nygrid
            elif iy0 >= nygrid:
                iy0 -= nygrid

            if iy1 < 0:
                iy1 += nygrid
            elif iy1 >= nygrid:
                iy1 -= nygrid

            if iy2 < 0:
                iy2 += nygrid
            elif iy2 >= nygrid:
                iy2 -= nygrid

            if iy3 < 0:
                iy3 += nygrid
            elif iy3 >= nygrid:
                iy3 -= nygrid

        if periodz:
            if iz0 < 0:
                iz0 += nzgrid
            elif iz0 >= nzgrid:
                iz0 -= nzgrid

            if iz1 < 0:
                iz1 += nzgrid
            elif iz1 >= nzgrid:
                iz1 -= nzgrid

            if iz2 < 0:
                iz2 += nzgrid
            elif iz2 >= nzgrid:
                iz2 -= nzgrid

            if iz3 < 0:
                iz3 += nzgrid
            elif iz3 >= nzgrid:
                iz3 -= nzgrid

        fp = norm

        ix = (ix0, ix1, ix2, ix3)
        iy = (iy0, iy1, iy2, iy3)
        iz = (iz0, iz1, iz2, iz3)

        wx = (wx0, wx1, wx2, wx3)
        wy = (wy0, wy1, wy2, wy3)
        wz = (wz0, wz1, wz2, wz3)

        for jx in range(4):
            if (ix[jx] < 0) or (ix[jx] >= nxgrid):
                continue

            for jy in range(4):
                if (iy[jy] < 0) or (iy[jy] >= nygrid):
                    continue

                wxy = fp * wx[jx] * wy[jy]

                for jz in range(4):
                    if (iz[jz] < 0) or (iz[jz] >= nzgrid):
                        continue

                    fgrid[ix[jx], iy[jy], iz[jz]] += wxy * wz[jz]

    return fgrid


def part2grid3D(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    boxsize: Union[float, List[float]],
    ngrid: Union[float, List[float]],
    f: np.ndarray = None,
    method: str = "TSC",
    periodic: Union[bool, List[bool]] = True,
    origin: Union[float, List[float]] = 0.0,
    dtype: np.dtype = np.float64
) -> np.ndarray:
    """
    Returns the density contrast for the nearest grid point grid assignment.

    Parameters
    ----------
    x : array
        X coordinates of the particle.
    y : array
        Y coordinates of the particle.
    z : array
        Z coordinates of the particle.
    boxsize : float or list
        Box size.
    ngrid : int or list
        Grid divisions across one axis.
    f : array, optional
        Value of each particle to be assigned to the grid.
    method : str, optional
        Grid assignment scheme, either 'NGP', 'CIC', 'TSC' or 'PCS'.
    periodic : bool or list, optional
        Assign particles with periodic boundaries.
    origin : float or list, optional
        Origin.
    dtype : numpy dtype, optional
        Data type used for the output grid. Default is ``np.float64``.

    Returns
    -------
    fgrid : array
        Grid assigned values.
    """
    if np.isscalar(boxsize):
        xlength, ylength, zlength = boxsize, boxsize, boxsize
    else:
        xlength, ylength, zlength = boxsize[0], boxsize[1], boxsize[2]
    if np.isscalar(origin):
        xmin = origin
        ymin = origin
        zmin = origin
    else:
        xmin, ymin, zmin = origin[0], origin[1], origin[2]
    if np.isscalar(ngrid):
        nxgrid, nygrid, nzgrid = int(ngrid), int(ngrid), int(ngrid)
    else:
        nxgrid, nygrid, nzgrid = int(ngrid[0]), int(ngrid[1]), int(ngrid[2])
    if np.isscalar(periodic):
        periodx = periodic
        periody = periodic
        periodz = periodic
    else:
        periodx, periody, periodz = periodic[0], periodic[1], periodic[2]
    if method == "NGP":
        if f is None:
            fgrid = part2grid_ngp_3d_unit(
                x, y, z, xlength, ylength, zlength, xmin, ymin, zmin, 
                nxgrid, nygrid, nzgrid, dtype=dtype
            )
        else:
            fgrid = part2grid_ngp_3d(
                x, y, z, f, xlength, ylength, zlength, xmin, ymin, zmin, 
                nxgrid, nygrid, nzgrid, dtype=dtype
            )
    elif method == "CIC":
        if f is None:
            fgrid = part2grid_cic_3d_unit(
                x, y, z, xlength, ylength, zlength, xmin, ymin, zmin,
                nxgrid, nygrid, nzgrid, periodx, periody, periodz, dtype=dtype
            )
        else:
            fgrid = part2grid_cic_3d(
                x, y, z, f, xlength, ylength, zlength, xmin, ymin, zmin,
                nxgrid, nygrid, nzgrid, periodx, periody, periodz, dtype=dtype
            )
    elif method == "TSC":
        if f is None:
            fgrid = part2grid_tsc_3d_unit(
                x, y, z, xlength, ylength, zlength, xmin, ymin, zmin,
                nxgrid, nygrid, nzgrid, periodx, periody, periodz, dtype=dtype
            )
        else:
            fgrid = part2grid_tsc_3d(
                x, y, z, f, xlength, ylength, zlength, xmin, ymin, zmin,
                nxgrid, nygrid, nzgrid, periodx, periody, periodz, dtype=dtype
            )
    elif method == "PCS":
        if f is None:
            fgrid = part2grid_pcs_3d_unit(
                x, y, z, xlength, ylength, zlength, xmin, ymin, zmin,
                nxgrid, nygrid, nzgrid, periodx, periody, periodz, dtype=dtype
            )
        else:
            fgrid = part2grid_pcs_3d(
                x, y, z, f, xlength, ylength, zlength, xmin, ymin, zmin,
                nxgrid, nygrid, nzgrid, periodx, periody, periodz, dtype=dtype
            )
    return fgrid
