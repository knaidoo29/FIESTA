import numpy as np
from numba import njit


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
) -> np.ndarray: # pragma: no cover
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
) -> np.ndarray: # pragma: no cover
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
) -> np.ndarray: # pragma: no cover
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
) -> np.ndarray: # pragma: no cover
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
