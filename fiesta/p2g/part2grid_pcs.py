import numpy as np
from numba import njit


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
) -> np.ndarray: # pragma: no cover
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
) -> np.ndarray: # pragma: no cover
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
) -> np.ndarray: # pragma: no cover
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
) -> np.ndarray: # pragma: no cover
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

