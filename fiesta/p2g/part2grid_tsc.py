import numpy as np
from numba import njit


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
) -> np.ndarray: # pragma: no cover
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
) -> np.ndarray: # pragma: no cover
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
) -> np.ndarray: # pragma: no cover
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
) -> np.ndarray: # pragma: no cover
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

