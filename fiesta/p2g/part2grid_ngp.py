import numpy as np
from numba import njit


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
) -> np.ndarray: # pragma: no cover
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
) -> np.ndarray: # pragma: no cover
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
) -> np.ndarray: # pragma: no cover
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
) -> np.ndarray: # pragma: no cover
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

