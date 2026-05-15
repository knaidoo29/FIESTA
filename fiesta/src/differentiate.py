import numpy as np
from numba import njit


@njit
def dfdx_1d_periodic(f: np.ndarray, boxsize: float) -> np.ndarray:  # pragma: no cover
    """
    Numerical differentiation assuming a uniform grid with periodic boundaries in 1D.

    Parameters
    ----------
    f : array
        Field values on the grid.
    boxsize : float
        Size of the grid.

    Returns
    -------
    df : array
        Differentiated function along x.
    """
    ngrid = len(f)
    df = np.zeros(ngrid)
    dl = boxsize / float(ngrid)
    for i in range(0, ngrid):
        if i == 0:
            i1 = ngrid - 1
        else:
            i1 = i - 1
        if i == ngrid - 1:
            i2 = 0
        else:
            i2 = i + 1
        df[i] = (f[i2] - f[i1]) / (2 * dl)
    return df


@njit
def dfdx_2d_periodic(
    f: np.ndarray, xboxsize: float, nxgrid: int, nygrid: int
) -> np.ndarray:  # pragma: no cover
    """
    Numerical differentiation in x assuming a uniform grid with periodic boundaries in 2D.

    Parameters
    ----------
    f : array
        Field values on the grid.
    xboxsize : float
        Size of the x-axis grid
    nxgrid : int
        Number of grid points along the x-axis grid.
    nygrid : int
        Number of grid points along the y-axis grid.

    Returns
    -------
    df : array
        Differentiated function along x.
    """
    dx = xboxsize / float(nxgrid)
    df = np.zeros(nxgrid * nygrid)
    for i in range(0, nxgrid):
        if i == 0:
            i1 = nxgrid - 1
        else:
            i1 = i - 1
        if i == nxgrid - 1:
            i2 = 0
        else:
            i2 = i + 1
        for j in range(0, nygrid):
            df[j + nygrid * i] = (f[j + nygrid * i2] - f[j + nygrid * i1]) / (2 * dx)
    return df


@njit
def dfdy_2d_periodic(
    f: np.ndarray, yboxsize: float, nxgrid: int, nygrid: int
):  # pragma: no cover
    """
    Numerical differentiation in y assuming a uniform grid with periodic boundaries in 2D.

    Parameters
    ----------
    f : array
        Field values on the grid.
    yboxsize : float
        Size of the y-axis grid
    nxgrid : int
        Number of grid points along the x-axis grid.
    nygrid : int
        Number of grid points along the y-axis grid.

    Returns
    -------
    df : array
        Differentiated function along x.
    """
    dy = yboxsize / float(nygrid)
    df = np.zeros(nxgrid * nygrid)
    for j in range(0, nygrid):
        if j == 0:
            j1 = nygrid - 1
        else:
            j1 = j - 1
        if j == nygrid - 1:
            j2 = 0
        else:
            j2 = j + 1
        for i in range(0, nxgrid):
            df[j + nygrid * i] = (f[j2 + nygrid * i] - f[j1 + nygrid * i]) / (2 * dy)
    return df


@njit
def dfdx_3d_periodic(
    f: np.ndarray, xboxsize: float, nxgrid: int, nygrid: int, nzgrid: int
) -> np.ndarray:  # pragma: no cover
    """
    Numerical differentiation in x assuming a uniform grid with periodic boundaries in 3D.

    Parameters
    ----------
    f : array
        Field values on the grid.
    xboxsize : float
        Size of the x-axis grid
    nxgrid : int
        Number of grid points along the x-axis grid.
    nygrid : int
        Number of grid points along the y-axis grid.
    nzgrid : int
        Number of grid points along the z-axis grid.

    Returns
    -------
    df : array
        Differentiated function along x.
    """
    dx = xboxsize / float(nxgrid)
    df = np.zeros(nxgrid * nygrid * nzgrid)
    for i in range(0, nxgrid):
        if i == 0:
            i1 = nxgrid - 1
        else:
            i1 = i - 1
        if i == nxgrid - 1:
            i2 = 0
        else:
            i2 = i + 1
        for j in range(0, nygrid):
            for k in range(0, nzgrid):
                df[k + nzgrid * (j + nygrid * i)] = (
                    f[k + nzgrid * (j + nygrid * i2)]
                    - f[k + nzgrid * (j + nygrid * i1)]
                ) / (2 * dx)
    return df


@njit
def dfdy_3d_periodic(
    f: np.ndarray, yboxsize: float, nxgrid: int, nygrid: int, nzgrid: int
) -> np.ndarray:  # pragma: no cover
    """
    Numerical differentiation in x assuming a uniform grid with periodic boundaries in 3D.

    Parameters
    ----------
    f : array
        Field values on the grid.
    yboxsize : float
        Size of the y-axis grid
    nxgrid : int
        Number of grid points along the x-axis grid.
    nygrid : int
        Number of grid points along the y-axis grid.
    nzgrid : int
        Number of grid points along the z-axis grid.

    Returns
    -------
    df : array
        Differentiated function along x.
    """
    dy = yboxsize / float(nygrid)
    df = np.zeros(nxgrid * nygrid * nzgrid)
    for j in range(0, nygrid):
        if j == 0:
            j1 = nygrid - 1
        else:
            j1 = j - 1
        if j == nygrid - 1:
            j2 = 0
        else:
            j2 = j + 1
        for i in range(0, nxgrid):
            for k in range(0, nzgrid):
                df[k + nzgrid * (j + nygrid * i)] = (
                    f[k + nzgrid * (j2 + nygrid * i)]
                    - f[k + nzgrid * (j1 + nygrid * i)]
                ) / (2 * dy)
    return df


@njit
def dfdz_3d_periodic(
    f: np.ndarray, zboxsize: float, nxgrid: int, nygrid: int, nzgrid: int
) -> np.ndarray:  # pragma: no cover
    """
    Numerical differentiation in x assuming a uniform grid with periodic boundaries in 3D.

    Parameters
    ----------
    f : array
        Field values on the grid.
    zboxsize : float
        Size of the z-axis grid
    nxgrid : int
        Number of grid points along the x-axis grid.
    nygrid : int
        Number of grid points along the y-axis grid.
    nzgrid : int
        Number of grid points along the z-axis grid.

    Returns
    -------
    df : array
        Differentiated function along x.
    """
    dz = zboxsize / float(nzgrid)
    df = np.zeros(nxgrid * nygrid * nzgrid)
    for k in range(0, nzgrid):
        if k == 0:
            k1 = nzgrid - 1
        else:
            k1 = k - 1
        if k == nzgrid - 1:
            k2 = 0
        else:
            k2 = k + 1
        for j in range(0, nygrid):
            for i in range(0, nxgrid):
                df[k + nzgrid * (j + nygrid * i)] = (
                    f[k2 + nzgrid * (j + nygrid * i)]
                    - f[k1 + nzgrid * (j + nygrid * i)]
                ) / (2 * dz)
    return df
