import numpy as np
from numba import njit


@njit
def dfdx_1d_periodic(f, boxsize):
    """
    Numerical differentiation assuming a uniform grid with periodic boundaries in 1D.

    Parameters
    ----------
    f : array
        Field values on the grid.
    boxsize : float
        Size of the grid.
    ngrid : int
        Number of points on the grid.

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
        df[i] = (f[i2] - f[i1])/(2*dl)
    return df


@njit
def dfdx_2d_periodic(f, boxsize, ngrid):
    """
    Numerical differentiation in x assuming a uniform grid with periodic boundaries in 2D.

    Parameters
    ----------
    f : array
        Field values on the grid.
    boxsize : float
        Size of the grid.
    ngrid : int
        Number of points on the grid.

    Returns
    -------
    df : array
        Differentiated function along x.
    """
    dl = boxsize / float(ngrid)
    df = np.zeros(ngrid*ngrid)
    for i in range(0, ngrid):
        if i == 0:
            i1 = ngrid - 1
        else:
            i1 = i - 1
        if i == ngrid - 1:
            i2 = 0
        else:
            i2 = i + 1
        for j in range(0, ngrid):
            df[j + ngrid*i] = (f[j + ngrid*i2] - f[j + ngrid*i1])/(2*dl)
    return df


@njit
def dfdy_2d_periodic(f, boxsize, ngrid):
    """
    Numerical differentiation in y assuming a uniform grid with periodic boundaries in 2D.

    Parameters
    ----------
    f : array
        Field values on the grid.
    boxsize : float
        Size of the grid.
    ngrid : int
        Number of points on the grid.
    
    Returns
    -------
    df : array
        Differentiated function along x.
    """
    dl = boxsize / float(ngrid)
    df = np.zeros(ngrid*ngrid)
    for j in range(0, ngrid):
        if j == 0:
            j1 = ngrid - 1
        else:
            j1 = j - 1
        if j == ngrid - 1:
            j2 = 0
        else:
            j2 = j + 1
        for i in range(0, ngrid):
            df[j + ngrid*i] = (f[j2 + ngrid*i] - f[j1 + ngrid*i])/(2*dl)
    return df


@njit
def dfdx_3d_periodic(f, boxsize, ngrid):
    """
    Numerical differentiation in x assuming a uniform grid with periodic boundaries in 3D.

    Parameters
    ----------
    f : array
        Field values on the grid.
    boxsize : float
        Size of the grid.
    ngrid : int
        Number of points on the grid.
    
    Returns
    -------
    df : array
        Differentiated function along x.
    """
    dl = boxsize / float(ngrid)
    df = np.zeros(ngrid*ngrid*ngrid)
    for i in range(0, ngrid):
        if i == 0:
            i1 = ngrid - 1
        else:
            i1 = i - 1
        if i == ngrid - 1:
            i2 = 0
        else:
            i2 = i + 1
        for j in range(0, ngrid):
            for k in range(0, ngrid):
                df[k + ngrid*(j + ngrid*i)] = (f[k + ngrid*(j + ngrid*i2)] - f[k + ngrid*(j + ngrid*i1)])/(2*dl)
    return df


@njit
def dfdy_3d_periodic(f, boxsize, ngrid):
    """
    Numerical differentiation in x assuming a uniform grid with periodic boundaries in 3D.

    Parameters
    ----------
    f : array
        Field values on the grid.
    boxsize : float
        Size of the grid.
    ngrid : int
        Number of points on the grid.
    
    Returns
    -------
    df : array
        Differentiated function along x.
    """
    dl = boxsize / float(ngrid)
    df = np.zeros(ngrid*ngrid*ngrid)
    for j in range(0, ngrid):
        if j == 0:
            j1 = ngrid - 1
        else:
            j1 = j - 1
        if j == ngrid - 1:
            j2 = 0
        else:
            j2 = j + 1
        for i in range(0, ngrid):
            for k in range(0, ngrid):
                df[k + ngrid*(j + ngrid*i)] = (f[k + ngrid*(j2 + ngrid*i)] - f[k + ngrid*(j1 + ngrid*i)])/(2*dl)
    return df


@njit
def dfdz_3d_periodic(f, boxsize, ngrid):
    """
    Numerical differentiation in x assuming a uniform grid with periodic boundaries in 3D.

    Parameters
    ----------
    f : array
        Field values on the grid.
    boxsize : float
        Size of the grid.
    ngrid : int
        Number of points on the grid.
    
    Returns
    -------
    df : array
        Differentiated function along x.
    """
    dl = boxsize / float(ngrid)
    df = np.zeros(ngrid*ngrid*ngrid)
    for k in range(0, ngrid):
        if k == 0:
            k1 = ngrid - 1
        else:
            k1 = k - 1
        if k == ngrid - 1:
            k2 = 0
        else:
            k2 = k + 1
        for j in range(0, ngrid):
            for i in range(0, ngrid):
                df[k + ngrid*(j + ngrid*i)] = (f[k2 + ngrid*(j + ngrid*i)] - f[k1 + ngrid*(j + ngrid*i)])/(2*dl)
    return df