import numpy as np


def _dfdx(xgrid: np.ndarray, f: np.ndarray) -> np.ndarray:
    """
    Differentiate using a symmetric two-point derivative and the non-symmetric
    three point derivative estimator.

    Parameters
    ----------
    xgrid : array
        X-axis.
    f : array
        Function values at x.

    Returns
    -------
    df : array
        Numerical differentiation values for f evaluated at points x.

    Notes
    -----
    For non-boundary values:

    df   f(x + dx) - f(x - dx)
    -- = ---------------------
    dx            2dx

    For boundary values:

    df   - f(x + 2dx) + 4f(x + dx) - 3f(x)
    -- = ---------------------------------
    dx                  2dx
    """
    dx = xgrid[1] - xgrid[0]
    dfdx = np.zeros(np.shape(f))
    # boundary differentials
    dfdx[0] = (-f[2] + 4 * f[1] - 3.0 * f[0]) / (2.0 * dx)
    dfdx[-1] = (f[-3] - 4 * f[-2] + 3.0 * f[-1]) / (2.0 * dx)
    # non-boundary differentials
    dfdx[1:-1] = (f[2:] - f[:-2]) / (2.0 * dx)
    return dfdx


def _dfdy(ygrid: np.ndarray, f: np.ndarray) -> np.ndarray:
    """
    Differentiate using a symmetric two-point derivative and the non-symmetric
    three point derivative estimator.

    Parameters
    ----------
    ygrid : array
        Y-axis.
    f : array
        Function values at y.

    Returns
    -------
    df : array
        Numerical differentiation values for f evaluated at points y.
    """
    dy = ygrid[1] - ygrid[0]
    dfdy = np.zeros(np.shape(f))
    # boundary differentials
    dfdy[:, 0] = (-f[:, 2] + 4 * f[:, 1] - 3.0 * f[:, 0]) / (2.0 * dy)
    dfdy[:, -1] = (f[:, -3] - 4 * f[:, -2] + 3.0 * f[:, -1]) / (2.0 * dy)
    # non-boundary differentials
    dfdy[:, 1:-1] = (f[:, 2:] - f[:, :-2]) / (2.0 * dy)
    return dfdy


def _dfdz(zgrid: np.ndarray, f: np.ndarray) -> np.ndarray:
    """
    Differentiate using a symmetric two-point derivative and the non-symmetric
    three point derivative estimator.

    Parameters
    ----------
    zgrid : array
        Z-axis.
    f : array
        Function values at z.

    Returns
    -------
    df : array
        Numerical differentiation values for f evaluated at points z.
    """
    dz = zgrid[1] - zgrid[0]
    dfdz = np.zeros(np.shape(f))
    # boundary differentials
    dfdz[:, :, 0] = (-f[:, :, 2] + 4 * f[:, :, 1] - 3.0 * f[:, :, 0]) / (2.0 * dz)
    dfdz[:, :, -1] = (f[:, :, -3] - 4 * f[:, :, -2] + 3.0 * f[:, :, -1]) / (2.0 * dz)
    # non-boundary differentials
    dfdz[:, :, 1:-1] = (f[:, :, 2:] - f[:, :, :-2]) / (2.0 * dz)
    return dfdz


def dfdx(xgrid: np.ndarray, f: np.ndarray, periodic: bool = False) -> np.ndarray:
    """
    Differentiate using a symmetric two-point derivative and the non-symmetric
    three point derivative estimator.

    Parameters
    ----------
    xgrid : array
        X-axis.
    f : array
        Function values at x.
    periodic : bool, optional
        Sets periodic boundary conditions.

    Returns
    -------
    dfdx : array
        Numerical differentiation values for f evaluated at points x.
    """
    if periodic is True:
        dx = xgrid[1] - xgrid[0]
        _xgrid = np.zeros(len(xgrid) + 2)
        shape = np.array(np.shape(f))
        shape[0] += 2
        _f = np.zeros(shape)
        _xgrid[1:-1] = xgrid
        _xgrid[0] = xgrid[0] - dx
        _xgrid[-1] = xgrid[-1] + dx
        _f[1:-1] = f
        _f[0] = f[-1]
        _f[-1] = f[0]
    else:
        _xgrid, _f = xgrid, f
    dfdx = _dfdx(_xgrid, _f)
    if periodic is True:
        dfdx = dfdx[1:-1]
    return dfdx


def dfdy(ygrid: np.ndarray, f: np.ndarray, periodic: bool = False) -> np.ndarray:
    """
    Differentiate using a symmetric two-point derivative and the non-symmetric
    three point derivative estimator.

    Parameters
    ----------
    ygrid : array
        Y-axis.
    f : array
        Function values at y.
    periodic : bool, optional
        Sets periodic boundary conditions.

    Returns
    -------
    dfdy : array
        Numerical differentiation values for f evaluated at points y.
    """
    if periodic is True:
        dy = ygrid[1] - ygrid[0]
        _ygrid = np.zeros(len(ygrid) + 2)
        shape = np.array(np.shape(f))
        shape[1] += 2
        _f = np.zeros(shape)
        _ygrid[1:-1] = ygrid
        _ygrid[0] = ygrid[0] - dy
        _ygrid[-1] = ygrid[-1] + dy
        _f[:, 1:-1] = f
        _f[:, 0] = f[:, -1]
        _f[:, -1] = f[:, 0]
    else:
        _ygrid, _f = ygrid, f
    dfdy = _dfdy(_ygrid, _f)
    if periodic is True:
        dfdy = dfdy[:, 1:-1]
    return dfdy


def dfdz(zgrid: np.ndarray, f: np.ndarray, periodic: bool = False) -> np.ndarray:
    """
    Differentiate using a symmetric two-point derivative and the non-symmetric
    three point derivative estimator.

    Parameters
    ----------
    zgrid : array
        Z-axis.
    f : array
        Function values at z.
    periodic : bool, optional
        Sets periodic boundary conditions.

    Returns
    -------
    dfdz : array
        Numerical differentiation values for f evaluated at points z.
    """
    if periodic is True:
        dz = zgrid[1] - zgrid[0]
        _zgrid = np.zeros(len(zgrid) + 2)
        shape = np.array(np.shape(f))
        shape[2] += 2
        _f = np.zeros(shape)
        _zgrid[1:-1] = zgrid
        _zgrid[0] = zgrid[0] - dz
        _zgrid[-1] = zgrid[-1] + dz
        _f[:, :, 1:-1] = f
        _f[:, :, 0] = f[:, :, -1]
        _f[:, :, -1] = f[:, :, 0]
    else:
        _zgrid, _f = zgrid, f
    dfdz = _dfdz(_zgrid, _f)
    if periodic is True:
        dfdz = dfdz[:, :, 1:-1]
    return dfdz
