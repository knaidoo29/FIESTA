import numpy as np
from numba import njit

from . import grid


@njit
def which_pix(x: float, dx: float, xmin: float) -> int:  # pragma: no cover
    """
    Find pixel along a defined grid.

    Parameters
    ----------
    x : float
        A point for which we would like to determine the pixel index.
    dx : float
        Pixel width.
    xmin : float
        Minimum along the grid.

    Returns
    -------
    pix: int
        The pixel the point corresponds to.
    """
    pix = int(np.floor((x - xmin) / dx))
    return pix


@njit
def which_pixs(x: np.ndarray, dx: float, xmin: float) -> np.ndarray:  # pragma: no cover
    """
    Find pixels along a defined grid for an array of points.

    Parameters
    ----------
    x : array
        Points for which we want to determine the pixel indices.
    dx : float
        Pixel width.
    xmin : float
        Minimum along the grid.

    Returns
    -------
    pix: ndarray
        Pixel indices corresponding to the points.
    """
    pix = np.floor((x - xmin) / dx).astype(np.int32)
    return pix


@njit
def pix1dto2d(
    xpix: np.ndarray, ypix: np.ndarray, ygrid: int
) -> np.ndarray:  # pragma: no cover
    """
    Maps pixels given along a single axis in x and y onto a flattened 2D grid.

    Parameters
    ----------
    xpix : int array
        Pixel indices along the x-axis.
    ypix : int array
        Pixel indices along the y-axis.
    ygrid : int
        Length of the y-axis grid.

    Returns
    -------
    pix: ndarray
        Flattened 2D grid pixel indices.
    """
    xlen, ylen = len(xpix), len(ypix)
    pix = np.full(xlen * ylen, -1, dtype=np.int32)
    idx = 0
    for i in range(xlen):
        for j in range(ylen):
            if xpix[i] != -1 and ypix[j] != -1:
                pix[idx] = ypix[j] + ygrid * xpix[i]
            idx += 1
    return pix


@njit
def pix1dto3d(
    xpix: np.ndarray, ypix: np.ndarray, zpix: np.ndarray, ygrid: int, zgrid: int
) -> np.ndarray:  # pragma: no cover
    """
    Maps pixels given along a single axis in x, y, and z onto a flattened 3D grid.

    Parameters
    ----------
    xpix : int array
        Pixel indices along the x-axis.
    ypix : int array
        Pixel indices along the y-axis.
    zpix : int array
        Pixel indices along the z-axis.
    ygrid : int
        Length of the y-axis grid.
    zgrid : int
        Length of the z-axis grid.

    Returns
    -------
    pix: ndarray
        Flattened 3D grid pixel indices.
    """
    xlen, ylen, zlen = len(xpix), len(ypix), len(zpix)
    pix = np.full(xlen * ylen * zlen, -1, dtype=np.int32)
    idx = 0
    for i in range(xlen):
        for j in range(ylen):
            for k in range(zlen):
                if xpix[i] != -1 and ypix[j] != -1 and zpix[k] != -1:
                    pix[idx] = zpix[k] + zgrid * (ypix[j] + ygrid * xpix[i])
                idx += 1
    return pix


@njit
def pix1dto2d_scalar(xpix: int, ypix: int, ygrid: int) -> int:  # pragma: no cover
    """
    Maps pixels given along a single axis in x and y onto a flattened 2D grid.

    Parameters
    ----------
    xpix : int
        Pixel indices along the x-axis.
    ypix : int
        Pixel indices along the y-axis.
    ygrid : int
        Length of the y-axis grid.

    Returns
    -------
    pix : int
        Flattened 2D grid pixel indices.
    """
    if xpix != -1 and ypix != -1:
        pix = ypix + ygrid * xpix
    else:
        pix = -1
    return pix


@njit
def pix1dto3d_scalar(
    xpix: np.ndarray, ypix: np.ndarray, zpix: np.ndarray, ygrid: int, zgrid: int
) -> np.ndarray:  # pragma: no cover
    """
    Maps pixels given along a single axis in x, y, and z onto a flattened 3D grid.

    Parameters
    ----------
    xpix : int array
        Pixel indices along the x-axis.
    ypix : int array
        Pixel indices along the y-axis.
    zpix : int array
        Pixel indices along the z-axis.
    ygrid : int
        Length of the y-axis grid.
    zgrid : int
        Length of the z-axis grid.

    Returns
    -------
    pix : int
        Flattened 3D grid pixel indices.
    """
    if xpix != -1 and ypix != -1 and zpix != -1:
        pix = zpix + zgrid * (ypix + ygrid * xpix)
    else:
        pix = -1
    return pix


@njit
def periodic_pix(pix: np.ndarray, ngrid: int) -> np.ndarray:  # pragma: no cover
    """
    Applies periodic boundary conditions to pixel indices.

    Parameters
    ----------
    pix : ndarray
        Pixel index array.
    ngrid : int
        Grid dimensions.

    Returns
    -------
    ndarray
        Periodic pixel index array.
    """
    for i in range(len(pix)):
        if pix[i] < 0:
            pix[i] += ngrid
        elif pix[i] >= ngrid:
            pix[i] -= ngrid
    return pix


@njit
def ngp_pix(x: float, dx: float, xmin: float) -> int:  # pragma: no cover
    """
    Nearest-grid-point pixel index.

    Parameters
    ----------
    x : float
        X-coordinate.
    dx : float
        Grid spacing.
    xmin : float
        Minimum grid value.

    Returns
    -------
    pix: int
        Nearest grid-point pixel index.
    """
    pix = which_pix(x, dx, xmin)
    return pix


@njit
def cic_pix(x: float, dx: float, xmin: float) -> np.ndarray:  # pragma: no cover
    """
    Cloud-in-cell pixel indices.

    Parameters
    ----------
    x : float
        X-coordinate.
    dx : float
        Grid spacing.
    xmin : float
        Minimum grid value.

    Returns
    -------
    idx: array
        Two closest pixel indices.
    """
    xpix = which_pix(x, dx, xmin)
    xg = grid.xgrid(xpix, dx, xmin)
    if x < xg:
        idx = np.array([xpix - 1, xpix])

    else:
        idx = np.array([xpix, xpix + 1])
    return idx


@njit
def tsc_pix(x: float, dx: float, xmin: float) -> np.ndarray:  # pragma: no cover
    """
    Triangular-shaped-cloud pixel indices.

    Parameters
    ----------
    x : float
        X-coordinate.
    dx : float
        Grid spacing.
    xmin : float
        Minimum grid value.

    Returns
    -------
    idx: array
        Three closest pixel indices.
    """
    xpix = which_pix(x, dx, xmin)
    idx = np.array([xpix - 1, xpix, xpix + 1])
    return idx


@njit
def pcs_pix(x: float, dx: float, xmin: float) -> np.ndarray:  # pragma: no cover
    """
    Piecewise-Cubic-Spline pixel indices.

    Parameters
    ----------
    x : float
        X-coordinate.
    dx : float
        Grid spacing.
    xmin : float
        Minimum grid value.

    Returns
    -------
    idx : array
        Four closest pixel indices.
    """
    xpix = which_pix(x, dx, xmin)
    xg = grid.xgrid(xpix, dx, xmin)
    if x < xg:
        idx = np.array([xpix - 2, xpix - 1, xpix, xpix + 1])
    else:
        idx = np.array([xpix - 1, xpix, xpix + 1, xpix + 2])
    return idx
