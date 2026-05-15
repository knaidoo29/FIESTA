import numpy as np
from numba import njit


@njit
def xgrid(xpix: int, dx: float, xmin: float) -> float:  # pragma: no cover
    """
    Returns the grid value along one axis.

    Parameters
    ----------
    xpix : int
        Index along the axis.
    dx : float
        Cell width.
    xmin : float
        Minimum x-coordinate.

    Returns
    -------
    xg : float
        Grid coordinate.
    """
    return dx / 2.0 + xpix * dx + xmin


@njit
def xgrids(xpix: np.ndarray, dx: int, xmin: int) -> np.ndarray:  # pragma: no cover
    """
    Returns the grid values along one axis for multiple indices.

    Parameters
    ----------
    xpix : array (int)
        Array of indices along the axis.
    dx : float
        Cell width.
    xmin : float
        Minimum x-coordinate.

    Returns
    -------
    xg : array (float)
        Grid coordinates.
    """
    pixlen = len(xpix)
    xg = np.zeros(pixlen, dtype=np.float64)

    for i in range(0, pixlen):
        xg[i] = xgrid(xpix[i], dx, xmin)

    return xg
