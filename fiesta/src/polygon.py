import numpy as np
from numba import njit


@njit
def triangle_area(xa: float, ya: float, xb: float, yb: float, xc: float, yc: float) -> float:
    """
    Determines the area of a triangle given its vertex coordinates.

    Parameters
    ----------
    xa : float
        X-axis coordinate for point A of the triangle.
    ya : float
        Y-axis coordinate for point A of the triangle.
    xb : float
        X-axis coordinate for point B of the triangle.
    yb : float
        Y-axis coordinate for point B of the triangle.
    xc : float
        X-axis coordinate for point C of the triangle.
    yc : float
        Y-axis coordinate for point C of the triangle.
    
    Yields
    ------
    Area of the triangle.
    """
    return 0.5 * abs(xa * (yb - yc) + xb * (yc - ya) + xc * (ya - yb))


@njit
def sum_triangle_area(
    xas: np.ndarray, yas: np.ndarray, xbs: np.ndarray, ybs: np.ndarray, xcs: np.ndarray, ycs: np.ndarray
) -> np.ndarray:
    """
    Computes the total area of multiple triangles.

    Parameters
    ----------
    xas : array
        X-axis coordinates for point A of the triangles.
    yas : array
        Y-axis coordinates for point A of the triangles.
    xbs : array
        X-axis coordinates for point B of the triangles.
    ybs : array
        Y-axis coordinates for point B of the triangles.
    xcs : array
        X-axis coordinates for point C of the triangles.
    ycs : array
        Y-axis coordinates for point C of the triangles.
    
    Returns
    -------
    total_area : array
        Area for each triangle.
    """
    ntri = len(xas)
    total_area = 0.0
    
    for i in range(ntri):
        total_area += triangle_area(xas[i], yas[i], xbs[i], ybs[i], xcs[i], ycs[i])
    
    return total_area
