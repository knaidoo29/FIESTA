import numpy as np
from numba import njit


@njit
def triangle_area(xa, ya, xb, yb, xc, yc):
    """
    Determines the area of a triangle given its vertex coordinates.
    """
    return 0.5 * abs(xa * (yb - yc) + xb * (yc - ya) + xc * (ya - yb))


@njit
def sum_triangle_area(xas, yas, xbs, ybs, xcs, ycs):
    """
    Computes the total area of multiple triangles.
    """
    ntri = len(xas)
    total_area = 0.0
    
    for i in range(ntri):
        total_area += triangle_area(xas[i], yas[i], xbs[i], ybs[i], xcs[i], ycs[i])
    
    return total_area
