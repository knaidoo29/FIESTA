import numpy as np
from numba import njit

from . import matrix
from . import polygon

@njit
def delaunay_area_2d(x, y, del_vert0, del_vert1, del_vert2):
    """
    Determines area for each simplex.
    
    Parameters
    ----------
    x : array
        X-coordinates.
    y : array
        Y-coordinates.
    del_vert0 : array
        Index for vertex 0 of each simplices.
    del_vert1 : array
        Index for vertex 1 of each simplices.
    del_vert2 : array
        Index for vertex 2 of each simplices.

    Returns
    -------
    areas : array
        The areas of each simplex.
    """
    nvert = len(del_vert0)
    areas = np.zeros(nvert)
    for i in range(0, nvert):
        i0 = del_vert0[i]
        i1 = del_vert1[i]
        i2 = del_vert2[i]
        areas[i] = polygon.triangle_area(x[i0], y[i0], x[i1], y[i1], x[i2], y[i2])
    return areas

@njit
def sum_delaunay_area_4_points_2d(delaunay_area, del_vert0, del_vert1, del_vert2, npart):
    """
    Finds the Delaunay area for each point.
    
    Parameters
    ----------
    delaunay_area : array
        Delaunay area.
    del_vert0 : array
        Index for vertex 0 of each simplices.
    del_vert1 : array
        Index for vertex 1 of each simplices.
    del_vert2 : array
        Index for vertex 2 of each simplices.
    npart : int
        Number of points.
    
    Returns
    -------
    point_area : array
        Delaunay area for each point.
    """
    point_area = np.zeros(npart)
    nvert = len(del_vert0)
    for i in range(0, nvert):
        i0 = del_vert0[i]
        i1 = del_vert1[i]
        i2 = del_vert2[i]
        point_area[i0] += delaunay_area[i]/3.
        point_area[i1] += delaunay_area[i]/3.
        point_area[i2]+= delaunay_area[i]/3.
    return point_area

@njit
def get_delf0_2d(x, y, f, del_vert0, del_vert1, del_vert2):
    """
    Determines delf0 for each simplices.

    Parameters
    ----------
    x : array
        X-coordinates.
    y : array
        Y-coordinates.
    f : array
        Field values at these points.
    del_vert0 : array
        Index for vertex 0 of each simplices.
    del_vert1 : array
        Index for vertex 1 of each simplices.
    del_vert2 : array
        Index for vertex 2 of each simplices.
    
    Returns
    -------
    delf0 : array
        The 2D difference in each simplices.
    """
    nvert = len(del_vert0)
    delf0 = np.zeros(2*nvert)
    m = np.zeros(4)
    for i in range(0, nvert):
        i0 = del_vert0[i]
        i1 = del_vert1[i]
        i2 = del_vert2[i]

        dx1 = x[i1] - x[i0]
        dx2 = x[i2] - x[i0]

        dy1 = y[i1] - y[i0]
        dy2 = y[i2] - y[i0]

        df1 = f[i1] - f[i0]
        df2 = f[i2] - f[i0]

        m[0] = dx1
        m[1] = dy1
        m[2] = dx2
        m[3] = dy2

        invm = matrix.inv2by2(m)

        delf0[2*i] = invm[0]*df1 + invm[1]*df2
        delf0[2*i+1] = invm[2]*df1 + invm[3]*df2
    return delf0

@njit
def delaunay_estimate_2d(simplices, x, y, x0, y0, f0, delf0):
    """
    Estimates a field from Delaunay tesselation.

    Parameters
    ----------
    simplices : array
    x : array
        X-coordinates for estimates.
    y : array
        Y-coordinates for estimates.
    x0 : array
        X-coordinate of vertex 0 of each simplices.
    y0 : array
        Y-coordinate of vertex 0 of each simplices.
    f0 : array
        Field values at vertex 0 of each simplices.
    delf0 : array
        The 2D difference in each simplices.
    
    Returns
    -------
    f_est : array
        Estimates of the field.
    """
    npart = len(x)
    f_est = np.zeros(npart)
    for i in range(0, npart):
        j = simplices[i]
        f_est[i] = f0[j] + delf0[2*j]*(x[i] - x0[j]) + delf0[2*j+1]*(y[i] - y0[j])
    return f_est