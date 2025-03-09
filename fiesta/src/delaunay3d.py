import numpy as np
from numba import njit

from . import matrix
from . import polyhedron


@njit
def delaunay_volume_3d(x, y, z, del_vert0, del_vert1, del_vert2, del_vert3):
    """
    Determines the volume for each simplex.

    Parameters
    ----------
    x : array
        X-coordinates.
    y : array
        Y-coordinates.
    z : array
        Z-coordinates.
    del_vert0 : array
        Index for vertex 0 of each simplices.
    del_vert1 : array
        Index for vertex 1 of each simplices.
    del_vert2 : array
        Index for vertex 2 of each simplices.
    del_vert3 : array
        Index for vertex 3 of each simplices.
    npart : int
        Number of points.
    nvert : int
        Number of vertices.
    
    Returns
    -------
    volumes : array
        The volumes of each simplex.
    """
    nvert = len(del_vert0)
    volumes = np.zeros(nvert)
    for i in range(0, nvert):
        i0 = del_vert0[i]
        i1 = del_vert1[i]
        i2 = del_vert2[i]
        i3 = del_vert3[i]
        volumes[i] = polyhedron.tetrahedron_volume(x[i0], y[i0], z[i0], x[i1], y[i1], z[i1], x[i2], y[i2], z[i2], x[i3], y[i3], z[i3])
    return volumes


@njit
def sum_delaunay_vol_4_points_3d(delaunay_vol, del_vert0, del_vert1, del_vert2, del_vert3, npart):
    """
    Finds the Delaunay volume for each point.
    
    Parameters
    ----------
    delaunay_vol : array
        Delaunay volumes.
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
    point_vol : array
        Sum of delaunay value for each point.
    """
    nvert = len(del_vert0)
    point_vol = np.zeros(npart)
    for i in range(0, nvert):
        i0 = del_vert0[i]
        i1 = del_vert1[i]
        i2 = del_vert2[i]
        i3 = del_vert3[i]
        point_vol[i0] = point_vol[i0] + delaunay_vol[i]/4.
        point_vol[i1] = point_vol[i1] + delaunay_vol[i]/4.
        point_vol[i2] = point_vol[i2] + delaunay_vol[i]/4.
        point_vol[i3] = point_vol[i3] + delaunay_vol[i]/4.
    return point_vol


@njit
def get_delf0_3d(x, y, z, f, del_vert0, del_vert1, del_vert2, del_vert3):
    """
    Determines delf0 for each simplices.

    Parameters
    ----------
    x : array
        X-coordinates.
    y : array
        Y-coordinates.
    z : array
        Z-coordinates.
    f : array
        Field values at these points.
    del_vert0 : array
        Index for vertex 0 of each simplices.
    del_vert1 : array
        Index for vertex 1 of each simplices.
    del_vert2 : array
        Index for vertex 2 of each simplices.
    del_vert3 : array
        Index for vertex 3 of each simplices.
    
    Returns
    -------
    delf0 : array
        The 3D difference in each simplices.
    """
    nvert = len(del_vert0)
    delf0 = np.zeros(3*nvert)
    m = np.zeros(9)
    for i in range(0, nvert):
        i0 = del_vert0[i]
        i1 = del_vert1[i]
        i2 = del_vert2[i]
        i3 = del_vert3[i]

        dx1 = x[i1] - x[i0]
        dx2 = x[i2] - x[i0]
        dx3 = x[i3] - x[i0]

        dy1 = y[i1] - y[i0]
        dy2 = y[i2] - y[i0]
        dy3 = y[i3] - y[i0]

        dz1 = z[i1] - z[i0]
        dz2 = z[i2] - z[i0]
        dz3 = z[i3] - z[i0]

        df1 = f[i1] - f[i0]
        df2 = f[i2] - f[i0]
        df3 = f[i3] - f[i0]

        m[0] = dx1
        m[1] = dy1
        m[2] = dz1

        m[3] = dx2
        m[4] = dy2
        m[5] = dz2

        m[6] = dx3
        m[7] = dy3
        m[8] = dz3

        invm = matrix.inv3by3(m)

        delf0[3*i]   = invm[0]*df1 + invm[1]*df2 + invm[2]*df3
        delf0[3*i+1] = invm[3]*df1 + invm[4]*df2 + invm[5]*df3
        delf0[3*i+2] = invm[6]*df1 + invm[7]*df2 + invm[8]*df3

    return delf0


@njit
def delaunay_estimate_3d(simplices, x, y, z, x0, y0, z0, f0, delf0):
    """
    Estimates a field from Delaunay tesselation.

    Parameters
    ----------
    x : array
        X-coordinates for estimates.
    y : array
        Y-coordinates for estimates.
    z : array
        Z-coordinates for estimates.
    x0 : array
        X-coordinate of vertex 0 of each simplices.
    y0 : array
        Y-coordinate of vertex 0 of each simplices.
    z0 : array
        Z-coordinate of vertex 0 of each simplices.
    f0 : array
        Field values at vertex 0 of each simplices.
    delf0 : array
        The 3D difference in each simplices.
    
    Returns
    -------
    f_est : array
        Estimates of the field.
    """
    npart = len(x)
    f_est = np.zeros(npart)
    for i in range(0, npart):
        j = simplices[i]
        f_est[i] = f0[j] + delf0[3*j]*(x[i] - x0[j]) + delf0[3*j+1]*(y[i] - y0[j]) + delf0[3*j+2]*(z[i] - z0[j])
    return f_est