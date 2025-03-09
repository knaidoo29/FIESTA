import numpy as np
from numba import njit

from . import polyhedron 


@njit 
def voronoi_3d_volume(xpoints, ypoints, zpoints, xverts, yverts, zverts, ridge_point1, ridge_point2, ridge_vertices, ridge_start, ridge_end):
    """
    Determines the volume of the voronoi cells.
    
    Parameters
    ----------
    xpoints : array
        X-coordinates for the points.
    ypoints : array
        Y-coordinates for the points.
    zpoints : array
        Z-coordinates for the points.
    xverts : array
        X-coordinates for the voronoi vertices.
    yverts : array
        Y-coordinates for the voronoi vertices.
    zverts : array
        Z-coordinates for the voronoi vertices.
    ridge_points1 : array
        The points on either side of a voronoi ridge.
    ridge_points2 : array
        The points on either side of a voronoi ridge.
    ridge_vertices : array
        The index of the vertices on the ridge.
    ridge_start : array
        Index in ridge vertex which specifies the start of vertices that make up a single ridge.
    ridge_end : array
        Index in ridge vertex which specifies the end of vertices that make up a single ridge.
    
    Returns
    -------
    volume : array
        The volume of each voronoi cell.
    """

    npoints = len(xpoints)
    # nvertices = len(xverts)
    nridge = len(ridge_point1)
    # nridge_vertices = len(ridge_vertices)
    
    volume = np.zeros(npoints)

    # computes the volume of each voronoi cell by breaking it into tetrahedron between points and ridges.

    for i in range(0, nridge):
        check = 1
        for j in range(ridge_start[i], ridge_end[i]):
            if ridge_vertices[j] == -1:
                check = 0
            if ridge_vertices[j+1] == -1:
                check = 0

        if check == 1:

            xa1 = xpoints[ridge_point1[i]]
            ya1 = ypoints[ridge_point1[i]]
            za1 = zpoints[ridge_point1[i]]

            xa2 = xpoints[ridge_point2[i]]
            ya2 = ypoints[ridge_point2[i]]
            za2 = zpoints[ridge_point2[i]]

            xb = xverts[ridge_vertices[ridge_start[i]]]
            yb = yverts[ridge_vertices[ridge_start[i]]]
            zb = zverts[ridge_vertices[ridge_start[i]]]
            
            for j in range(ridge_start[i], ridge_end[i]-1):

                xc = xverts[ridge_vertices[j+1]]
                yc = yverts[ridge_vertices[j+1]]
                zc = zverts[ridge_vertices[j+1]]

                xd = xverts[ridge_vertices[j+2]]
                yd = yverts[ridge_vertices[j+2]]
                zd = zverts[ridge_vertices[j+2]]

                vol1 = polyhedron.tetrahedron_volume(xa1, ya1, za1, xb, yb, zb, xc, yc, zc, xd, yd, zd)
                vol2 = polyhedron.tetrahedron_volume(xa2, ya2, za2, xb, yb, zb, xc, yc, zc, xd, yd, zd)

                if volume[ridge_point1[i]] != -1.:
                    volume[ridge_point1[i]] = volume[ridge_point1[i]] + vol1

                if volume[ridge_point2[i]] != -1.:
                    volume[ridge_point2[i]] = volume[ridge_point2[i]] + vol2
        else:
            volume[ridge_point1[i]] = -1.
            volume[ridge_point2[i]] = -1.
    return volume
