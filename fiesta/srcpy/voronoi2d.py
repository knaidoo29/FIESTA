import numpy as np
from numba import njit

from . import polygon


@njit
def voronoi_2d_area(xpoints, ypoints, xverts, yverts, ridge_point1, ridge_point2, ridge_vertices, ridge_start, ridge_end):
    """
    Determines the area of the voronoi cells.

    Parameters
    ----------
    xpoints : array
        X-coordinates for the points.
    ypoints : array
        Y-coordinates for the points.
    xverts : array
        X-coordinates for the voronoi vertices.
    yverts : array
        Y-coordinates for the voronoi vertices.
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
    area : array
        The area of each voronoi cell.
    """
    npoints = len(xpoints)
    # nvertices = len(xverts)
    nridge = len(ridge_point1)
    # nridge_vertices = len(ridge_vertices)
    
    area = np.zeros(npoints)
    
    # computes the area of each voronoi cell by breaking it into triangles between points and ridges.

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

            xa2 = xpoints[ridge_point2[i]]
            ya2 = ypoints[ridge_point2[i]]

            for j in range(ridge_start[i], ridge_end[i]):

                xb = xverts[ridge_vertices[j]]
                yb = yverts[ridge_vertices[j]]

                xc = xverts[ridge_vertices[j+1]]
                yc = yverts[ridge_vertices[j+1]]

                area1 = polygon.triangle_area(xa1, ya1, xb, yb, xc, yc)
                area2 = polygon.triangle_area(xa2, ya2, xb, yb, xc, yc)

                if area[ridge_point1[i]] != -1.:
                    area[ridge_point1[i]] = area[ridge_point1[i]] + area1

                if area[ridge_point2[i]] != -1.:
                    area[ridge_point2[i]] = area[ridge_point2[i]]+ area2

        else:
            area[ridge_point1[i]] = -1.
            area[ridge_point2[i]] = -1.
    
    return area
