import numpy as np


def grid1D(boxsize, ngrid):
    """Returns the x coordinates of a cartesian grid.

    Parameters
    ----------
    boxsize : float
        Box size.
    ngrid : int
        Grid division along one axis.

    Returns
    -------
    xedges : array
        x coordinate bin edges.
    x : array
        X coordinates bin centers.
    """
    xedges = np.linspace(0., boxsize, ngrid + 1)
    x = 0.5*(xedges[1:] + xedges[:-1])
    return xedges, x


def grid2D(boxsize, ngrid):
    """Returns the x, y coordinates of a cartesian grid.

    Parameters
    ----------
    boxsize : float
        Box size.
    ngrid : int
        Grid division along one axis.

    Returns
    -------
    x2D : array
        X coordinates on a 2D cartesian grid.
    y2D : array
        Y coordinates on a 2D cartesian grid.
    """
    xedges, x = grid1D(boxsize, ngrid)
    x2D, y2D = np.meshgrid(x, x, indexing='ij')
    return x2D, y2D


def grid3D(boxsize, ngrid):
    """Returns the x, y, z coordinates of a cartesian grid.

    Parameters
    ----------
    boxsize : float
        Box size.
    ngrid : int
        Grid division along one axis.

    Returns
    -------
    x3D : array
        X coordinates on a 3D cartesian grid.
    y3D : array
        Y coordinates on a 3D cartesian grid.
    z3D : array
        Z coordinates on a 3D cartesian grid.
    """
    xedges, x = grid1D(boxsize, ngrid)
    x3D, y3D, z3D = np.meshgrid(x, x, x, indexing='ij')
    return x3D, y3D, z3D
