import numpy as np


def grid2d(boxsize, ngrid):
    """Returns the x, y, z coordinates of a cartesian grid.

    Parameters
    ----------
    boxsize : float
        Box size.
    ngrid : int
        Grid division along one axis.

    Returns
    -------
    x2d : array
        X coordinates on a 2D cartesian grid.
    y2d : array
        Y coordinates on a 2D cartesian grid.
    """
    xedges = np.linspace(0., boxsize, ngrid + 1)
    x = 0.5*(xedges[1:] + xedges[:-1])
    x2d, y2d = np.meshgrid(x, x, indexing='ij')
    return x2d, y2d


def grid3d(boxsize, ngrid):
    """Returns the x, y, z coordinates of a cartesian grid.

    Parameters
    ----------
    boxsize : float
        Box size.
    ngrid : int
        Grid division along one axis.

    Returns
    -------
    x3d : array
        X coordinates on a 3D cartesian grid.
    y3d : array
        Y coordinates on a 3D cartesian grid.
    z3d : array
        Z coordinates on a 3D cartesian grid.
    """
    xedges = np.linspace(0., boxsize, ngrid + 1)
    x = 0.5*(xedges[1:] + xedges[:-1])
    x3d, y3d, z3d = np.meshgrid(x, x, x, indexing='ij')
    return x3d, y3d, z3d
