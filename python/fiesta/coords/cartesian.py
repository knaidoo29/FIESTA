import numpy as np


def xy2points(x, y):
    """Column stacks input coordinates.

    Parameters
    ----------
    x : array
        X coordinates.
    y : array
        Y coordinates.

    Return
    ------
    points : 2darray
        Column stacked array.
    """
    if np.isscalar(x) == True:
        points = np.array([[x, y]])
    else:
        points = np.column_stack((x, y))
    return points


def points2xy(points):
    """Unstacks input coordinates.

    Parameters
    ----------
    points : 2darray
        Column stacked array.

    Return
    ------
    x : array
        X coordinates.
    y : array
        Y coordinates.
    """
    x, y = points[:, 0], points[:, 1]
    return x, y


def xyz2points(x, y, z):
    """Column stacks input coordinates.

    Parameters
    ----------
    x : array
        X coordinates.
    y : array
        Y coordinates.
    z : array
        Z coordinates.

    Return
    ------
    points : 2darray
        Column stacked array.
    """
    if np.isscalar(x) == True:
        points = np.array([[x, y, z]])
    else:
        points = np.column_stack((x, y, z))
    return points


def points2xyz(points):
    """Unstacks input coordinates.

    Parameters
    ----------
    points : 2darray
        Column stacked array.

    Return
    ------
    x : array
        X coordinates.
    y : array
        Y coordinates.
    z : array
        Z coordinates.
    """
    x, y, z = points[:, 0], points[:, 1], points[:, 2]
    return x, y, z
