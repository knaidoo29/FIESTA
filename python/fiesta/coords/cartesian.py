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
