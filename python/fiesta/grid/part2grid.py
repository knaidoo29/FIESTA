import numpy as np
from .. import src


def part2grid2d(x, y, f, boxsize, ngrid, method='TSC'):
    """Returns the density contrast for the nearest grid point grid assignment.

    Parameters
    ----------
    x : array
        X coordinates of the particle.
    y : array
        Y coordinates of the particle.
    f : array
        Value of each particle to be assigned to the grid.
    boxsize : float
        Box size.
    ngrid : int
        Grid divisions across one axis.
    method : str, optional
        Grid assignment scheme, either 'NGP', 'CIC' or 'TSC'.

    Returns
    -------
    grid_value : array
        Grid assigned values.
    """
    if method == 'NGP':
        grid_value = src.p2g_ngp_2d(x, y, f, boxsize, ngrid, len(x))
    elif method == 'CIC':
        grid_value = src.p2g_cic_2d(x, y, f, boxsize, ngrid, len(x))
    elif method == 'TSC':
        grid_value = src.p2g_tsc_2d(x, y, f, boxsize, ngrid, len(x))
    return grid_value


def part2grid3d(x, y, z, f, boxsize, ngrid, method='TSC'):
    """Returns the density contrast for the nearest grid point grid assignment.

    Parameters
    ----------
    x : array
        X coordinates of the particle.
    y : array
        Y coordinates of the particle.
    z : array
        Z coordinates of the particle.
    f : array
        Value of each particle to be assigned to the grid.
    boxsize : float
        Box size.
    ngrid : int
        Grid divisions across one axis.
    method : str, optional
        Grid assignment scheme, either 'NGP', 'CIC' or 'TSC'.

    Returns
    -------
    grid_value : array
        Grid assigned values.
    """
    if method == 'NGP':
        grid_value = src.p2g_ngp_3d(x, y, z, f, boxsize, ngrid, len(x))
    elif method == 'CIC':
        grid_value = src.p2g_cic_3d(x, y, z, f, boxsize, ngrid, len(x))
    elif method == 'TSC':
        grid_value = src.p2g_tsc_3d(x, y, z, f, boxsize, ngrid, len(x))
    return grid_value
