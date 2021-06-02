import numpy as np
from .. import src


def part2grid2d(x, y, f, boxsize, ngrid, method='TSC', periodic=True):
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
    periodic : bool, optional
        Assign particles with periodic boundaries.

    Returns
    -------
    fgrid : array
        Grid assigned values.
    """
    if method == 'NGP':
        fgrid = src.p2g_ngp_2d(x, y, f, boxsize, ngrid, len(x))
    elif method == 'CIC':
        if periodic == True:
            fgrid = src.p2g_cic_2d_periodic(x, y, f, boxsize, ngrid, len(x))
        else:
            fgrid = src.p2g_cic_2d_nonperiodic(x, y, f, boxsize, ngrid, len(x))
    elif method == 'TSC':
        if periodic == True:
            fgrid = src.p2g_tsc_2d_periodic(x, y, f, boxsize, ngrid, len(x))
        else:
            fgrid = src.p2g_tsc_2d_nonperiodic(x, y, f, boxsize, ngrid, len(x))
    return fgrid


def part2grid3d(x, y, z, f, boxsize, ngrid, method='TSC', periodic=True):
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
    periodic : bool, optional
        Assign particles with periodic boundaries.

    Returns
    -------
    fgrid : array
        Grid assigned values.
    """
    if method == 'NGP':
        fgrid = src.p2g_ngp_3d(x, y, z, f, boxsize, ngrid, len(x))
    elif method == 'CIC':
        if periodic == True:
            fgrid = src.p2g_cic_3d_periodic(x, y, z, f, boxsize, ngrid, len(x))
        else:
            fgrid = src.p2g_cic_3d_nonperiodic(x, y, z, f, boxsize, ngrid, len(x))
    elif method == 'TSC':
        if periodic == True:
            fgrid = src.p2g_tsc_3d_periodic(x, y, z, f, boxsize, ngrid, len(x))
        else:
            fgrid = src.p2g_tsc_3d_nonperiodic(x, y, z, f, boxsize, ngrid, len(x))
    return fgrid
