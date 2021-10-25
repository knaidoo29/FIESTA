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
        fgrid = src.part2grid_ngp_2d(x=x, y=y, f=f, boxsize=boxsize, npart=len(x), ngrid=ngrid)
    elif method == 'CIC':
        if periodic == True:
            fgrid = src.part2grid_cic_2d(x=x, y=y, f=f, boxsize=boxsize, npart=len(x), ngrid=ngrid, periodic=True)
        else:
            fgrid = src.part2grid_cic_2d(x=x, y=y, f=f, boxsize=boxsize, npart=len(x), ngrid=ngrid, periodic=False)
    elif method == 'TSC':
        if periodic == True:
            fgrid = src.part2grid_tsc_2d(x=x, y=y, f=f, boxsize=boxsize, npart=len(x), ngrid=ngrid, periodic=True)
        else:
            fgrid = src.part2grid_tsc_2d(x=x, y=y, f=f, boxsize=boxsize, npart=len(x), ngrid=ngrid, periodic=False)
    return fgrid.reshape(ngrid, ngrid)


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
        fgrid = src.part2grid_ngp_3d(x=x, y=y, z=z, f=f, boxsize=boxsize, npart=len(x), ngrid=ngrid)
    elif method == 'CIC':
        if periodic == True:
            fgrid = src.part2grid_cic_3d(x=x, y=y, z=z, f=f, boxsize=boxsize, npart=len(x), ngrid=ngrid, periodic=True)
        else:
            fgrid = src.part2grid_cic_3d(x=x, y=y, z=z, f=f, boxsize=boxsize, npart=len(x), ngrid=ngrid, periodic=False)
    elif method == 'TSC':
        if periodic == True:
            fgrid = src.part2grid_tsc_3d(x=x, y=y, z=z, f=f, boxsize=boxsize, npart=len(x), ngrid=ngrid, periodic=True)
        else:
            fgrid = src.part2grid_tsc_3d(x=x, y=y, z=z, f=f, boxsize=boxsize, npart=len(x), ngrid=ngrid, periodic=False)
    return fgrid.reshape(ngrid, ngrid, ngrid)
