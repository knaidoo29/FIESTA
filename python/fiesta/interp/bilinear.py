import numpy as np
from .. import src


def bilinear(fgrid, boxsize, x, y, fill_value=np.nan, periodic=True):
    """ Bilinear interpolation from a 2D grid defined in box of [0., boxsize].

    Parameter
    ---------
    fgrid : array
        Field values on a 2D grid.
    boxsize : float
        Box size.
    x : array
        x coordinate values.
    y : array
        y coordinate values.
    fill_value : float, optional
        Fill outside boundary values.
    periodic : bool, optional
        Determines whether to interpolate on a periodic grid.

    Returns
    -------
    f : array
        Field interpolation values.
    """
    # determine ngrid from fgrid
    ngrid = int(np.sqrt(len(fgrid.flatten())))
    # correct ngrid if its slightly off
    while ngrid*ngrid != len(fgrid.flatten()):
        if ngrid*ngrid < len(fgrid.flatten()):
            ngrid += 1
        else:
            ngrid -= 1
    # check if particles are inside the box
    condition = np.where((x >= 0.) & (x < boxsize) & (y >= 0.) & (y < boxsize))[0]
    if len(condition) == len(x):
        # All particles are within the boundaries so no boundary management is necessary.
        npart = len(x)
        if periodic == True:
            f = src.bilinear_periodic(fgrid.flatten(), x[condition], y[condition], boxsize, ngrid, npart)
        else:
            f = src.bilinear_nonperiodic(fgrid.flatten(), x[condition], y[condition], boxsize, ngrid, npart)
    else:
        # Some particles are outside the boundary.
        # create a mask for in and outside the boxmask = np.zeros(len(x))
        mask = np.zeros(len(x))
        # assign particles in the boundary a binary mask of 1
        mask[condition] = 1.
        # find bilinear interpolation for points inside the boundary.
        npart = len(x[condition])
        f = np.zeros(len(x))
        if periodic == True:
            f[condition] = src.bilinear_periodic(fgrid.flatten(), x[condition], y[condition], boxsize, ngrid, npart)
        else:
            f[condition] = src.bilinear_nonperiodic(fgrid.flatten(), x[condition], y[condition], boxsize, ngrid, npart)
        # fill outside boundary with fill values.
        condition = np.where(mask == 0.)[0]
        f[condition] = fill_value
    return f
