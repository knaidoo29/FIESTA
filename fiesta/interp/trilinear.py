import numpy as np
from .. import src


def trilinear(fgrid, boxsize, x, y, z, fill_value=np.nan, periodic=True):
    """ Trilinear interpolation from a 3D grid defined in box of [0., boxsize].

    Parameter
    ---------
    fgrid : array
        Field values on a 3D grid.
    boxsize : float
        Box size.
    x : array
        x coordinate values.
    y : array
        y coordinate values.
    z : array
        z coordinate values.
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
    ngrid = int((len(fgrid.flatten()))**(1./3.))
    # correct ngrid if its slightly off
    while ngrid*ngrid*ngrid != len(fgrid.flatten()):
        if ngrid*ngrid*ngrid < len(fgrid.flatten()):
            ngrid += 1
        else:
            ngrid -= 1
    # check if particles are inside the box
    condition = np.where((x >= 0.) & (x < boxsize) & (y >= 0.) & (y < boxsize) & (z >= 0.) & (z < boxsize))[0]
    if len(condition) == len(x):
        # All particles are within the boundaries so no boundary management is necessary.
        npart = len(x)
        if periodic == True:
            f = src.trilinear_periodic(fgrid=fgrid.flatten(), x=x, y=y, z=z, boxsize=boxsize, ngrid=ngrid, npart=npart)
        else:
            f = src.trilinear_nonperiodic(fgrid=fgrid.flatten(), x=x, y=y, z=z, boxsize=boxsize, ngrid=ngrid, npart=npart)
    else:
        # Some particles are outside the boundary.
        # create a mask for in and outside the box
        mask = np.zeros(len(x))
        # assign particles in the boundary a binary mask of 1.
        mask[condition] = 1.
        # find trilinear interpolation for points inside the boundary.
        npart = len(x[condition])
        f = np.zeros(len(x))
        if periodic == True:
            f[condition] = src.trilinear_periodic(fgrid=fgrid.flatten(), x=x[condition], y=y[condition], z=z[condition], boxsize=boxsize, ngrid=ngrid, npart=npart)
        else:
            f[condition] = src.trilinear_nonperiodic(fgrid=fgrid.flatten(), x=x[condition], y=y[condition], z=z[condition], boxsize=boxsize, ngrid=ngrid, npart=npart)
        # fill outside boundary with fill values.
        condition = np.where(mask == 0.)[0]
        f[condition] = fill_value
    return f
