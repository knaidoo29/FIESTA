import numpy as np
from .. import src


def bilinear(fgrid, boxsize, x, y, origin=0., fill_value=np.nan, periodic=True):
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
    ngrids = np.shape(fgrid)
    if np.isscalar(boxsize):
        xbox = boxsize
        ybox = boxsize
    else:
        xbox, ybox = boxsize[0], boxsize[1]
    if np.isscalar(origin):
        _x = np.copy(x) - origin
        _y = np.copy(y) - origin
    else:
        _x = np.copy(x) - origin[0]
        _y = np.copy(y) - origin[1]
    # check if particles are inside the box
    cond = np.where((_x >= 0.) & (_x < xbox) & (_y >= 0.) & (_y < ybox))[0]
    if len(cond) == len(_x):
        # All particles are within the boundaries so no boundary management is necessary.
        npart = len(_x)
        if np.isscalar(periodic):
            if periodic == True:
                f = src.bilinear_periodic(fgrid=fgrid.flatten(), x=_x, y=_y,
                                          xbox=xbox, ybox=ybox, ngridx=ngrids[0],
                                          ngridy=ngrids[1], npart=npart)
            else:
                f = src.bilinear_nonperiodic(fgrid=fgrid.flatten(), x=_x, y=_y,
                                          xbox=xbox, ybox=ybox, ngridx=ngrids[0],
                                          ngridy=ngrids[1], npart=npart)
        else:
            if periodic[0] is True:
                perix = 1
            else:
                perix = 0
            if periodic[1] is True:
                periy = 1
            else:
                periy = 0
            f = src.bilinear_axisperiodic(fgrid=fgrid.flatten(), x=_x, y=_y,
                                          xbox=xbox, ybox=ybox, perix=perix,
                                          periy=periy, ngridx=ngrids[0],
                                          ngridy=ngrids[1], npart=npart)
    else:
        # Some particles are outside the boundary.
        # create a mask for in and outside the boxmask = np.zeros(len(x))
        mask = np.zeros(len(_x))
        # assign particles in the boundary a binary mask of 1
        mask[cond] = 1.
        # find bilinear interpolation for points inside the boundary.
        npart = len(x[cond])
        f = np.zeros(len(_x))
        if np.isscalar(periodic):
            if periodic == True:
                f[cond] = src.bilinear_periodic(fgrid=fgrid.flatten(), x=_x[cond],
                                                y=_y[cond], xbox=xbox, ybox=ybox,
                                                ngridx=ngrids[0], ngridy=ngrids[1],
                                                npart=npart)
            else:
                f[cond] = src.bilinear_nonperiodic(fgrid=fgrid.flatten(), x=_x[cond],
                                                y=_y[cond], xbox=xbox, ybox=ybox,
                                                ngridx=ngrids[0], ngridy=ngrids[1],
                                                npart=npart)
        else:
            if periodic[0] is True:
                perix = 1
            else:
                perix = 0
            if periodic[1] is True:
                periy = 1
            else:
                periy = 0
            f[cond] = src.bilinear_axisperiodic(fgrid=fgrid.flatten(), x=_x[cond],
                                                y=_y[cond], xbox=xbox, ybox=ybox,
                                                perix=perix, periy=periy, ngridx=ngrids[0],
                                                ngridy=ngrids[1], npart=npart)
        # fill outside boundary with fill values.
        cond = np.where(mask == 0.)[0]
        f[cond] = fill_value
    return f
