import numpy as np
from . import _trilinear


def trilinear(fgrid, boxsize, x, y, z):
    ngrid = int((len(fgrid.flatten()))**(1./3.))
    while ngrid*ngrid*ngrid != len(fgrid.flatten()):
        if ngrid*ngrid*ngrid < len(fgrid.flatten()):
            ngrid += 1
        else:
            ngrid -= 1
    npart = len(x)
    return _trilinear.trilinear(fgrid.flatten(), x, y, z, boxsize, ngrid, npart)
