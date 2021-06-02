import numpy as np
from . import _bilinear


def bilinear(fgrid, boxsize, x, y):
    ngrid = int(np.sqrt(len(fgrid.flatten())))
    npart = len(x)
    return _bilinear.bilinear(fgrid.flatten(), x, y, boxsize, ngrid, npart)
