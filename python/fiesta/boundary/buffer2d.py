import numpy as np

from .. import randoms


def buffer_random_particles_2d(npart, boxsize, buffer_length):
    """Generates random buffer particles around a 2D box.

    Parameters
    ----------
    npart : int
        Number of particles in the 2D box.
    boxsize : float
        Box size.
    buffer_length : float
        Length of the buffer region.

    Returns
    -------
    x : array
        Random x-values.
    y : array
        Random y-values.
    """

    # Determine number of particles needed in the padded region to maintain
    # the average density.

    part_dens = npart / (boxsize ** 2.)
    buffer_area = ((boxsize + 2.*buffer_length)**2.) - (boxsize**2.)
    npart_buffer = part_dens * buffer_area

    # Divide this into four regions.

    size_quarter = int(npart_buffer/4.)

    # box edge 1
    _xmin = 0. - buffer_length
    _xmax = boxsize
    _ymin = 0. - buffer_length
    _ymax = 0.
    x1, y1 = randoms.random_box(size_quarter, _xmin, _xmax, _ymin, _ymax)

    # box edge 2
    _xmin = boxsize
    _xmax = boxsize + buffer_length
    _ymin = 0. - buffer_length
    _ymax = boxsize
    x2, y2 = randoms.random_box(size_quarter, _xmin, _xmax, _ymin, _ymax)

    # box edge 3
    _xmin = 0.
    _xmax = boxsize + buffer_length
    _ymin = boxsize
    _ymax = boxsize + buffer_length
    x3, y3 = randoms.random_box(size_quarter, _xmin, _xmax, _ymin, _ymax)

    # box edge 4
    _xmin = 0. - buffer_length
    _xmax = 0.
    _ymin = 0.
    _ymax = boxsize + buffer_length
    x4, y4 = randoms.random_box(size_quarter, _xmin, _xmax, _ymin, _ymax)

    x = np.concatenate([x1, x2, x3, x4])
    y = np.concatenate([y1, y2, y3, y4])

    return x, y


def buffer_periodic_particles_2d(x, y, boxsize, buffer_length):
    """Generates random buffer particles around a 2D box.

    Parameters
    ----------
    x : array
        X-coordinates.
    y : array
        Y-coordinates.
    boxsize : float
        Box size.
    buffer_length : float
        Length of the buffer region.

    Returns
    -------
    xp : array
        Periodic x-values.
    yp : array
        Periodic y-values.
    """
    assert buffer_length < boxsize, "buffer_length must be smaller than the boxsize."
    ix = np.array([-1, 0, 1])
    ixs, iys = np.meshgrid(ix, ix)
    ixs = ixs.flatten()
    iys = iys.flatten()
    for i in range(0, len(ixs)):
        ix, iy = ixs[i], iys[i]
        if ix != 0 or iy != 0:
            xnew = np.copy(x) + ix*boxsize
            ynew = np.copy(y) + iy*boxsize
            cond = np.where((xnew >= -buffer_length) & (xnew <= boxsize+buffer_length) &
                            (ynew >= -buffer_length) & (ynew <= boxsize+buffer_length))[0]
            if ix == -1 and iy == -1:
                xp = xnew[cond]
                yp = ynew[cond]
            else:
                xp = np.concatenate([xp, xnew[cond]])
                yp = np.concatenate([yp, ynew[cond]])
        else:
            pass
    return xp, yp
