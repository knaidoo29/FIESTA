import numpy as np

from .. import randoms


def buffer_random_2D(npart, boxsize, buffer_length, origin=0.):
    """Generates random buffer particles around a 2D box.

    Parameters
    ----------
    npart : int
        Number of particles in the 2D box.
    boxsize : float
        Box size.
    buffer_length : float
        Length of the buffer region.
    origin : float, optional
        Origin point.

    Returns
    -------
    x : array
        Random x-values.
    y : array
        Random y-values.
    """
    if np.isscalar(origin):
        xorigin, yorigin = origin, origin
    else:
        xorigin, yorigin = origin[0], origin[1]
    # Determine number of particles needed in the padded region to maintain
    # the average density.
    part_dens = npart / (boxsize ** 2.)
    buffer_area = ((boxsize + 2.*buffer_length)**2.) - (boxsize**2.)
    npart_buffer = part_dens * buffer_area
    # Divide this into four regions.
    size_quarter = int(npart_buffer/4.)
    # box edge 1
    _xmin = xorigin - buffer_length
    _xmax = xorigin + boxsize
    _ymin = yorigin - buffer_length
    _ymax = yorigin
    x1, y1 = randoms.random_box(size_quarter, _xmin, _xmax, _ymin, _ymax)
    # box edge 2
    _xmin = xorigin + boxsize
    _xmax = xorigin + boxsize + buffer_length
    _ymin = yorigin - buffer_length
    _ymax = yorigin + boxsize
    x2, y2 = randoms.random_box(size_quarter, _xmin, _xmax, _ymin, _ymax)
    # box edge 3
    _xmin = xorigin
    _xmax = xorigin + boxsize + buffer_length
    _ymin = yorigin + boxsize
    _ymax = yorigin + boxsize + buffer_length
    x3, y3 = randoms.random_box(size_quarter, _xmin, _xmax, _ymin, _ymax)
    # box edge 4
    _xmin = xorigin - buffer_length
    _xmax = xorigin
    _ymin = yorigin
    _ymax = yorigin + boxsize + buffer_length
    x4, y4 = randoms.random_box(size_quarter, _xmin, _xmax, _ymin, _ymax)
    x = np.concatenate([x1, x2, x3, x4])
    y = np.concatenate([y1, y2, y3, y4])
    return x, y


def buffer_periodic_2D(data, boxsize, buffer_length, origin=0.):
    """Generates random buffer particles around a 2D box.

    Parameters
    ----------
    data : array
        Where columns 0 and 1 are x and y.
    boxsize : float
        Box size.
    buffer_length : float
        Length of the buffer region.
    origin : float, optional
        Origin point.

    Returns
    -------
    datap : array
        Periodic data values.
    """
    if np.isscalar(origin):
        xorigin, yorigin = origin, origin
    else:
        xorigin, yorigin = origin[0], origin[1]
    assert buffer_length < boxsize, "buffer_length must be smaller than the boxsize."
    ix = np.array([-1, 0, 1])
    ixs, iys = np.meshgrid(ix, ix)
    ixs = ixs.flatten()
    iys = iys.flatten()
    for i in range(0, len(ixs)):
        ix, iy = ixs[i], iys[i]
        if ix != 0 or iy != 0:
            xnew = np.copy(data[:,0]) + ix*boxsize
            ynew = np.copy(data[:,1]) + iy*boxsize
            cond = np.where((xnew >= xorigin-buffer_length) &
                            (xnew <= xorigin+boxsize+buffer_length) &
                            (ynew >= yorigin-buffer_length) &
                            (ynew <= yorigin+boxsize+buffer_length))[0]
            datanew = np.copy(data[cond])
            datanew[:,0] += ix*boxsize
            datanew[:,1] += iy*boxsize
            if ix == -1 and iy == -1:
                datap = datanew
            else:
                datap = np.vstack([datap, datanew])
        else:
            pass
    return datap
