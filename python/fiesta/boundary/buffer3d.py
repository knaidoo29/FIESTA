import numpy as np

from .. import randoms


def buffer_random_particles_3d(npart, boxsize, buffer_length):
    """Generates random buffer particles around a 3D box.

    Parameters
    ----------
    npart : int
        Number of particles in the 3D box.
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
    z : array
        Random z-values.
    """

    box_vol = boxsize**3.
    part_dens = npart / box_vol
    bufferbox_vol = (boxsize + 2.*buffer_length)**3.
    buffer_vol = bufferbox_vol - box_vol

    # Divide into three regions along z-axis
    # Region 1: zrange = [-buffer_length, 0.]
    # Region 2: zrange = [0., boxsize]
    # Region 3: zrange = [boxsize, boxsize+buffer_length]

    vol1 = (boxsize + 2.*buffer_length)*buffer_length
    vol3 = vol1
    vol2 = buffer_vol - vol1 - vol3

    npart_buffer1 = int(part_dens * vol1)
    npart_buffer2 = int(part_dens * vol2)
    npart_buffer3 = int(part_dens * vol3)

    # Splitting region 2, into four regions.
    # Subregion 21: xrange = [-buffer_length, boxsize], yrange = [-buffer_length, 0]
    # Subregion 22: xrange = [boxsize, boxsize+buffer_length], yrange =[-buffer_length, boxsize]
    # Subregion 23: xrange = [0, boxsize+buffer_length], yrange=[boxsize, boxsize+buffer_length]
    # Subregion 24: xrange = [-buffer_length, 0], yrange = [0, boxsize+buffer_length]

    subvol = vol2/4.
    npart_subbuffer2 = int(part_dens * subvol)

    # Region 1

    xmin, xmax = -buffer_length, boxsize+buffer_length
    ymin, ymax = -buffer_length, boxsize+buffer_length
    zmin, zmax = -buffer_length, 0.
    x1, y1, z1 = randoms.random_cube(npart_buffer1, xmin, xmax, ymin, ymax, zmin, zmax)

    # Region 3

    xmin, xmax = -buffer_length, boxsize+buffer_length
    ymin, ymax = -buffer_length, boxsize+buffer_length
    zmin, zmax = boxsize, boxsize+buffer_length
    x3, y3, z3 = randoms.random_cube(npart_buffer3, xmin, xmax, ymin, ymax, zmin, zmax)

    # Subregion 21

    xmin, xmax = -buffer_length, boxsize
    ymin, ymax = -buffer_length, 0.
    zmin, zmax = 0., boxsize
    x21, y21, z21 = randoms.random_cube(npart_subbuffer2, xmin, xmax, ymin, ymax, zmin, zmax)

    # Subregion 22

    xmin, xmax = boxsize, boxsize+buffer_length
    ymin, ymax = -buffer_length, boxsize
    zmin, zmax = 0., boxsize
    x22, y22, z22 = randoms.random_cube(npart_subbuffer2, xmin, xmax, ymin, ymax, zmin, zmax)

    # Subregion 23

    xmin, xmax = 0., boxsize+buffer_length
    ymin, ymax = boxsize, boxsize+buffer_length
    zmin, zmax = 0., boxsize
    x23, y23, z23 = randoms.random_cube(npart_subbuffer2, xmin, xmax, ymin, ymax, zmin, zmax)

    # Subregion 24

    xmin, xmax = -buffer_length, 0.
    ymin, ymax = 0., boxsize+buffer_length
    zmin, zmax = 0., boxsize
    x24, y24, z24 = randoms.random_cube(npart_subbuffer2, xmin, xmax, ymin, ymax, zmin, zmax)

    x = np.concatenate([x1, x21, x22, x23, x24, x3])
    y = np.concatenate([y1, y21, y22, y23, y24, y3])
    z = np.concatenate([z1, z21, z22, z23, z24, z3])

    return x, y, z


def buffer_periodic_particles_3d(x, y, z, boxsize, buffer_length):
    """Generates random buffer particles around a 3D box.

    Parameters
    ----------
    x : array
        X-coordinates.
    y : array
        Y-coordinates.
    z : array
        Z-coordinates.
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
    zp : array
        Periodic z-values.
    """
    assert buffer_length < boxsize, "buffer_length must be smaller than the boxsize."
    ix = np.array([-1, 0, 1])
    ixs, iys, izs = np.meshgrid(ix, ix, ix)
    ixs = ixs.flatten()
    iys = iys.flatten()
    izs = izs.flatten()
    for i in range(0, len(ixs)):
        ix, iy, iz = ixs[i], iys[i], izs[i]
        if ix != 0 or iy != 0 or iz != 0:
            xnew = np.copy(x) + ix*boxsize
            ynew = np.copy(y) + iy*boxsize
            znew = np.copy(z) + iz*boxsize
            cond = np.where((xnew >= -buffer_length) & (xnew <= boxsize+buffer_length) &
                            (ynew >= -buffer_length) & (ynew <= boxsize+buffer_length) &
                            (znew >= -buffer_length) & (znew <= boxsize+buffer_length))[0]
            if ix == -1 and iy == -1 and iz == -1:
                xp = xnew[cond]
                yp = ynew[cond]
                zp = znew[cond]
            else:
                xp = np.concatenate([xp, xnew[cond]])
                yp = np.concatenate([yp, ynew[cond]])
                zp = np.concatenate([zp, znew[cond]])
        else:
            pass
    return xp, yp, zp
