import numpy as np
from numba import njit

from . import grid
from . import part2grid_pix
from . import part2grid_wei


@njit
def part2grid_ngp_3d(x, y, z, f, xlength, ylength, zlength, xmin, ymin, zmin, nxgrid, nygrid, nzgrid):
    """
    Nearest-grid-point assignment in 3D.

    Parameters
    ----------
    x, y, z : array
        Cartesian coordinate system.
    f : array
        Field values are x, y and z coordinates.
    xlength, ylength, zlength : float
        Length of the box along the x, y and z coordinates.
    xmin, ymin, zmin : float
        Minimum values along the x, y and z axis.
    npart : int
        Number of x, y and z coordinates.
    nxgrid, nygrid, nzgrid : int
        Number of grids along x, y and z coordinates.
    
    Returns
    -------
    fgrid : array
        NGP field assignments.
    """
    npart = len(x)
    # Nearest-grid-point assignment in 3D.
    dx = xlength / float(nxgrid)
    dy = ylength / float(nygrid)
    dz = zlength / float(nzgrid)
    wngp = 1. / (dx * dy * dz)

    fgrid = np.zeros(nxgrid * nygrid * nzgrid)

    for i in range(npart):
        xp = x[i]
        yp = y[i]
        zp = z[i]
        fp = f[i]

        xpix = part2grid_pix.ngp_pix(xp, dx, xmin)
        ypix = part2grid_pix.ngp_pix(yp, dy, ymin)
        zpix = part2grid_pix.ngp_pix(zp, dz, zmin)

        if (0 <= xpix < nxgrid) & (0 <= ypix < nygrid) & (0 <= zpix < nzgrid):
            pix = part2grid_pix.pix1dto3d_scalar(xpix, ypix, zpix, nygrid, nzgrid)
            fgrid[pix] += fp * wngp
    return fgrid


@njit
def part2grid_cic_3d(x, y, z, f, xlength, ylength, zlength, xmin, ymin, zmin, nxgrid, nygrid, nzgrid, 
    periodx, periody, periodz):
    """
    Cloud-in-cell assignment in 3D.
    
    Parameters
    ----------
    x, y, z : array
        Cartesian coordinate system.
    f : array
        Field values are x, y & z coordinates.
    xlength, ylength, zlength : float
        Length of the box along the x, y & z coordinates.
    xmin, ymin : float
        Minimum values along the x, y & z axis.
    nxgrid, nygrid : int
        Number of grids along x, y & z coordinates.
    periodx, periody, periodz : bool
        Periodic boundary conditions.
    
    Returns
    -------
    fgrid : array
        CIC field assignments.
    """
    npart = len(x)
    # Cloud-in-cell assignment in 3D.
    dx = xlength / float(nxgrid)
    dy = ylength / float(nygrid)
    dz = zlength / float(nzgrid)

    fgrid = np.zeros(nxgrid * nygrid * nzgrid)

    for i in range(npart):
        xp = x[i]
        yp = y[i]
        zp = z[i]
        fp = f[i]

        xpix = part2grid_pix.cic_pix(xp, dx, xmin)
        ypix = part2grid_pix.cic_pix(yp, dy, ymin)
        zpix = part2grid_pix.cic_pix(zp, dz, zmin)

        xg = grid.xgrids(xpix, dx, xmin)
        yg = grid.xgrids(ypix, dy, ymin)
        zg = grid.xgrids(zpix, dz, zmin)

        if periodx:
            xpix = part2grid_pix.periodic_pix(xpix, nxgrid)
        if periody:
            ypix = part2grid_pix.periodic_pix(ypix, nygrid)
        if periodz:
            zpix = part2grid_pix.periodic_pix(zpix, nzgrid)

        for j1 in range(2):
            for j2 in range(2):
                for j3 in range(2):
                    if (xpix[j1] >= 0) and (xpix[j1] < nxgrid) and (ypix[j2] >= 0) and (ypix[j2] < nygrid) and (zpix[j3] >= 0) and (zpix[j3] < nygrid):
                        pix = part2grid_pix.pix1dto3d_scalar(xpix[j1], ypix[j2], zpix[j3], nygrid, nzgrid)
                        wx = part2grid_wei.weight_cic(xp, xg[j1], dx)
                        wy = part2grid_wei.weight_cic(yp, yg[j2], dy)
                        wz = part2grid_wei.weight_cic(zp, zg[j3], dz)
                        fgrid[pix] += fp * wx * wy * wz
    return fgrid


@njit
def part2grid_tsc_3d(x, y, z, f, xlength, ylength, zlength, xmin, ymin, zmin, nxgrid, nygrid, nzgrid, 
    periodx, periody, periodz):
    """
    Triangular-shaped-cloud assignment in 3D.

    Parameters
    ----------
    x, y, z : array
        Cartesian coordinate system.
    f : array
        Field values are x, y & z coordinates.
    xlength, ylength, zlength : float
        Length of the box along the x, y & z coordinates.
    xmin, ymin, zmin : float
        Minimum values along the x, y & z axis.
    nxgrid, nygrid, nzgrid : int
        Number of grids along x, y & z coordinates.
    periodx, periody, periodz : bool
        Periodic boundary conditions.
    
    Returns
    -------
    fgrid : array
        TSC field assignments.
    """
    npart = len(x)
    # Triangular-shaped-cloud assignment in 3D.
    dx = xlength / float(nxgrid)
    dy = ylength / float(nygrid)
    dz = zlength / float(nzgrid)

    fgrid = np.zeros(nxgrid * nygrid * nzgrid)

    for i in range(npart):
        xp = x[i]
        yp = y[i]
        zp = z[i]
        fp = f[i]

        xpix = part2grid_pix.tsc_pix(xp, dx, xmin)
        ypix = part2grid_pix.tsc_pix(yp, dy, ymin)
        zpix = part2grid_pix.tsc_pix(zp, dz, zmin)

        xg = grid.xgrids(xpix, dx, xmin)
        yg = grid.xgrids(ypix, dy, ymin)
        zg = grid.xgrids(zpix, dz, zmin)

        if periodx:
            xpix = part2grid_pix.periodic_pix(xpix, nxgrid)
        if periody:
            ypix = part2grid_pix.periodic_pix(ypix, nygrid)
        if periodz:
            zpix = part2grid_pix.periodic_pix(zpix, nzgrid)

        for j1 in range(3):
            for j2 in range(3):
                for j3 in range(3):
                    if (xpix[j1] >= 0) and (xpix[j1] < nxgrid) and (ypix[j2] >= 0) and (ypix[j2] < nygrid) and (zpix[j3] >= 0) and (zpix[j3] < nygrid):
                        pix = part2grid_pix.pix1dto3d_scalar(xpix[j1], ypix[j2], zpix[j3], nygrid, nzgrid)
                        wx = part2grid_wei.weight_tsc(xp, xg[j1], dx)
                        wy = part2grid_wei.weight_tsc(yp, yg[j2], dy)
                        wz = part2grid_wei.weight_tsc(zp, zg[j3], dz)
                        fgrid[pix] += fp * wx * wy * wz
    return fgrid


@njit
def part2grid_pcs_3d(x, y, z, f, xlength, ylength, zlength, xmin, ymin, zmin, nxgrid, nygrid, nzgrid,
    periodx, periody, periodz):
    """
    Piecewise-Cubic-Spline assignment in 3D.
    
    Parameters
    ----------
    x, y, z : array
        Cartesian coordinate system.
    f : array
        Field values are x, y & z coordinates.
    xlength, ylength, zlength : float
        Length of the box along the x, y & z coordinates.
    xmin, ymin, zmin : float
        Minimum values along the x, y & z axis.
    nxgrid, nygrid, nzgrid : int
        Number of grids along x, y & z coordinates.
    periodx, periody, periodz : bool
        Periodic boundary conditions.
    
    Returns
    -------
    fgrid : array
        PCS field assignments.
    """
    npart = len(x)
    # Piecewise-Cubic-Spline assignment in 3D.
    dx = xlength / float(nxgrid)
    dy = ylength / float(nygrid)
    dz = zlength / float(nzgrid)

    fgrid = np.zeros(nxgrid * nygrid * nzgrid)

    for i in range(npart):
        xp = x[i]
        yp = y[i]
        zp = z[i]
        fp = f[i]

        xpix = part2grid_pix.pcs_pix(xp, dx, xmin)
        ypix = part2grid_pix.pcs_pix(yp, dy, ymin)
        zpix = part2grid_pix.pcs_pix(zp, dz, zmin)

        xg = grid.xgrids(xpix, dx, xmin)
        yg = grid.xgrids(ypix, dy, ymin)
        zg = grid.xgrids(zpix, dz, zmin)

        if periodx:
            xpix = part2grid_pix.periodic_pix(xpix, nxgrid)
        if periody:
            ypix = part2grid_pix.periodic_pix(ypix, nygrid)
        if periodz:
            zpix = part2grid_pix.periodic_pix(zpix, nzgrid)

        for j1 in range(4):
            for j2 in range(4):
                for j3 in range(4):
                    if (xpix[j1] >= 0) and (xpix[j1] < nxgrid) and (ypix[j2] >= 0) and (ypix[j2] < nygrid) and (zpix[j3] >= 0) and (zpix[j3] < nygrid):
                        pix = part2grid_pix.pix1dto3d_scalar(xpix[j1], ypix[j2], zpix[j3], nygrid, nzgrid)
                        wx = part2grid_wei.weight_pcs(xp, xg[j1], dx)
                        wy = part2grid_wei.weight_pcs(yp, yg[j2], dy)
                        wz = part2grid_wei.weight_pcs(zp, zg[j3], dz)
                        fgrid[pix] += fp * wx * wy * wz
    return fgrid
