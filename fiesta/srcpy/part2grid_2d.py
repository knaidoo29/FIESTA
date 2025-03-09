import numpy as np
from numba import njit

from . import grid
from . import part2grid_pix
from . import part2grid_wei


@njit
def part2grid_ngp_2d(x, y, f, xlength, ylength, xmin, ymin, nxgrid, nygrid):
    """
    Nearest-grid-point assignment in 2D.
    
    Parameters
    ----------
    x, y : array
        Cartesian coordinate system.
    f : array
        Field values are x & y coordinates.
    xlength, ylength : float
        Length of the box along the x and y coordinates.
    xmin, ymin : float
        Minimum values along the x and y axis.
    nxgrid, nygrid : int
        Number of grids along x and y coordinates.
    
    Returns
    -------
    fgrid : array
        NGP field assignments.
    """
    npart = len(x)
    dx = xlength / nxgrid
    dy = ylength / nygrid
    wngp = 1. / (dx * dy)

    # Initialize fgrid to zero
    fgrid = np.zeros(nxgrid*nygrid)

    for i in range(0, npart):
        xp, yp, fp = x[i], y[i], f[i]
        
        xpix = part2grid_pix.ngp_pix(xp, dx, xmin)
        ypix = part2grid_pix.ngp_pix(yp, dy, ymin)

        if (xpix >= 0) and (xpix < nxgrid) and (ypix >= 0) and (ypix < nygrid):
            pix = part2grid_pix.pix1dto2d_scalar(xpix, ypix, nygrid)
            fgrid[pix] += fp * wngp
    return fgrid


@njit
def part2grid_cic_2d(x, y, f, xlength, ylength, xmin, ymin, nxgrid, nygrid, periodx, periody):
    """
    Cloud-in-cell assignment in 2D.

    Parameters
    ----------
    x, y : array
        Cartesian coordinate system.
    f : array
        Field values are x & y coordinates.
    xlength, ylength : float
        Length of the box along the x and y coordinates.
    xmin, ymin : float
        Minimum values along the x and y axis.
    npart : int
        Number of x and y coordinates.
    nxgrid, nygrid : int
        Number of grids along x and y coordinates.
    periodx, periody, periodz : bool
        Periodic boundary conditions.
    
    Returns
    -------
    fgrid : array
        CIC field assignments.
    """
    npart = len(x)
    dx = xlength / nxgrid
    dy = ylength / nygrid

    # Initialize fgrid to zero
    fgrid = np.zeros(nxgrid*nygrid)

    for i in range(npart):
        xp, yp, fp = x[i], y[i], f[i]

        # Call cic_pix to get xpix and ypix (not implemented here, assume it exists)
        xpix = part2grid_pix.cic_pix(xp, dx, xmin)
        ypix = part2grid_pix.cic_pix(yp, dy, ymin)
        
        # Call xgrids to compute the grid positions (not implemented here, assume it exists)
        xg = grid.xgrids(xpix, dx, xmin)
        yg = grid.xgrids(ypix, dy, ymin)

        if periodx:
            xpix = part2grid_pix.periodic_pix(xpix, nxgrid)
        if periody:
            ypix = part2grid_pix.periodic_pix(ypix, nygrid)

        for j1 in range(2):
            for j2 in range(2):
                if (xpix[j1] >= 0) and (xpix[j1] < nxgrid) and (ypix[j2] >= 0) and (ypix[j2] < nygrid):
                    pix = part2grid_pix.pix1dto2d_scalar(xpix[j1], ypix[j2], nygrid)
                    wx = part2grid_wei.weight_cic(xp, xg[j1], dx)
                    wy = part2grid_wei.weight_cic(yp, yg[j2], dy)
                    fgrid[pix] += fp * wx * wy
    return fgrid


@njit
def part2grid_tsc_2d(x, y, f, xlength, ylength, xmin, ymin, nxgrid, nygrid, periodx, periody):
    """
    Triangular-shaped-cloud assignment in 2D.
    
    Parameters
    ----------
    x, y : array
        Cartesian coordinate system.
    f : array
        Field values are x & y coordinates.
    xlength, ylength : float
        Length of the box along the x and y coordinates.
    xmin, ymin : float
        Minimum values along the x and y axis.
    npart : int
        Number of x and y coordinates.
    nxgrid, nygrid : int
        Number of grids along x and y coordinates.
    periodx, periody : bool
        Periodic boundary conditions.
    
    Returns
    -------
    fgrid : array
        TSC field assignments.
    """
    npart = len(x)
    dx = xlength / nxgrid
    dy = ylength / nygrid

    # Initialize fgrid to zero
    fgrid = np.zeros(nxgrid*nygrid)

    for i in range(npart):
        xp, yp, fp = x[i], y[i], f[i]

        # Call tsc_pix to get xpix and ypix (not implemented here, assume it exists)
        xpix = part2grid_pix.tsc_pix(xp, dx, xmin)
        ypix = part2grid_pix.tsc_pix(yp, dy, ymin)
        
        # Call xgrids to compute the grid positions (not implemented here, assume it exists)
        xg = grid.xgrids(xpix, dx, xmin)
        yg = grid.xgrids(ypix, dy, ymin)

        if periodx:
            xpix = part2grid_pix.periodic_pix(xpix, nxgrid)
        if periody:
            ypix = part2grid_pix.periodic_pix(ypix, nygrid)

        for j1 in range(3):
            for j2 in range(3):
                if (xpix[j1] >= 0) and (xpix[j1] < nxgrid) and (ypix[j2] >= 0) and (ypix[j2] < nygrid):
                    pix = part2grid_pix.pix1dto2d_scalar(xpix[j1], ypix[j2], nygrid)
                    wx = part2grid_wei.weight_tsc(xp, xg[j1], dx)
                    wy = part2grid_wei.weight_tsc(yp, yg[j2], dy)
                    fgrid[pix] += fp * wx * wy
    return fgrid


@njit
def part2grid_pcs_2d(x, y, f, xlength, ylength, xmin, ymin, nxgrid, nygrid, periodx, periody):
    """
    Piecewise-Cubic-Spline assignment in 2D.

    Parameters
    ----------
    x, y : array
        Cartesian coordinate system.
    f : array
        Field values are x & y coordinates.
    xlength, ylength : float
        Length of the box along the x and y coordinates.
    xmin, ymin : float
        Minimum values along the x and y axis.
    npart : int
        Number of x and y coordinates.
    nxgrid, nygrid : int
        Number of grids along x and y coordinates.
    periodx, periody : bool
        Periodic boundary conditions.
    
    Returns
    -------
    fgrid : array
        TSC field assignments.
    """
    npart = len(x)
    dx = xlength / nxgrid
    dy = ylength / nygrid

    # Initialize fgrid to zero
    fgrid = np.zeros(nxgrid*nygrid)

    for i in range(npart):
        xp, yp, fp = x[i], y[i], f[i]

        # Call pcs_pix to get xpix and ypix (not implemented here, assume it exists)
        xpix = part2grid_pix.pcs_pix(xp, dx, xmin) 
        ypix = part2grid_pix.pcs_pix(yp, dy, ymin)
        
        # Call xgrids to compute the grid positions (not implemented here, assume it exists)
        xg = grid.xgrids(xpix, dx, xmin)
        yg = grid.xgrids(ypix, dy, ymin)

        if periodx:
            xpix = part2grid_pix.periodic_pix(xpix, nxgrid)
        if periody:
            ypix = part2grid_pix.periodic_pix(ypix, nygrid)

        for j1 in range(4):
            for j2 in range(4):
                if (xpix[j1] >= 0) and (xpix[j1] < nxgrid) and (ypix[j2] >= 0) and (ypix[j2] < nygrid):
                    pix = part2grid_pix.pix1dto2d_scalar(xpix[j1], ypix[j2], nygrid)
                    wx = part2grid_wei.weight_pcs(xp, xg[j1], dx)
                    wy = part2grid_wei.weight_pcs(yp, yg[j2], dy)
                    fgrid[pix] += fp * wx * wy
    return fgrid
