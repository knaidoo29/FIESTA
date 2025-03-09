import numpy as np
from numba import jit, njit

from . import grid

@njit
def bilinear_periodic(fgrid, x, y, xbox, ybox, ngridx, ngridy):
    """
    Bilinear interpolation of field defined on a grid.

    Parameters
    ----------
    fgrid : array
        Field values on the grid.
    x : array
        X coordinates where we need interpolated values.
    y : array
        Y coordinates where we need interpolated values.
    xbox, ybox : float
        Size of the box.
    ngridx, ngridy : int
        Size of the grid alone each axes.
    
    Returns
    -------
    f : array
        Interpolated field values.
    """
    npart = len(x)
    f = np.zeros(npart, dtype=np.float64)
    dx = xbox / float(ngridx)
    dy = ybox / float(ngridy)

    minx = 0.
    
    for i in range(0, npart):
        xp = x[i]
        yp = y[i]
        
        if xp - (dx / 2.0) < 0.0:
            xp += xbox
        if yp - (dy / 2.0) < 0.0:
            yp += ybox
        
        ix1 = int((xp - (dx / 2.0)) / dx)
        xg1 = grid.xgrid(ix1, dx, minx)
        ix2 = ix1 + 1
        xg2 = grid.xgrid(ix2, dx, minx)

        if ix2 == ngridx:
            ix2 = ix2 - ngridx
        
        iy1 = int((yp - (dy / 2.0)) / dy)
        yg1 = grid.xgrid(iy1, dy, minx)
        iy2 = iy1 + 1
        yg2 = grid.xgrid(iy2, dx, minx)

        if iy2 == ngridy:
            iy2 = iy2 - ngridy
        
        q11 = iy1 + ngridy * ix1
        q12 = iy1 + ngridy * ix2
        q21 = iy2 + ngridy * ix1
        q22 = iy2 + ngridy * ix2
        
        f11 = fgrid[q11]
        f12 = fgrid[q12]
        f21 = fgrid[q21]
        f22 = fgrid[q22]
        
        f1 = ((xg2 - xp) / (xg2 - xg1)) * f11 + ((xp - xg1) / (xg2 - xg1)) * f12
        f2 = ((xg2 - xp) / (xg2 - xg1)) * f21 + ((xp - xg1) / (xg2 - xg1)) * f22
        
        f[i] = ((yg2 - yp) / (yg2 - yg1)) * f1 + ((yp - yg1) / (yg2 - yg1)) * f2
    
    return f

@njit
def bilinear_nonperiodic(fgrid, x, y, xbox, ybox, ngridx, ngridy):
    """
    Bilinear interpolation of field defined on a grid.
    
    Parameters
    ----------
    fgrid : array
        Field values on the grid.
    x : array
        X coordinates where we need interpolated values.
    y : array
        Y coordinates where we need interpolated values.
    xbox, ybox : float
        Size of the box.
    ngridx, ngridy : int
        Size of the grid alone each axis.
    
    Returns
    -------
    f : array
        Interpolated field values.
    """
    npart = len(x)
    f = np.zeros(npart, dtype=np.float64)
    dx = xbox / float(ngridx)
    dy = ybox / float(ngridy)

    minx = 0.
    
    for i in range(0, npart):
        xp = x[i]
        yp = y[i]
        
        if xp - dx/2. < 0.:
            ix1 = -1
            ix2 = 0
            xg1 = grid.xgrid(ix1, dx, minx)
            xg2 = grid.xgrid(ix2, dx, minx)
            ix1 = 0
        elif xp > xbox - dx/2.:
            ix1 = ngridx - 1
            ix2 = ngridx
            xg1 = grid.xgrid(ix1, dx, minx)
            xg1 = grid.xgrid(ix2, dx, minx)
            ix2 = ngridx - 1
        else:
            ix1 = int((xp - dx/2.) / dx)
            ix2 = ix1 + 1
            xg1 = grid.xgrid(ix1, dx, minx)
            xg2 = grid.xgrid(ix2, dx, minx)

        if yp - dy/2. < 0.:
            iy1 = -1
            iy2 = 0
            yg1 = grid.xgrid(iy1, dy, minx)
            yg2 = grid.xgrid(iy2, dy, minx)
            iy1 = 0
        elif yp > ybox - dy/2.:
            iy1 = ngridy - 1
            iy2 = ngridy
            yg1 = grid.xgrid(iy1, dy, minx)
            yg2 = grid.xgrid(iy2, dy, minx)
            iy2 = ngridy - 1
        else:
            iy1 = int((yp - dy/2.) / dy)
            iy2 = iy1 + 1
            yg1 = grid.xgrid(iy1, dy, minx)
            yg2 = grid.xgrid(iy2, dy, minx)
        
        q11 = iy1 + ngridy * ix1
        q12 = iy1 + ngridy * ix2
        q21 = iy2 + ngridy * ix1
        q22 = iy2 + ngridy * ix2
        
        f11 = fgrid[q11]
        f12 = fgrid[q12]
        f21 = fgrid[q21]
        f22 = fgrid[q22]
        
        f1 = ((xg2 - xp) / (xg2 - xg1)) * f11 + ((xp - xg1) / (xg2 - xg1)) * f12
        f2 = ((xg2 - xp) / (xg2 - xg1)) * f21 + ((xp - xg1) / (xg2 - xg1)) * f22
        
        f[i] = ((yg2 - yp) / (yg2 - yg1)) * f1 + ((yp - yg1) / (yg2 - yg1)) * f2
    
    return f


@njit
def bilinear_axisperiodic(fgrid, x, y, xbox, ybox, perix, periy, ngridx, ngridy):
    """
    Bilinear interpolation of a field defined on a grid with axis-periodic conditions.
    
    Parameters
    ----------
    fgrid : array
        Field values on the grid (1D flattened array).
    x, y : arrays
        Coordinates where interpolation is needed.
    xbox, ybox : float
        Size of the domain in x and y directions.
    perix, periy : int
        Flags indicating periodicity along each axis (1 = periodic, 0 = non-periodic).
    ngridx, ngridy : int
        Grid resolution along each axis.

    Returns
    -------
    f : array
        Interpolated field values at given (x, y) locations.
    """

    npart = len(x)
    f = np.zeros(npart, dtype=np.float64)
    dx = xbox / float(ngridx)
    dy = ybox / float(ngridy)
    minx = 0.

    for i in range(0, npart):
        xp = x[i]
        yp = y[i]

        # Handle periodicity along x-axis
        if perix == 1:

            if xp - (dx / 2.0) < 0.0:
                xp += xbox
                
            ix1 = int((xp - (dx / 2.0)) / dx)
            xg1 = grid.xgrid(ix1, dx, minx)
            ix2 = ix1 + 1
            xg2 = grid.xgrid(ix2, dx, minx)

            if ix2 == ngridx:
                ix2 = ix2 - ngridx
        else:
            if xp - dx/2. < 0.:
                ix1 = -1
                ix2 = 0
                xg1 = grid.xgrid(ix1, dx, minx)
                xg2 = grid.xgrid(ix2, dx, minx)
                ix1 = 0
            elif xp > xbox - dx/2.:
                ix1 = ngridx - 1
                ix2 = ngridx
                xg1 = grid.xgrid(ix1, dx, minx)
                xg2 = grid.xgrid(ix2, dx, minx)
                ix2 = ngridx - 1
            else:
                ix1 = int((xp - dx/2.) / dx)
                ix2 = ix1 + 1
                xg1 = grid.xgrid(ix1, dx, minx)
                xg2 = grid.xgrid(ix2, dx, minx)

        # Handle periodicity along y-axis
        if periy == 1:
            if yp - (dy / 2.0) < 0.0:
                yp += ybox
                
            iy1 = int((yp - (dy / 2.0)) / dy)
            yg1 = grid.xgrid(iy1, dy, minx)
            iy2 = iy1 + 1
            yg2 = grid.xgrid(iy2, dx, minx)

            if iy2 == ngridy:
                iy2 = iy2 - ngridy
        else:
            if yp - dy/2. < 0.:
                iy1 = -1
                iy2 = 0
                yg1 = grid.xgrid(iy1, dy, minx)
                yg2 = grid.xgrid(iy2, dy, minx)
                iy1 = 0
            elif yp > ybox - dy/2.:
                iy1 = ngridy - 1
                iy2 = ngridy
                yg1 = grid.xgrid(iy1, dy, minx)
                yg2 = grid.xgrid(iy2, dy, minx)
                iy2 = ngridy - 1
            else:
                iy1 = int((yp - dy/2.) / dy)
                iy2 = iy1 + 1
                yg1 = grid.xgrid(iy1, dy, minx)
                yg2 = grid.xgrid(iy2, dy, minx)
        # Flattened 1D indices in the fgrid array
        q11 = iy1 + ngridy * ix1
        q12 = iy1 + ngridy * ix2
        q21 = iy2 + ngridy * ix1
        q22 = iy2 + ngridy * ix2

        f11 = fgrid[q11]
        f12 = fgrid[q12]
        f21 = fgrid[q21]
        f22 = fgrid[q22]

        # Bilinear interpolation
        f1 = ((xg2 - xp) / (xg2 - xg1)) * f11 + ((xp - xg1) / (xg2 - xg1)) * f12
        f2 = ((xg2 - xp) / (xg2 - xg1)) * f21 + ((xp - xg1) / (xg2 - xg1)) * f22

        f[i] = ((yg2 - yp) / (yg2 - yg1)) * f1 + ((yp - yg1) / (yg2 - yg1)) * f2

    return f
