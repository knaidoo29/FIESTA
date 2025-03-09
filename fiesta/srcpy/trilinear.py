import numpy as np
from numba import njit

from . import grid


@njit
def trilinear_periodic(fgrid, x, y, z, xbox, ybox, zbox, ngridx, ngridy, ngridz):
    """
    Trilinear interpolation of field defined on a grid.

    Parameters
    ----------
    fgrid : array
        Field values on the grid.
    x : array
        X coordinates where we need interpolated values.
    y : array
        Y coordinates where we need interpolated values.
    z : array
        Z coordinates where we need interpolated values.
    xbox, ybox, zbox : float
        Size of the box.
    ngridx, ngridy, ngridz : int
        Size of the grid along each axis.
    npart : int
        Number of particles.

    Returns
    -------
    f : array
        Interpolated field values.
    """
    npart = len(x)
    f = np.zeros(npart)
    minx = 0.
    dx = xbox / float(ngridx)
    dy = ybox / float(ngridy)
    dz = zbox / float(ngridz)
    for i in range(0, npart):
        xp, yp, zp = x[i], y[i], z[i]

        if xp - dx/2. < 0.:
            xp = xp + xbox
        if yp - dy/2. < 0.:
            yp = yp + ybox
        if zp - dz/2. < 0.:
            zp = zp + zbox

        ix1 = int((xp - dx/2.) / dx)
        xg1 = grid.xgrid(ix1, dx, minx)

        ix2 = ix1 + 1
        xg2 = grid.xgrid(ix2, dx, minx)

        if ix2 == ngridx:
            ix2 = ix2 - ngridx

        iy1 = int((yp - dy/2.) / dy)
        yg1 = grid.xgrid(iy1, dy, minx)

        iy2 = iy1 + 1
        yg2 = grid.xgrid(iy2, dy, minx)

        if iy2 == ngridy:
            iy2 = iy2 - ngridy

        iz1 = int((zp - dz/2.) / dz)
        zg1 = grid.xgrid(iz1, dz, minx)

        iz2 = iz1 + 1
        zg2 = grid.xgrid(iz2, dz, minx)

        if iz2 == ngridz:
            iz2 = iz2 - ngridz

        # surround points in the grid of a single point for interpolation.

        q111 = iz1 + ngridz*(iy1 + ngridy*ix1)
        q112 = iz1 + ngridz*(iy1 + ngridy*ix2)
        q121 = iz1 + ngridz*(iy2 + ngridy*ix1) 
        q122 = iz1 + ngridz*(iy2 + ngridy*ix2) 
        q211 = iz2 + ngridz*(iy1 + ngridy*ix1) 
        q212 = iz2 + ngridz*(iy1 + ngridy*ix2) 
        q221 = iz2 + ngridz*(iy2 + ngridy*ix1)
        q222 = iz2 + ngridz*(iy2 + ngridy*ix2)

        f111 = fgrid[q111]
        f112 = fgrid[q112]
        f121 = fgrid[q121]
        f122 = fgrid[q122]
        f211 = fgrid[q211]
        f212 = fgrid[q212]
        f221 = fgrid[q221]
        f222 = fgrid[q222]

        xd = (xp - xg1) / (xg2 - xg1)
        yd = (yp - yg1) / (yg2 - yg1)
        zd = (zp - zg1) / (zg2 - zg1)

        f11 = f111*(1-xd) + f112*xd
        f21 = f211*(1-xd) + f212*xd
        f12 = f121*(1-xd) + f122*xd
        f22 = f221*(1-xd) + f222*xd

        f1 = f11*(1-yd) + f12*yd
        f2 = f21*(1-yd) + f22*yd

        f[i] = f1*(1-zd) + f2*zd

    return f


@njit
def trilinear_nonperiodic(fgrid, x, y, z, xbox, ybox, zbox, ngridx, ngridy, ngridz):
    """
    Trilinear interpolation of field defined on a grid.
  
    Parameters
    ----------
    fgrid : array
        Field values on the grid.
    x : array
        X coordinates where we need interpolated values.
    y : array
        Y coordinates where we need interpolated values.
    z : array
        Z coordinates where we need interpolated values.
    xbox, ybox, zbox : float
        Size of the box.
    ngridx, ngridy, ngridz : int
        Size of the grid along each axis.
    
    Returns
    -------
    f : array
        Interpolated field values.
    """
    npart = len(x)
    f = np.zeros(npart)
    minx = 0.
    dx = xbox / float(ngridx)
    dy = ybox / float(ngridy)
    dz = zbox / float(ngridz)
    for i in range(0, npart):
        xp, yp, zp = x[i], y[i], z[i]

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
            xg1 = grid.xgrid(ix1, dy, minx)
            xg2 = grid.xgrid(ix2, dy, minx)

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

        if zp - dz/2. < 0.:
            iz1 = -1
            iz2 = 0
            zg1 = grid.xgrid(iz1, dz, minx)
            zg2 = grid.xgrid(iz2, dz, minx)
            iz1 = 0
        elif zp > zbox - dz/2.:
            iz1 = ngridz - 1
            iz2 = ngridz
            zg1 = grid.xgrid(iz1, dz, minx)
            zg2 = grid.xgrid(iz2, dz, minx)
            iz2 = ngridz - 1
        else:
            iz1 = int((zp - dz/2.) / dz)
            iz2 = iz1 + 1
            zg1 = grid.xgrid(iz1, dz, minx)
            zg2 = grid.xgrid(iz2, dz, minx)

        # surround points in the grid of a single point for interpolation.

        q111 = iz1 + ngridz*(iy1 + ngridy*ix1)
        q112 = iz1 + ngridz*(iy1 + ngridy*ix2)        
        q121 = iz1 + ngridz*(iy2 + ngridy*ix1) 
        q122 = iz1 + ngridz*(iy2 + ngridy*ix2) 
        q211 = iz2 + ngridz*(iy1 + ngridy*ix1) 
        q212 = iz2 + ngridz*(iy1 + ngridy*ix2) 
        q221 = iz2 + ngridz*(iy2 + ngridy*ix1) 
        q222 = iz2 + ngridz*(iy2 + ngridy*ix2)

        f111 = fgrid[q111]
        f112 = fgrid[q112]
        f121 = fgrid[q121]
        f122 = fgrid[q122]
        f211 = fgrid[q211]
        f212 = fgrid[q212]
        f221 = fgrid[q221]
        f222 = fgrid[q222]

        xd = (xp - xg1) / (xg2 - xg1)
        yd = (yp - yg1) / (yg2 - yg1)
        zd = (zp - zg1) / (zg2 - zg1)

        f11 = f111*(1-xd) + f112*xd
        f21 = f211*(1-xd) + f212*xd
        f12 = f121*(1-xd) + f122*xd
        f22 = f221*(1-xd) + f222*xd

        f1 = f11*(1-yd) + f12*yd
        f2 = f21*(1-yd) + f22*yd

        f[i] = f1*(1-zd) + f2*zd

    return f


@njit
def trilinear_axisperiodic(fgrid, x, y, z, xbox, ybox, zbox, perix, periy, periz, ngridx, ngridy, ngridz):
    """
    Trilinear interpolation of field defined on a grid.

    Parameters
    ----------
    fgrid : array
        Field values on the grid.
    x : array
        X coordinates where we need interpolated values.
    y : array
        Y coordinates where we need interpolated values.
    z : array
        Z coordinates where we need interpolated values.
    xbox, ybox, zbox : float
        Size of the box.
    perix, periy, periz : int
        0 = non-periodic, 1 = periodic
    ngridx, ngridy, ngridz : int
        Size of the grid along each axis.
    npart : int
        Number of particles.

    Returns
    -------
    f : array
        Interpolated field values.
    """
    npart = len(x)
    f = np.zeros(npart)
    minx = 0.
    dx = xbox / float(ngridx)
    dy = ybox / float(ngridy)
    dz = zbox / float(ngridz)
    for i in range(0, npart):
        xp, yp, zp = x[i], y[i], z[i]

        if perix == 1:
            if xp - dx/2. < 0.:
                xp = xp + xbox
            ix1 = int((xp - dx/2.) / dx)
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
                xg1 = grid.xgrid(ix1, dy, minx)
                xg2 = grid.xgrid(ix2, dy, minx)

        if periy == 1:
            if yp - dy/2. < 0.:
                yp = yp + ybox
            iy1 = int((yp - dy/2.) / dy)
            yg1 = grid.xgrid(iy1, dy, minx)
            iy2 = iy1 + 1
            yg2 = grid.xgrid(iy2, dy, minx)
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

        if periz == 1:
            if zp - dz/2. < 0.:
                zp = zp + zbox
            iz1 = int((zp - dz/2.) / dz)
            zg1 = grid.xgrid(iz1, dz, minx)
            iz2 = iz1 + 1
            zg2 = grid.xgrid(iz2, dz, minx)
            if iz2 == ngridz:
                iz2 = iz2 - ngridz
        else:
            if zp - dz/2. < 0.:
                iz1 = -1
                iz2 = 0
                zg1 = grid.xgrid(iz1, dz, minx)
                zg2 = grid.xgrid(iz2, dz, minx)
                iz1 = 0
            elif zp > zbox - dz/2.:
                iz1 = ngridz - 1
                iz2 = ngridz
                zg1 = grid.xgrid(iz1, dz, minx)
                zg2 = grid.xgrid(iz2, dz, minx)
                iz2 = ngridz - 1
            else:
                iz1 = int((zp - dz/2.) / dz)
                iz2 = iz1 + 1
                zg1 = grid.xgrid(iz1, dz, minx)
                zg2 = grid.xgrid(iz2, dz, minx)

        # surround points in the grid of a single point for interpolation.

        q111 = iz1 + ngridz*(iy1 + ngridy*ix1)
        q112 = iz1 + ngridz*(iy1 + ngridy*ix2)
        q121 = iz1 + ngridz*(iy2 + ngridy*ix1)
        q122 = iz1 + ngridz*(iy2 + ngridy*ix2)
        q211 = iz2 + ngridz*(iy1 + ngridy*ix1)
        q212 = iz2 + ngridz*(iy1 + ngridy*ix2)
        q221 = iz2 + ngridz*(iy2 + ngridy*ix1)
        q222 = iz2 + ngridz*(iy2 + ngridy*ix2)

        f111 = fgrid[q111]
        f112 = fgrid[q112]
        f121 = fgrid[q121]
        f122 = fgrid[q122]
        f211 = fgrid[q211]
        f212 = fgrid[q212]
        f221 = fgrid[q221]
        f222 = fgrid[q222]

        xd = (xp - xg1) / (xg2 - xg1)
        yd = (yp - yg1) / (yg2 - yg1)
        zd = (zp - zg1) / (zg2 - zg1)

        f11 = f111*(1-xd) + f112*xd
        f21 = f211*(1-xd) + f212*xd
        f12 = f121*(1-xd) + f122*xd
        f22 = f221*(1-xd) + f222*xd

        f1 = f11*(1-yd) + f12*yd
        f2 = f21*(1-yd) + f22*yd

        f[i] = f1*(1-zd) + f2*zd

    return f