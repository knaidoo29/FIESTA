import numpy as np
from numba import njit

from . import integral_image
from .. import p2g


@njit
def sum_from_integral_image_2D(igrid, ixmin, ixmax, iymin, iymax, periodic=True):
    """
    Quick summations from an integral image in 2D.

    Parameters
    ----------
    igrid : array_like
        Integral image in 2D.
    ixmin : int
        Minimum grid value along the x grid.
    ixmax : int
        Maximum grid value along the x grid.
    iymin : int
        Minimum grid value along the y grid.
    iymax : int
        Maximum grid value along the y grid.
    periodic : bool, optional
        Periodic boundary condition.
    
    Yields
    ------
    The sum of the integral image within a rectangle.
    """
    if ixmin >= 0 and ixmax < len(igrid[0]) and iymin >= 0 and iymax < len(igrid[:,0]):
        return igrid[ixmax,iymax] - igrid[ixmin,iymax] - igrid[ixmax,iymin] + igrid[ixmin,iymin]
    else:
        if periodic == False:
            if ixmin < 0:
                ixmin = 0
            if ixmax >= len(igrid[0]):
                ixmax = len(igrid[0])-1
            if iymin < 0:
                iymin = 0
            if iymax >= len(igrid[:,0]):
                iymax = len(igrid[:,0])-1
            return igrid[ixmax,iymax] - igrid[ixmin,iymax] - igrid[ixmax,iymin] + igrid[ixmin,iymin]
        else:
            if ixmin < 0:
                ixmin1 = ixmin+len(igrid[0])
                ixmax1 = len(igrid[0])-1
                ixmin2 = 0
                ixmax2 = ixmax
            elif ixmax > len(igrid[0])-1:
                ixmin1 = ixmin
                ixmax1 = len(igrid[0])-1
                ixmin2 = 0
                ixmax2 = ixmax-len(igrid[0])
            else:
                ixmin1 = -2*len(igrid[0])
                ixmax1 = -2*len(igrid[0])
                ixmin2 = -2*len(igrid[0])
                ixmax2 = -2*len(igrid[0])
            if iymin < 0:
                iymin1 = iymin+len(igrid[:,0])
                iymax1 = len(igrid[:,0])-1
                iymin2 = 0
                iymax2 = iymax
            elif iymax > len(igrid[:,0])-1:
                iymin1 = iymin
                iymax1 = len(igrid[:,0])-1
                iymin2 = 0
                iymax2 = iymax-len(igrid[:,0])
            else:
                iymin1 = -2*len(igrid[:,0])
                iymax1 = -2*len(igrid[:,0])
                iymin2 = -2*len(igrid[:,0])
                iymax2 = -2*len(igrid[:,0])
            if ixmin1 is -2*len(igrid[0]):
                ixmins = [ixmin, ixmin]
                ixmaxs = [ixmax, ixmax]
                iymins = [iymin1, iymin2]
                iymaxs = [iymax1, iymax2]
            elif iymin1 is -2*len(igrid[:,0]):
                ixmins = [ixmin1, ixmin2]
                ixmaxs = [ixmax1, ixmax2]
                iymins = [iymin, iymin]
                iymaxs = [iymax, iymax]
            else:
                ixmins = [ixmin1, ixmin1, ixmin2, ixmin2]
                ixmaxs = [ixmax1, ixmax1, ixmax2, ixmax2]
                iymins = [iymin1, iymin2, iymin1, iymin2]
                iymaxs = [iymax1, iymax2, iymax1, iymax2]
            isum = 0
            for i in range(0, len(ixmins)):
                isum += igrid[ixmaxs[i],iymaxs[i]] - igrid[ixmins[i],iymaxs[i]] - igrid[ixmaxs[i],iymins[i]] + igrid[ixmins[i],iymins[i]]
            return isum


@njit
def get_volume_enclosing_box(boxsize, ngrid, dgrid, idgrid, minpart, periodic=True, ifgrid=None):
    """
    Get the volume for an enclosing box containing minpart number of particles.

    Parameters
    ----------
    boxsize : float
        Size of the 2D grid.
    ngrid : int
       Grid size.
    dgrid : array_like
        The 2D density field grid.
    idgrid : array_like
        The 2D density field integral image.
    minpart : int
        Minimum number of particles.
    periodic : bool, optional
        Periodic boundary condition. 
    ifgrid : array_like
        Field integral image. If not None, then the field will be estimate via 
        the volume enclosing box method.
    
    Returns
    -------
    vgrid : array_like
        The 2D volume enclosing box for density estimation.
    """
    ngrid = len(dgrid)
    dx = boxsize / ngrid
    vgrid = np.zeros(np.shape(dgrid), dtype=np.float64)
    if ifgrid is not None:
        fgridVEB = np.zeros(np.shape(dgrid), dtype=np.float64)
    for i in range(0, np.shape(dgrid)[0]):
        for j in range(0, np.shape(dgrid)[1]):
            counts = dgrid[i,j]*dx*dx
            if counts >= minpart:
                vol = minpart/(counts/(dx**2))
                if ifgrid is not None:
                    iadd = 1
                    isub = 0
                    fgridVEB[i,j] = sum_from_integral_image_2D(ifgrid, i+isub, i+iadd, j+isub, j+iadd, periodic=periodic)/dgrid[i,j]
            else:
                iadd = 1
                isub = 0
                # smaller enclosing volume
                IcountS = sum_from_integral_image_2D(idgrid, i+isub, i+iadd, j+isub, j+iadd, periodic=periodic)
                IcountS *= dx*dx
                if ifgrid is not None:
                    fsumS = sum_from_integral_image_2D(ifgrid, i+isub, i+iadd, j+isub, j+iadd, periodic=periodic)
                    fsumS *= dx*dx
                # larger enclosing volume
                iadd += 1
                isub -= 1
                IcountL = sum_from_integral_image_2D(idgrid, i+isub, i+iadd, j+isub, j+iadd, periodic=periodic)
                IcountL *= dx*dx
                if ifgrid is not None:
                    fsumL = sum_from_integral_image_2D(ifgrid, i+isub, i+iadd, j+isub, j+iadd, periodic=periodic)
                    fsumL *= dx*dx
                while IcountL < minpart:
                    IcountS = IcountL
                    if ifgrid is not None:
                        fsumS = fsumL
                    iadd += 1
                    isub -= 1
                    IcountL = sum_from_integral_image_2D(idgrid, i+isub, i+iadd, j+isub, j+iadd, periodic=periodic)
                    IcountL *= dx*dx
                    if ifgrid is not None:
                        fsumL = sum_from_integral_image_2D(ifgrid, i+isub, i+iadd, j+isub, j+iadd, periodic=periodic)
                        fsumL *= dx*dx
                voxelvolS = (iadd - isub - 2)**2
                voxelvolL = (iadd - isub)**2 - voxelvolS
                volS = voxelvolS*(dx**2)
                volL = voxelvolL*(dx**2)
                densL = (IcountL-IcountS)/volL
                inS = IcountS
                inL = minpart - inS
                vol = volS + inL/densL
                if ifgrid is not None:
                    weiS = inS/minpart
                    weiL = inL/minpart
                    if inS == 0.:
                        fgridVEB[i,j] = (fsumL-fsumS)/(IcountL-IcountS)
                    else:
                        fgridVEB[i,j] = weiS*fsumS/IcountS
                        fgridVEB[i,j] += weiL*(fsumL-fsumS)/(IcountL-IcountS)
            vgrid[i,j] = vol
    if ifgrid is None:
        return vgrid
    else:
        return fgridVEB


def vebfe2D(boxsize, ngrid, x, y, minpart=1, w=None, f=None, periodic=True):
    """
    The Volume Enclosing Box Field Estimation (VEBFE) method. This method is similar to a k-Nearest Neighbour 
    method although performed on a grid for speed.

    Parameters
    ----------
    boxsize : float
        Size of the 2D grid.
    ngrid : int
       Grid size.
    x, y : array_like
        X and Y coordinates of the points.
    minpart : int, optional
        Minimum number of particles.
    w : array_like, optional
        Weights for the points, if None assumed to be unitary for all.
    f : array_like
        Field values for the points, if None it is assumed the intention is to compute density.
    periodic : bool, optional
        Periodic boundary condition.
        
    Returns
    -------
    dgridVEB : array_like
        VEB density estimation.
    fgridVEB : array_like
        If f is not None then the field is estimated via VEB.
    """
    if w is None:
        w = np.ones(len(x))
    dgrid = p2g.part2grid2D(x, y, w, boxsize, ngrid, method='NGP', periodic=True, origin=0.)
    if f is not None:
        fgrid = p2g.part2grid2D(x, y, f, boxsize, ngrid, method='NGP', periodic=True, origin=0.)
    idgrid = integral_image.get_integral_image2D(dgrid)
    idgrid = idgrid.astype(np.float64)
    if f is not None:
        ifgrid = integral_image.get_integral_image2D(fgrid)
        ifgrid = ifgrid.astype(np.float64)
    if f is None:
        vgrid = get_volume_enclosing_box(boxsize, ngrid, dgrid, idgrid, minpart, periodic=periodic)
        dgridVEB = minpart/vgrid
        return dgridVEB
    else:
        fgridVEB = get_volume_enclosing_box(boxsize, ngrid, dgrid, idgrid, minpart, periodic=periodic, ifgrid=ifgrid)
        return fgridVEB
