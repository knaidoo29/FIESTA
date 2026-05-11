import numpy as np
from numba import njit

from . import integral_image
from .. import p2g

from typing import Optional, Tuple


@njit
def sum_from_integral_image_3D(
    igrid: np.ndarray, ixmin: int, ixmax: int, iymin: int, iymax: int, izmin: int, izmax: int, periodic: bool = True
) -> np.ndarray:
    """
    Quick summations from an integral image in 3D.

    Parameters
    ----------
    igrid : array
        Integral image in 3D.
    ixmin : int
        Minimum grid value along the x grid.
    ixmax : int
        Maximum grid value along the x grid.
    iymin : int
        Minimum grid value along the y grid.
    iymax : int
        Maximum grid value along the y grid.
    izmin : int
        Minimum grid value along the z grid.
    izmax : int
        Maximum grid value along the z grid.
    periodic : bool, optional
        Periodic boundary condition.
    
    Yields
    ------
    The sum of the integral image within a rectangle.
    """
    if ixmin >= 0 and ixmax < len(igrid[0]) and iymin >= 0 and iymax < len(igrid[:,0]) and izmin >= 0 and izmax < len(igrid[:, :, 0]):
        isum = igrid[ixmax,iymax,izmax] - igrid[ixmax,iymax,izmin]
        isum += - igrid[ixmin, iymax, izmax] + igrid[ixmin, iymax, izmin]
        isum += - igrid[ixmax, iymin, izmax] + igrid[ixmax, iymin, izmin] + igrid[ixmin, iymin, izmax] - igrid[ixmin, iymin, izmin]
        return isum    
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
            if izmin < 0:
                izmin = 0
            if izmax >= len(igrid[:, :, 0]):
                izmax = len(igrid[:, :, 0])-1
            isum = igrid[ixmax,iymax,izmax] - igrid[ixmax,iymax,izmin]
            isum += - igrid[ixmin, iymax, izmax] + igrid[ixmin, iymax, izmin]
            isum += - igrid[ixmax, iymin, ixmax] + igrid[ixmax, iymin, izmin] + igrid[ixmin, iymin, izmax] - igrid[ixmin, iymin, izmin]
            return isum    
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
            
            if izmin < 0:
                izmin1 = izmin+len(igrid[:, :, 0])
                izmax1 = len(igrid[:, :, 0])-1
                izmin2 = 0
                izmax2 = izmax
            elif izmax > len(igrid[:, :, 0])-1:
                izmin1 = izmin
                izmax1 = len(igrid[:, :, 0])-1
                izmin2 = 0
                izmax2 = izmax-len(igrid[:, :, 0])
            else:
                izmin1 = -2*len(igrid[:, :, 0])
                izmax1 = -2*len(igrid[:, :, 0])
                izmin2 = -2*len(igrid[:, :, 0])
                izmax2 = -2*len(igrid[:, :, 0])

            x_split = ixmin1 != -2*len(igrid[0])
            y_split = iymin1 != -2*len(igrid[:,0])
            z_split = izmin1 != -2*len(igrid[:, :, 0])

            if not x_split and not y_split and not z_split:
                ixmins = [ixmin]
                ixmaxs = [ixmax]
                iymins = [iymin]
                iymaxs = [iymax]
                izmins = [izmin]
                izmaxs = [izmax]
            elif x_split and not y_split and not z_split:
                ixmins = [ixmin1, ixmin2]
                ixmaxs = [ixmax1, ixmax2]
                iymins = [iymin, iymin]
                iymaxs = [iymax, iymax]
                izmins = [izmin, izmin]
                izmaxs = [izmax, izmax]
            elif not x_split and y_split and not z_split:
                ixmins = [ixmin, ixmin]
                ixmaxs = [ixmax, ixmax]
                iymins = [iymin1, iymin2]
                iymaxs = [iymax1, iymax2]
                izmins = [izmin, izmin]
                izmaxs = [izmax, izmax]
            elif not x_split and not y_split and z_split:
                ixmins = [ixmin, ixmin]
                ixmaxs = [ixmax, ixmax]
                iymins = [iymin, iymin]
                iymaxs = [iymax, iymax]
                izmins = [izmin1, izmin2]
                izmaxs = [izmax1, izmax2]
            elif x_split and y_split and not z_split:
                ixmins = [ixmin1, ixmin1, ixmin2, ixmin2]
                ixmaxs = [ixmax1, ixmax1, ixmax2, ixmax2]
                iymins = [iymin1, iymin2, iymin1, iymin2]
                iymaxs = [iymax1, iymax2, iymax1, iymax2]
                izmins = [izmin, izmin, izmin, izmin]
                izmaxs = [izmax, izmax, izmax, izmax]
            elif x_split and not y_split and z_split:
                ixmins = [ixmin1, ixmin1, ixmin2, ixmin2]
                ixmaxs = [ixmax1, ixmax1, ixmax2, ixmax2]
                iymins = [iymin, iymin, iymin, iymin]
                iymaxs = [iymax, iymax, iymax, iymax]
                izmins = [izmin1, izmin2, izmin1, izmin2]
                izmaxs = [izmax1, izmax2, izmax1, izmax2]
            elif not x_split and y_split and z_split:
                ixmins = [ixmin, ixmin, ixmin, ixmin]
                ixmaxs = [ixmax, ixmax, ixmax, ixmax]
                iymins = [iymin1, iymin1, iymin2, iymin2]
                iymaxs = [iymax1, iymax1, iymax2, iymax2]
                izmins = [izmin1, izmin2, izmin1, izmin2]
                izmaxs = [izmax1, izmax2, izmax1, izmax2]
            else:  # x_split and y_split and z_split
                ixmins = [ixmin1, ixmin1, ixmin1, ixmin1, ixmin2, ixmin2, ixmin2, ixmin2]
                ixmaxs = [ixmax1, ixmax1, ixmax1, ixmax1, ixmax2, ixmax2, ixmax2, ixmax2]
                iymins = [iymin1, iymin1, iymin2, iymin2, iymin1, iymin1, iymin2, iymin2]
                iymaxs = [iymax1, iymax1, iymax2, iymax2, iymax1, iymax1, iymax2, iymax2]
                izmins = [izmin1, izmin2, izmin1, izmin2, izmin1, izmin2, izmin1, izmin2]
                izmaxs = [izmax1, izmax2, izmax1, izmax2, izmax1, izmax2, izmax1, izmax2]
            isum = 0
            for i in range(0, len(ixmins)):
                isum += igrid[ixmaxs[i],iymaxs[i],izmaxs[i]] - igrid[ixmaxs[i],iymaxs[i],izmins[i]]
                isum += - igrid[ixmins[i], iymaxs[i], izmaxs[i]] + igrid[ixmins[i], iymaxs[i], izmins[i]]
                isum += - igrid[ixmaxs[i], iymins[i], izmaxs[i]] + igrid[ixmaxs[i], iymins[i], izmins[i]] + igrid[ixmins[i], iymins[i], izmaxs[i]] - igrid[ixmins[i], iymins[i], izmins[i]]
            return isum


@njit
def get_volume_enclosing_box_3D(
    boxsize: float, ngrid: int, dgrid: np.ndarray, idgrid: np.ndarray, minpart: int, periodic: bool = True, 
    ifgrid: Optional[np.ndarray] = None
) -> np.ndarray:
    """
    Get the volume for an enclosing box containing minpart number of particles.

    Parameters
    ----------
    boxsize : float
        Size of the 3D grid.
    ngrid : int
       Grid size.
    dgrid : array
        The 3D density field grid.
    idgrid : array
        The 3D density field integral image.
    minpart : int
        Minimum number of particles.
    periodic : bool, optional
        Periodic boundary condition. 
    ifgrid : array
        Field integral image. If not None, then the field will be estimate via 
        the volume enclosing box method.
    
    Returns
    -------
    vgrid : array
        The 3D volume enclosing box for density estimation.
    """
    ngrid = len(dgrid)
    dx = boxsize / ngrid
    vgrid = np.zeros(np.shape(dgrid), dtype=np.float64)
    if ifgrid is not None:
        fgridVEB = np.zeros(np.shape(dgrid), dtype=np.float64)
    for i in range(0, np.shape(dgrid)[0]):
        for j in range(0, np.shape(dgrid)[1]):
            for k in range(0, np.shape(dgrid)[2]):
                counts = dgrid[i,j,k]*dx*dx*dx
                if counts >= minpart:
                    vol = minpart/(counts/(dx**3))
                    if ifgrid is not None:
                        iadd = 1
                        isub = 0
                        fgridVEB[i,j,k] = sum_from_integral_image_3D(ifgrid, i+isub, i+iadd, j+isub, j+iadd, k+isub, k+iadd, periodic=periodic)/dgrid[i,j]
                else:
                    iadd = 1
                    isub = 0
                    # smaller enclosing volume
                    IcountS = sum_from_integral_image_3D(idgrid, i+isub, i+iadd, j+isub, j+iadd, k+isub, k+iadd, periodic=periodic)
                    IcountS *= dx*dx*dx
                    if ifgrid is not None:
                        fsumS = sum_from_integral_image_3D(ifgrid, i+isub, i+iadd, j+isub, j+iadd, k+isub, k+iadd, periodic=periodic)
                        fsumS *= dx*dx*dx
                    # larger enclosing volume
                    iadd += 1
                    isub -= 1
                    IcountL = sum_from_integral_image_3D(idgrid, i+isub, i+iadd, j+isub, j+iadd, k+isub, k+iadd, periodic=periodic)
                    IcountL *= dx*dx*dx
                    if ifgrid is not None:
                        fsumL = sum_from_integral_image_3D(ifgrid, i+isub, i+iadd, j+isub, j+iadd, k+isub, k+iadd, periodic=periodic)
                        fsumL *= dx*dx*dx
                    while IcountL < minpart:
                        IcountS = IcountL
                        if ifgrid is not None:
                            fsumS = fsumL
                        iadd += 1
                        isub -= 1
                        IcountL = sum_from_integral_image_3D(idgrid, i+isub, i+iadd, j+isub, j+iadd, k+isub, k+iadd, periodic=periodic)
                        IcountL *= dx*dx*dx
                        if ifgrid is not None:
                            fsumL = sum_from_integral_image_3D(ifgrid, i+isub, i+iadd, j+isub, j+iadd, k+isub, k+iadd, periodic=periodic)
                            fsumL *= dx*dx*dx
                    voxelvolS = (iadd - isub - 2)**3
                    voxelvolL = (iadd - isub)**3 - voxelvolS
                    volS = voxelvolS*(dx**3)
                    volL = voxelvolL*(dx**3)
                    densL = (IcountL-IcountS)/volL
                    inS = IcountS
                    inL = minpart - inS
                    vol = volS + inL/densL
                    if ifgrid is not None:
                        weiS = inS/minpart
                        weiL = inL/minpart
                        if inS == 0.:
                            fgridVEB[i,j,k] = (fsumL-fsumS)/(IcountL-IcountS)
                        else:
                            fgridVEB[i,j,k] = weiS*fsumS/IcountS
                            fgridVEB[i,j,k] += weiL*(fsumL-fsumS)/(IcountL-IcountS)
                vgrid[i,j,k] = vol
    if ifgrid is None:
        return vgrid
    else:
        return fgridVEB


def vebfe3D(
    boxsize: float, ngrid: int, x: np.ndarray, y: np.ndarray, z: np.ndarray,minpart: int = 1, w: Optional[np.ndarray] = None, 
    f: Optional[np.ndarray] = None, periodic: bool = True
) -> Tuple[np.ndarray, np.ndarray]:
    """
    The Volume Enclosing Box Field Estimation (VEBFE) method. This method is similar to a k-Nearest Neighbour 
    method although performed on a grid for speed.

    Parameters
    ----------
    boxsize : float
        Size of the 3D grid.
    ngrid : int
       Grid size.
    x, y, z : array
        X, Y, and Z coordinates of the points.
    minpart : int, optional
        Minimum number of particles.
    w : array, optional
        Weights for the points, if None assumed to be unitary for all.
    f : array, optional
        Field values for the points, if None it is assumed the intention is to compute density.
    periodic : bool, optional
        Periodic boundary condition.
        
    Returns
    -------
    dgridVEB : array
        VEB density estimation.
    fgridVEB : array
        If f is not None then the field is estimated via VEB.
    """
    if w is None:
        w = np.ones(len(x))
    dgrid = p2g.part2grid3D(x, y, z, w, boxsize, ngrid, method='NGP', periodic=True, origin=0.)
    if f is not None:
        fgrid = p2g.part2grid3D(x, y, z, f, boxsize, ngrid, method='NGP', periodic=True, origin=0.)
    idgrid = integral_image.get_integral_image3D(dgrid)
    idgrid = idgrid.astype(np.float64)
    if f is not None:
        ifgrid = integral_image.get_integral_image3D(fgrid)
        ifgrid = ifgrid.astype(np.float64)
    if f is None:
        vgrid = get_volume_enclosing_box_3D(boxsize, ngrid, dgrid, idgrid, minpart, periodic=periodic)
        dgridVEB = minpart/vgrid
        return dgridVEB
    else:
        fgridVEB = get_volume_enclosing_box_3D(boxsize, ngrid, dgrid, idgrid, minpart, periodic=periodic, ifgrid=ifgrid)
        return fgridVEB
