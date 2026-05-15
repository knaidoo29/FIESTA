import numpy as np
from numba import njit

from typing import Optional


@njit
def sum_from_integral_image_2D(
    igrid: np.ndarray,
    ixmin: int,
    ixmax: int,
    iymin: int,
    iymax: int,
    xperiodic: bool = True,
    yperiodic: bool = True,
) -> np.ndarray:  # pragma: no cover
    """
    Quick summations from an integral image in 2D.

    Parameters
    ----------
    igrid : array
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
    lenx = len(igrid[:, 0])
    leny = len(igrid[0])
    if ixmin >= 0 and ixmax < lenx and iymin >= 0 and iymax < leny:
        return (
            igrid[ixmax, iymax]
            - igrid[ixmin, iymax]
            - igrid[ixmax, iymin]
            + igrid[ixmin, iymin]
        )
    else:
        if xperiodic == False and yperiodic == False:
            if ixmin < 0:
                ixmin = 0
            if ixmax >= lenx:
                ixmax = lenx - 1
            if iymin < 0:
                iymin = 0
            if iymax >= leny:
                iymax = leny - 1
            return (
                igrid[ixmax, iymax]
                - igrid[ixmin, iymax]
                - igrid[ixmax, iymin]
                + igrid[ixmin, iymin]
            )
        else:
            if xperiodic == True:
                if ixmin < 0:
                    ixmin1 = ixmin + lenx
                    ixmax1 = lenx - 1
                    ixmin2 = 0
                    ixmax2 = ixmax
                elif ixmax > lenx - 1:
                    ixmin1 = ixmin
                    ixmax1 = lenx - 1
                    ixmin2 = 0
                    ixmax2 = ixmax - lenx
                else:
                    ixmin1 = -2 * lenx
                    ixmax1 = -2 * lenx
                    ixmin2 = -2 * lenx
                    ixmax2 = -2 * lenx
            else:
                if ixmin < 0:
                    ixmin1 = -2 * lenx
                    ixmax1 = -2 * lenx
                    ixmin2 = 0
                    ixmax2 = ixmax
                elif ixmax > lenx - 1:
                    ixmin1 = ixmin
                    ixmax1 = lenx - 1
                    ixmin2 = -2 * lenx
                    ixmax2 = -2 * lenx
                else:
                    ixmin1 = -2 * lenx
                    ixmax1 = -2 * lenx
                    ixmin2 = -2 * lenx
                    ixmax2 = -2 * lenx

            if yperiodic == True:
                if iymin < 0:
                    iymin1 = iymin + leny
                    iymax1 = leny - 1
                    iymin2 = 0
                    iymax2 = iymax
                elif iymax > leny - 1:
                    iymin1 = iymin
                    iymax1 = leny - 1
                    iymin2 = 0
                    iymax2 = iymax - leny
                else:
                    iymin1 = -2 * leny
                    iymax1 = -2 * leny
                    iymin2 = -2 * leny
                    iymax2 = -2 * leny
            else:
                if iymin < 0:
                    iymin1 = -2 * leny
                    iymax1 = -2 * leny
                    iymin2 = 0
                    iymax2 = iymax
                elif iymax > leny - 1:
                    iymin1 = iymin
                    iymax1 = leny - 1
                    iymin2 = -2 * leny
                    iymax2 = -2 * leny
                else:
                    iymin1 = -2 * leny
                    iymax1 = -2 * leny
                    iymin2 = -2 * leny
                    iymax2 = -2 * leny

            if ixmin1 == -2 * lenx:
                xbox1 = False
            else:
                xbox1 = True
            if ixmax2 == -2 * lenx:
                xbox2 = False            
            else:
                xbox2 = True
            
            if iymin1 == -2 * leny:
                ybox1 = False
            else:                
                ybox1 = True
            if iymax2 == -2 * leny:
                ybox2 = False
            else:                
                ybox2 = True

            if not xbox1 and not xbox2:
                _ixmins = [ixmin]
                _ixmaxs = [ixmax]
            elif xbox1 and not xbox2:
                _ixmins = [ixmin1]
                _ixmaxs = [ixmax1]
            elif not xbox1 and xbox2:
                _ixmins = [ixmin2]
                _ixmaxs = [ixmax2]
            else:
                _ixmins = [ixmin1, ixmin2]
                _ixmaxs = [ixmax1, ixmax2]

            if not ybox1 and not ybox2:
                _iymins = [iymin]
                _iymaxs = [iymax]
            elif ybox1 and not ybox2:
                _iymins = [iymin1]
                _iymaxs = [iymax1]
            elif not ybox1 and ybox2:
                _iymins = [iymin2]
                _iymaxs = [iymax2]
            else:
                _iymins = [iymin1, iymin2]
                _iymaxs = [iymax1, iymax2]

            ixmins = []
            ixmaxs = []
            iymins = []
            iymaxs = []

            for i in range(0, len(_ixmins)):
                for j in range(0, len(_iymins)):
                    ixmins.append(_ixmins[i])
                    ixmaxs.append(_ixmaxs[i])
                    iymins.append(_iymins[j])
                    iymaxs.append(_iymaxs[j])         
            
            isum = 0
            for i in range(0, len(ixmins)):
                isum += (
                    igrid[ixmaxs[i], iymaxs[i]]
                    - igrid[ixmins[i], iymaxs[i]]
                    - igrid[ixmaxs[i], iymins[i]]
                    + igrid[ixmins[i], iymins[i]]
                )
            return isum


@njit
def get_volume_enclosing_box_2D(
    xboxsize: float,
    yboxsize: float,
    dgrid: np.ndarray,
    idgrid: np.ndarray,
    minpart: int,
    xperiodic: bool = True,
    yperiodic: bool = True,
    ifgrid: Optional[np.ndarray] = None,
) -> np.ndarray:  # pragma: no cover
    """
    Get the volume for an enclosing box containing minpart number of particles.

    Parameters
    ----------
    boxsize : float
        Size of the 2D grid.
    ngrid : int
       Grid size.
    dgrid : array
        The 2D density field grid.
    idgrid : array
        The 2D density field integral image.
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
        The 2D volume enclosing box for density estimation.
    """
    ngridx, ngridy = np.shape(dgrid)[0], np.shape(dgrid)[1]
    dx = xboxsize / ngridx
    dy = yboxsize / ngridy
    vgrid = np.zeros(np.shape(dgrid), dtype=np.float64)
    if ifgrid is not None:
        fgridVEB = np.zeros(np.shape(dgrid), dtype=np.float64)
    for i in range(0, np.shape(dgrid)[0]):
        for j in range(0, np.shape(dgrid)[1]):
            counts = dgrid[i, j] * dx * dy
            if counts >= minpart:
                vol = minpart / (counts / (dx * dy))
                if ifgrid is not None:
                    iadd = 1
                    isub = 0
                    fgridVEB[i, j] = (
                        sum_from_integral_image_2D(
                            ifgrid,
                            i + isub,
                            i + iadd,
                            j + isub,
                            j + iadd,
                            xperiodic=xperiodic,
                            yperiodic=yperiodic,
                        )
                        / dgrid[i, j]
                    )
            else:
                iadd = 1
                isub = 0
                # smaller enclosing volume
                IcountS = sum_from_integral_image_2D(
                    idgrid,
                    i + isub,
                    i + iadd,
                    j + isub,
                    j + iadd,
                    xperiodic=xperiodic,
                    yperiodic=yperiodic,
                )
                IcountS *= dx * dy
                if ifgrid is not None:
                    fsumS = sum_from_integral_image_2D(
                        ifgrid,
                        i + isub,
                        i + iadd,
                        j + isub,
                        j + iadd,
                        xperiodic=xperiodic,
                        yperiodic=yperiodic,
                    )
                    fsumS *= dx * dy
                # larger enclosing volume
                iadd += 1
                isub -= 1
                IcountL = sum_from_integral_image_2D(
                    idgrid,
                    i + isub,
                    i + iadd,
                    j + isub,
                    j + iadd,
                    xperiodic=xperiodic,
                    yperiodic=yperiodic,
                )
                IcountL *= dx * dy
                if ifgrid is not None:
                    fsumL = sum_from_integral_image_2D(
                        ifgrid,
                        i + isub,
                        i + iadd,
                        j + isub,
                        j + iadd,
                        xperiodic=xperiodic,
                        yperiodic=yperiodic,
                    )
                    fsumL *= dx * dy
                while IcountL < minpart:
                    IcountS = IcountL
                    if ifgrid is not None:
                        fsumS = fsumL
                    iadd += 1
                    isub -= 1
                    IcountL = sum_from_integral_image_2D(
                        idgrid,
                        i + isub,
                        i + iadd,
                        j + isub,
                        j + iadd,
                        xperiodic=xperiodic,
                        yperiodic=yperiodic,
                    )
                    IcountL *= dx * dy
                    if ifgrid is not None:
                        fsumL = sum_from_integral_image_2D(
                            ifgrid,
                            i + isub,
                            i + iadd,
                            j + isub,
                            j + iadd,
                            xperiodic=xperiodic,
                            yperiodic=yperiodic,
                        )
                        fsumL *= dx * dy
                voxelvolS = (iadd - isub - 2) ** 2
                voxelvolL = (iadd - isub) ** 2 - voxelvolS
                volS = voxelvolS * (dx * dy)
                volL = voxelvolL * (dx * dy)
                densL = (IcountL - IcountS) / volL
                inS = IcountS
                inL = minpart - inS
                vol = volS + inL / densL
                if ifgrid is not None:
                    weiS = inS / minpart
                    weiL = inL / minpart
                    if inS == 0.0:
                        fgridVEB[i, j] = (fsumL - fsumS) / (IcountL - IcountS)
                    else:
                        fgridVEB[i, j] = weiS * fsumS / IcountS
                        fgridVEB[i, j] += weiL * (fsumL - fsumS) / (IcountL - IcountS)
            vgrid[i, j] = vol
    if ifgrid is None:
        return vgrid
    else:
        return fgridVEB


@njit
def sum_from_integral_image_3D(
    igrid: np.ndarray,
    ixmin: int,
    ixmax: int,
    iymin: int,
    iymax: int,
    izmin: int,
    izmax: int,
    xperiodic: bool = True,
    yperiodic: bool = True,
    zperiodic: bool = True,
) -> np.ndarray:  # pragma: no cover
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
    xperiodic : bool, optional
        Periodic boundary condition along the x grid.
    yperiodic : bool, optional
        Periodic boundary condition along the y grid.
    zperiodic : bool, optional
        Periodic boundary condition along the z grid.

    Yields
    ------
    The sum of the integral image within a rectangle.
    """
    (lenx, leny, lenz) = np.shape(igrid)
    if (
        ixmin >= 0
        and ixmax < lenx
        and iymin >= 0
        and iymax < leny
        and izmin >= 0
        and izmax < lenz
    ):
        isum = igrid[ixmax, iymax, izmax] - igrid[ixmax, iymax, izmin]
        isum += -igrid[ixmin, iymax, izmax] + igrid[ixmin, iymax, izmin]
        isum += (
            -igrid[ixmax, iymin, izmax]
            + igrid[ixmax, iymin, izmin]
            + igrid[ixmin, iymin, izmax]
            - igrid[ixmin, iymin, izmin]
        )
        return isum
    else:
        if xperiodic == False and yperiodic == False and zperiodic == False:
            if ixmin < 0:
                ixmin = 0
            if ixmax >= lenx:
                ixmax = lenx - 1
            if iymin < 0:
                iymin = 0
            if iymax >= leny:
                iymax = leny - 1
            if izmin < 0:
                izmin = 0
            if izmax >= lenz:
                izmax = lenz - 1
            isum = igrid[ixmax, iymax, izmax] - igrid[ixmax, iymax, izmin]
            isum += -igrid[ixmin, iymax, izmax] + igrid[ixmin, iymax, izmin]
            isum += (
                -igrid[ixmax, iymin, izmax]
                + igrid[ixmax, iymin, izmin]
                + igrid[ixmin, iymin, izmax]
                - igrid[ixmin, iymin, izmin]
            )
            return isum
        else:
            if xperiodic == True:
                if ixmin < 0:
                    ixmin1 = ixmin + lenx
                    ixmax1 = lenx - 1
                    ixmin2 = 0
                    ixmax2 = ixmax
                elif ixmax > lenx - 1:
                    ixmin1 = ixmin
                    ixmax1 = lenx - 1
                    ixmin2 = 0
                    ixmax2 = ixmax - lenx
                else:
                    ixmin1 = -2 * lenx
                    ixmax1 = -2 * lenx
                    ixmin2 = -2 * lenx
                    ixmax2 = -2 * lenx
            else:
                if ixmin < 0:
                    ixmin1 = -2 * lenx
                    ixmax1 = -2 * lenx
                    ixmin2 = 0
                    ixmax2 = ixmax
                elif ixmax > lenx - 1:
                    ixmin1 = ixmin
                    ixmax1 = lenx - 1
                    ixmin2 = -2 * lenx
                    ixmax2 = -2 * lenx
                else:
                    ixmin1 = -2 * lenx
                    ixmax1 = -2 * lenx
                    ixmin2 = -2 * lenx
                    ixmax2 = -2 * lenx

            if yperiodic == True:
                if iymin < 0:
                    iymin1 = iymin + leny
                    iymax1 = leny - 1
                    iymin2 = 0
                    iymax2 = iymax
                elif iymax > leny - 1:
                    iymin1 = iymin
                    iymax1 = leny - 1
                    iymin2 = 0
                    iymax2 = iymax - leny
                else:
                    iymin1 = -2 * leny
                    iymax1 = -2 * leny
                    iymin2 = -2 * leny
                    iymax2 = -2 * leny
            else:
                if iymin < 0:
                    iymin1 = -2 * leny
                    iymax1 = -2 * leny
                    iymin2 = 0
                    iymax2 = iymax
                elif iymax > leny - 1:
                    iymin1 = iymin
                    iymax1 = leny - 1
                    iymin2 = -2 * leny
                    iymax2 = -2 * leny
                else:
                    iymin1 = -2 * leny
                    iymax1 = -2 * leny
                    iymin2 = -2 * leny
                    iymax2 = -2 * leny

            if zperiodic == True:
                if izmin < 0:
                    izmin1 = izmin + lenz
                    izmax1 = lenz - 1
                    izmin2 = 0
                    izmax2 = izmax
                elif izmax > lenz - 1:
                    izmin1 = izmin
                    izmax1 = lenz - 1
                    izmin2 = 0
                    izmax2 = izmax - lenz
                else:
                    izmin1 = -2 * lenz
                    izmax1 = -2 * lenz
                    izmin2 = -2 * lenz
                    izmax2 = -2 * lenz
            else:
                if izmin < 0:
                    izmin1 = -2 * lenz
                    izmax1 = -2 * lenz
                    izmin2 = 0
                    izmax2 = izmax
                elif izmax > lenz - 1:
                    izmin1 = izmin
                    izmax1 = lenz - 1
                    izmin2 = -2 * lenz
                    izmax2 = -2 * lenz
                else:
                    izmin1 = -2 * lenz
                    izmax1 = -2 * lenz
                    izmin2 = -2 * lenz
                    izmax2 = -2 * lenz

            if ixmin1 == -2 * lenx:
                xbox1 = False
            else:
                xbox1 = True
            if ixmax2 == -2 * lenx:
                xbox2 = False            
            else:
                xbox2 = True
            
            if iymin1 == -2 * leny:
                ybox1 = False
            else:                
                ybox1 = True
            if iymax2 == -2 * leny:
                ybox2 = False
            else:                
                ybox2 = True

            if izmin1 == -2 * lenz:
                zbox1 = False
            else:         
                zbox1 = True
            if izmax2 == -2 * lenz:
                zbox2 = False
            else:                
                zbox2 = True

            if not xbox1 and not xbox2:
                _ixmins = [ixmin]
                _ixmaxs = [ixmax]
            elif xbox1 and not xbox2:
                _ixmins = [ixmin1]
                _ixmaxs = [ixmax1]
            elif not xbox1 and xbox2:
                _ixmins = [ixmin2]
                _ixmaxs = [ixmax2]
            else:
                _ixmins = [ixmin1, ixmin2]
                _ixmaxs = [ixmax1, ixmax2]

            if not ybox1 and not ybox2:
                _iymins = [iymin]
                _iymaxs = [iymax]
            elif ybox1 and not ybox2:
                _iymins = [iymin1]
                _iymaxs = [iymax1]
            elif not ybox1 and ybox2:
                _iymins = [iymin2]
                _iymaxs = [iymax2]
            else:
                _iymins = [iymin1, iymin2]
                _iymaxs = [iymax1, iymax2]

            if not zbox1 and not zbox2:  
                _izmins = [izmin]
                _izmaxs = [izmax]
            elif zbox1 and not zbox2:
                _izmins = [izmin1]
                _izmaxs = [izmax1]
            elif not zbox1 and zbox2:
                _izmins = [izmin2]
                _izmaxs = [izmax2]
            else:
                _izmins = [izmin1, izmin2]
                _izmaxs = [izmax1, izmax2]

            ixmins = []
            ixmaxs = []
            iymins = []
            iymaxs = []
            izmins = []
            izmaxs = []

            for i in range(0, len(_ixmins)):
                for j in range(0, len(_iymins)):
                    for k in range(0, len(_izmins)):
                        ixmins.append(_ixmins[i])
                        ixmaxs.append(_ixmaxs[i])
                        iymins.append(_iymins[j])
                        iymaxs.append(_iymaxs[j])
                        izmins.append(_izmins[k])
                        izmaxs.append(_izmaxs[k])            
            isum = 0
            for i in range(0, len(ixmins)):
                isum += (
                    igrid[ixmaxs[i], iymaxs[i], izmaxs[i]]
                    - igrid[ixmaxs[i], iymaxs[i], izmins[i]]
                )
                isum += (
                    -igrid[ixmins[i], iymaxs[i], izmaxs[i]]
                    + igrid[ixmins[i], iymaxs[i], izmins[i]]
                )
                isum += (
                    -igrid[ixmaxs[i], iymins[i], izmaxs[i]]
                    + igrid[ixmaxs[i], iymins[i], izmins[i]]
                )
                isum += (
                    igrid[ixmins[i], iymins[i], izmaxs[i]]
                    - igrid[ixmins[i], iymins[i], izmins[i]]
                )
            return isum


@njit
def get_volume_enclosing_box_3D(
    xboxsize: float,
    yboxsize: float,
    zboxsize: float,
    dgrid: np.ndarray,
    idgrid: np.ndarray,
    minpart: int,
    xperiodic: bool = True,
    yperiodic: bool = True,
    zperiodic: bool = True,
    ifgrid: Optional[np.ndarray] = None,
) -> np.ndarray:  # pragma: no cover
    """
    Get the volume for an enclosing box containing minpart number of particles.

    Parameters
    ----------
    boxsize : float
        Size of the 3D grid.
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
    ngridx, ngridy, ngridz = np.shape(dgrid)[0], np.shape(dgrid)[1], np.shape(dgrid)[2]
    dx = xboxsize / ngridx
    dy = yboxsize / ngridy
    dz = zboxsize / ngridz
    vgrid = np.zeros(np.shape(dgrid), dtype=np.float64)
    if ifgrid is not None:
        fgridVEB = np.zeros(np.shape(dgrid), dtype=np.float64)
    for i in range(0, np.shape(dgrid)[0]):
        for j in range(0, np.shape(dgrid)[1]):
            for k in range(0, np.shape(dgrid)[2]):
                counts = dgrid[i, j, k] * dx * dy * dz
                if counts >= minpart:
                    vol = minpart / (counts / (dx * dy * dz))
                    if ifgrid is not None:
                        iadd = 1
                        isub = 0
                        fgridVEB[i, j, k] = (
                            sum_from_integral_image_3D(
                                ifgrid,
                                i + isub,
                                i + iadd,
                                j + isub,
                                j + iadd,
                                k + isub,
                                k + iadd,
                                xperiodic=xperiodic,
                                yperiodic=yperiodic,
                                zperiodic=zperiodic,
                            )
                            / dgrid[i, j, k]
                        )
                else:
                    iadd = 1
                    isub = 0
                    # smaller enclosing volume
                    IcountS = sum_from_integral_image_3D(
                        idgrid,
                        i + isub,
                        i + iadd,
                        j + isub,
                        j + iadd,
                        k + isub,
                        k + iadd,
                        xperiodic=xperiodic,
                        yperiodic=yperiodic,
                        zperiodic=zperiodic,
                    )
                    IcountS *= dx * dy * dz
                    if ifgrid is not None:
                        fsumS = sum_from_integral_image_3D(
                            ifgrid,
                            i + isub,
                            i + iadd,
                            j + isub,
                            j + iadd,
                            k + isub,
                            k + iadd,
                            xperiodic=xperiodic,
                            yperiodic=yperiodic,
                            zperiodic=zperiodic,
                        )
                        fsumS *= dx * dy * dz
                    # larger enclosing volume
                    iadd += 1
                    isub -= 1
                    IcountL = sum_from_integral_image_3D(
                        idgrid,
                        i + isub,
                        i + iadd,
                        j + isub,
                        j + iadd,
                        k + isub,
                        k + iadd,
                        xperiodic=xperiodic,
                        yperiodic=yperiodic,
                        zperiodic=zperiodic,
                    )
                    IcountL *= dx * dy * dz
                    if ifgrid is not None:
                        fsumL = sum_from_integral_image_3D(
                            ifgrid,
                            i + isub,
                            i + iadd,
                            j + isub,
                            j + iadd,
                            k + isub,
                            k + iadd,
                            xperiodic=xperiodic,
                            yperiodic=yperiodic,
                            zperiodic=zperiodic,
                        )
                        fsumL *= dx * dy * dz
                    while IcountL < minpart:
                        IcountS = IcountL
                        if ifgrid is not None:
                            fsumS = fsumL
                        iadd += 1
                        isub -= 1
                        IcountL = sum_from_integral_image_3D(
                            idgrid,
                            i + isub,
                            i + iadd,
                            j + isub,
                            j + iadd,
                            k + isub,
                            k + iadd,
                            xperiodic=xperiodic,
                            yperiodic=yperiodic,
                            zperiodic=zperiodic,
                        )
                        IcountL *= dx * dy * dz
                        if ifgrid is not None:
                            fsumL = sum_from_integral_image_3D(
                                ifgrid,
                                i + isub,
                                i + iadd,
                                j + isub,
                                j + iadd,
                                k + isub,
                                k + iadd,
                                xperiodic=xperiodic,
                                yperiodic=yperiodic,
                                zperiodic=zperiodic,
                            )
                            fsumL *= dx * dy * dz
                    voxelvolS = (iadd - isub - 2) ** 3
                    voxelvolL = (iadd - isub) ** 3 - voxelvolS
                    volS = voxelvolS * (dx * dy * dz)
                    volL = voxelvolL * (dx * dy * dz)
                    densL = (IcountL - IcountS) / volL
                    inS = IcountS
                    inL = minpart - inS
                    vol = volS + inL / densL
                    if ifgrid is not None:
                        weiS = inS / minpart
                        weiL = inL / minpart
                        if inS == 0.0:
                            fgridVEB[i, j, k] = (fsumL - fsumS) / (IcountL - IcountS)
                        else:
                            fgridVEB[i, j, k] = weiS * fsumS / IcountS
                            fgridVEB[i, j, k] += (
                                weiL * (fsumL - fsumS) / (IcountL - IcountS)
                            )
                vgrid[i, j, k] = vol
    if ifgrid is None:
        return vgrid
    else:
        return fgridVEB
