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
    nx = lenx - 1
    ny = leny - 1
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
                    ixmin1 = ixmin + nx
                    ixmax1 = nx
                    ixmin2 = 0
                    ixmax2 = ixmax

                elif ixmax > nx:
                    ixmin1 = ixmin
                    ixmax1 = nx
                    ixmin2 = 0
                    ixmax2 = ixmax - nx
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
                    iymin1 = iymin + ny
                    iymax1 = ny
                    iymin2 = 0
                    iymax2 = iymax
                elif iymax > ny:
                    iymin1 = iymin
                    iymax1 = ny
                    iymin2 = 0
                    iymax2 = iymax - ny
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
    Get the volume for an enclosing box containing minpart particles.

    Parameters
    ----------
    xboxsize, yboxsize : float
        Physical box size along x and y.
    dgrid : ndarray
        2D density grid.
    idgrid : ndarray
        Integral image of the density grid.
    minpart : int
        Minimum number of particles to enclose.
    xperiodic, yperiodic : bool, optional
        Periodic boundary conditions along x and y.
    ifgrid : ndarray, optional
        Integral image of an auxiliary field. If supplied, return the
        volume-enclosing-box estimate of that field instead of the volume.

    Returns
    -------
    ndarray
        Enclosing volume grid, or auxiliary field estimate if ``ifgrid``
        is supplied.
    """

    tol = 1e-12 * max(1.0, float(minpart))

    ngridx, ngridy = dgrid.shape

    dx = xboxsize / ngridx
    dy = yboxsize / ngridy

    cellvol = dx * dy

    vgrid = np.zeros(
        dgrid.shape,
        dtype=np.float64,
    )

    if ifgrid is not None:
        fgridVEB = np.zeros(
            dgrid.shape,
            dtype=np.float64,
        )

    for i in range(ngridx):

        for j in range(ngridy):

            # Number of particles represented by the central cell.
            counts = dgrid[i, j] * cellvol

            # --------------------------------------------------------
            # Central cell already contains minpart particles.
            # --------------------------------------------------------

            if counts >= minpart - tol:

                # If counts is infinitesimally below minpart because of
                # roundoff, regard the complete cell as the enclosing
                # volume rather than allowing a volume > one cell.
                if counts < minpart:
                    vol = cellvol
                else:
                    vol = minpart / (counts / cellvol)

                if ifgrid is not None:

                    fsum = sum_from_integral_image_2D(
                        ifgrid,
                        i,
                        i + 1,
                        j,
                        j + 1,
                        xperiodic=xperiodic,
                        yperiodic=yperiodic,
                    )

                    fgridVEB[i, j] = (
                        fsum / dgrid[i, j]
                    )

            # --------------------------------------------------------
            # Need to expand the enclosing box.
            # --------------------------------------------------------

            else:

                iadd = 1
                isub = 0

                # ----------------------------------------------------
                # Smaller enclosing box.
                # ----------------------------------------------------

                IcountS = sum_from_integral_image_2D(
                    idgrid,
                    i + isub,
                    i + iadd,
                    j + isub,
                    j + iadd,
                    xperiodic=xperiodic,
                    yperiodic=yperiodic,
                )

                IcountS *= cellvol

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

                    fsumS *= cellvol

                # ----------------------------------------------------
                # Larger enclosing box.
                # ----------------------------------------------------

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

                IcountL *= cellvol

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

                    fsumL *= cellvol

                # ----------------------------------------------------
                # Grow the box until it contains minpart particles,
                # allowing for floating-point roundoff.
                # ----------------------------------------------------

                while IcountL < minpart - tol:

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

                    IcountL *= cellvol

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

                        fsumL *= cellvol

                # ----------------------------------------------------
                # Volume of the smaller box and newly added shell.
                # ----------------------------------------------------

                voxelvolS = (
                    iadd - isub - 2
                ) ** 2

                voxelvolL = (
                    iadd - isub
                ) ** 2 - voxelvolS

                volS = voxelvolS * cellvol
                volL = voxelvolL * cellvol

                # ----------------------------------------------------
                # Determine how much of the outer shell is required.
                # ----------------------------------------------------

                inS = IcountS
                inL = minpart - IcountS

                shell_count = (
                    IcountL - IcountS
                )

                remaining_L = (
                    minpart - IcountL
                )

                # ----------------------------------------------------
                # Case 1:
                # The smaller box itself already contains minpart
                # to numerical precision.
                # ----------------------------------------------------

                if abs(inL) <= tol:

                    vol = volS

                    if ifgrid is not None:

                        if IcountS > tol:
                            fgridVEB[i, j] = (
                                fsumS / IcountS
                            )
                        else:
                            fgridVEB[i, j] = 0.0

                # ----------------------------------------------------
                # Case 2:
                # The larger box contains minpart to numerical
                # precision. Use the complete newly added shell.
                #
                # This avoids trying to interpolate through a shell
                # whose measured count may be comparable to numerical
                # roundoff.
                # ----------------------------------------------------

                elif abs(remaining_L) <= tol:

                    vol = volS + volL

                    if ifgrid is not None:

                        if IcountL > tol:
                            fgridVEB[i, j] = (
                                fsumL / IcountL
                            )
                        else:
                            fgridVEB[i, j] = 0.0

                # ----------------------------------------------------
                # Case 3:
                # Genuine threshold crossing. Interpolate through
                # the newly added shell.
                # ----------------------------------------------------

                else:

                    if shell_count <= 0.0:

                        print(
                            "GRIDSPH SHELL ERROR:",
                            "i =", i,
                            "j =", j,
                            "isub =", isub,
                            "iadd =", iadd,
                            "IcountS =", IcountS,
                            "IcountL =", IcountL,
                            "inL =", inL,
                            "remaining_L =", remaining_L,
                            "shell_count =", shell_count,
                            "tol =", tol,
                            "minpart =", minpart,
                        )

                        raise ValueError(
                            "GridSPH encountered a non-positive "
                            "shell particle count."
                        )

                    densL = (
                        shell_count / volL
                    )

                    vol = (
                        volS
                        + inL / densL
                    )

                    if ifgrid is not None:

                        weiS = (
                            inS / minpart
                        )

                        weiL = (
                            inL / minpart
                        )

                        if inS <= tol:

                            fgridVEB[i, j] = (
                                (fsumL - fsumS)
                                / shell_count
                            )

                        else:

                            fgridVEB[i, j] = (
                                weiS
                                * fsumS
                                / IcountS
                            )

                            fgridVEB[i, j] += (
                                weiL
                                * (fsumL - fsumS)
                                / shell_count
                            )

            vgrid[i, j] = vol

    if ifgrid is None:
        return vgrid

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
    lenx, leny, lenz = np.shape(igrid)
    nx = lenx - 1
    ny = leny - 1
    nz = lenz - 1
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
                    ixmin1 = ixmin + nx
                    ixmax1 = nx
                    ixmin2 = 0
                    ixmax2 = ixmax

                elif ixmax > nx:
                    ixmin1 = ixmin
                    ixmax1 = nx
                    ixmin2 = 0
                    ixmax2 = ixmax - nx
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
                    iymin1 = iymin + ny
                    iymax1 = ny
                    iymin2 = 0
                    iymax2 = iymax
                elif iymax > ny:
                    iymin1 = iymin
                    iymax1 = ny
                    iymin2 = 0
                    iymax2 = iymax - ny
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
                    izmin1 = izmin + nz
                    izmax1 = nz
                    izmin2 = 0
                    izmax2 = izmax
                elif izmax > nz:
                    izmin1 = izmin
                    izmax1 = nz
                    izmin2 = 0
                    izmax2 = izmax - nz
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
    Get the volume for an enclosing box containing minpart particles.

    Parameters
    ----------
    xboxsize, yboxsize, zboxsize : float
        Physical box size along x, y and z.
    dgrid : ndarray
        3D density grid.
    idgrid : ndarray
        Integral image of the density grid.
    minpart : int
        Minimum number of particles to enclose.
    xperiodic, yperiodic, zperiodic : bool, optional
        Periodic boundary conditions along x, y and z.
    ifgrid : ndarray, optional
        Integral image of an auxiliary field. If supplied, return the
        volume-enclosing-box estimate of that field instead of the volume.

    Returns
    -------
    ndarray
        Enclosing volume grid, or auxiliary field estimate if ``ifgrid``
        is supplied.
    """

    tol = 1e-12 * max(1.0, float(minpart))

    ngridx, ngridy, ngridz = dgrid.shape

    dx = xboxsize / ngridx
    dy = yboxsize / ngridy
    dz = zboxsize / ngridz

    cellvol = dx * dy * dz

    vgrid = np.zeros(
        dgrid.shape,
        dtype=np.float64,
    )

    if ifgrid is not None:
        fgridVEB = np.zeros(
            dgrid.shape,
            dtype=np.float64,
        )

    for i in range(ngridx):

        for j in range(ngridy):

            for k in range(ngridz):

                # Number of particles represented by the central cell.
                counts = dgrid[i, j, k] * cellvol

                # ----------------------------------------------------
                # Central cell already contains minpart particles.
                # ----------------------------------------------------

                if counts >= minpart - tol:

                    # If counts is infinitesimally below minpart because
                    # of roundoff, use the complete cell rather than
                    # allowing a volume slightly larger than one cell.
                    if counts < minpart:
                        vol = cellvol
                    else:
                        vol = (
                            minpart
                            / (counts / cellvol)
                        )

                    if ifgrid is not None:

                        fsum = sum_from_integral_image_3D(
                            ifgrid,
                            i,
                            i + 1,
                            j,
                            j + 1,
                            k,
                            k + 1,
                            xperiodic=xperiodic,
                            yperiodic=yperiodic,
                            zperiodic=zperiodic,
                        )

                        fgridVEB[i, j, k] = (
                            fsum
                            / dgrid[i, j, k]
                        )

                # ----------------------------------------------------
                # Need to expand the enclosing box.
                # ----------------------------------------------------

                else:

                    iadd = 1
                    isub = 0

                    # ------------------------------------------------
                    # Smaller enclosing box.
                    # ------------------------------------------------

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

                    IcountS *= cellvol

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

                        fsumS *= cellvol

                    # ------------------------------------------------
                    # Larger enclosing box.
                    # ------------------------------------------------

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

                    IcountL *= cellvol

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

                        fsumL *= cellvol

                    # ------------------------------------------------
                    # Grow the box until it contains minpart particles,
                    # allowing for floating-point roundoff.
                    # ------------------------------------------------

                    while IcountL < minpart - tol:

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

                        IcountL *= cellvol

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

                            fsumL *= cellvol

                    # ------------------------------------------------
                    # Volume of the smaller box and newly added shell.
                    # ------------------------------------------------

                    voxelvolS = (
                        iadd - isub - 2
                    ) ** 3

                    voxelvolL = (
                        iadd - isub
                    ) ** 3 - voxelvolS

                    volS = (
                        voxelvolS
                        * cellvol
                    )

                    volL = (
                        voxelvolL
                        * cellvol
                    )

                    # ------------------------------------------------
                    # Determine how much of the newly added shell is
                    # required to reach minpart.
                    # ------------------------------------------------

                    inS = IcountS
                    inL = minpart - IcountS

                    shell_count = (
                        IcountL - IcountS
                    )

                    remaining_L = (
                        minpart - IcountL
                    )

                    # ------------------------------------------------
                    # Case 1:
                    # Smaller box already contains minpart to
                    # numerical precision.
                    # ------------------------------------------------

                    if abs(inL) <= tol:

                        vol = volS

                        if ifgrid is not None:

                            if IcountS > tol:

                                fgridVEB[i, j, k] = (
                                    fsumS
                                    / IcountS
                                )

                            else:

                                fgridVEB[i, j, k] = 0.0

                    # ------------------------------------------------
                    # Case 2:
                    # Larger box contains minpart to numerical
                    # precision. Use the complete outer shell.
                    # ------------------------------------------------

                    elif abs(remaining_L) <= tol:

                        vol = (
                            volS
                            + volL
                        )

                        if ifgrid is not None:

                            if IcountL > tol:

                                fgridVEB[i, j, k] = (
                                    fsumL
                                    / IcountL
                                )

                            else:

                                fgridVEB[i, j, k] = 0.0

                    # ------------------------------------------------
                    # Case 3:
                    # Genuine threshold crossing. Interpolate through
                    # the newly added shell.
                    # ------------------------------------------------

                    else:

                        if shell_count <= 0.0:

                            print(
                                "GRIDSPH SHELL ERROR:",
                                "i =", i,
                                "j =", j,
                                "k =", k,
                                "isub =", isub,
                                "iadd =", iadd,
                                "IcountS =", IcountS,
                                "IcountL =", IcountL,
                                "inL =", inL,
                                "remaining_L =", remaining_L,
                                "shell_count =", shell_count,
                                "tol =", tol,
                                "minpart =", minpart,
                            )

                            raise ValueError(
                                "GridSPH encountered a non-positive "
                                "shell particle count."
                            )

                        densL = (
                            shell_count
                            / volL
                        )

                        vol = (
                            volS
                            + inL / densL
                        )

                        if ifgrid is not None:

                            weiS = (
                                inS
                                / minpart
                            )

                            weiL = (
                                inL
                                / minpart
                            )

                            if inS <= tol:

                                fgridVEB[i, j, k] = (
                                    (fsumL - fsumS)
                                    / shell_count
                                )

                            else:

                                fgridVEB[i, j, k] = (
                                    weiS
                                    * fsumS
                                    / IcountS
                                )

                                fgridVEB[i, j, k] += (
                                    weiL
                                    * (fsumL - fsumS)
                                    / shell_count
                                )

                vgrid[i, j, k] = vol

    if ifgrid is None:
        return vgrid

    return fgridVEB
