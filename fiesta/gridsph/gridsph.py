import numpy as np

from . import integral_image
from .. import p2g

from .. import src

from typing import Optional, Tuple, Union, List


def gridSPH2D(
    boxsize: Union[float, List[float]],
    ngrid: int,
    x: np.ndarray,
    y: np.ndarray,
    minpart: int = 1,
    w: Optional[np.ndarray] = None,
    f: Optional[np.ndarray] = None,
    periodic: Union[bool, List[bool]] = True,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    The grid SPH method. This method is similar to a k-Nearest Neighbour method although
    performed on a grid for speed.

    Parameters
    ----------
    boxsize : float or list
        Size of the 2D grid.
    ngrid : int or list
       Grid size.
    x, y : array
        X and Y coordinates of the points.
    minpart : int, optional
        Minimum number of particles.
    w : array, optional
        Weights for the points, if None assumed to be unitary for all.
    f : array, optional
        Field values for the points, if None it is assumed the intention is to compute density.
    periodic : bool or list, optional
        Periodic boundary condition.

    Returns
    -------
    dgridSPH : array
        Grid SPH density estimation.
    fgridSPH : array
        If f is not None then the field is estimated via SPH.
    """
    if w is None:
        w = np.ones(len(x))
    dgrid = p2g.part2grid2D(
        x, y, w, boxsize, ngrid, method="NGP", periodic=True, origin=0.0
    )
    if f is not None:
        fgrid = p2g.part2grid2D(
            x, y, f, boxsize, ngrid, method="NGP", periodic=True, origin=0.0
        )
    idgrid = integral_image.get_integral_image2D(dgrid)
    idgrid = idgrid.astype(np.float64)
    if f is not None:
        ifgrid = integral_image.get_integral_image2D(fgrid)
        ifgrid = ifgrid.astype(np.float64)
    if np.isscalar(boxsize):
        xboxsize = boxsize
        yboxsize = boxsize
    else:
        xboxsize = boxsize[0]
        yboxsize = boxsize[1]
    if np.isscalar(periodic):
        xperiodic = periodic
        yperiodic = periodic
    else:
        xperiodic = periodic[0]
        yperiodic = periodic[1]
    if f is None:
        vgrid = src.get_volume_enclosing_box_2D(
            xboxsize,
            yboxsize,
            dgrid,
            idgrid,
            minpart,
            xperiodic=xperiodic,
            yperiodic=yperiodic,
        )
        dgridSPH = minpart / vgrid
        return dgridSPH
    else:
        fgridSPH = src.get_volume_enclosing_box_2D(
            xboxsize,
            yboxsize,
            dgrid,
            idgrid,
            minpart,
            xperiodic=xperiodic,
            yperiodic=yperiodic,
            ifgrid=ifgrid,
        )
        return fgridSPH


def gridSPH3D(
    boxsize: Union[float, List[float]],
    ngrid: int,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    minpart: int = 1,
    w: Optional[np.ndarray] = None,
    f: Optional[np.ndarray] = None,
    periodic: Union[bool, List[bool]] = True,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    The grid SPH method. This method is similar to a k-Nearest Neighbour method although 
    performed on a grid for speed.

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
    dgridSPH : array
        Grid SPH density estimation.
    fgridSPH : array
        If f is not None then the field is estimated via SPH.
    """
    if w is None:
        w = np.ones(len(x))
    dgrid = p2g.part2grid3D(
        x, y, z, w, boxsize, ngrid, method="NGP", periodic=True, origin=0.0
    )
    if f is not None:
        fgrid = p2g.part2grid3D(
            x, y, z, f, boxsize, ngrid, method="NGP", periodic=True, origin=0.0
        )
    idgrid = integral_image.get_integral_image3D(dgrid)
    idgrid = idgrid.astype(np.float64)
    if f is not None:
        ifgrid = integral_image.get_integral_image3D(fgrid)
        ifgrid = ifgrid.astype(np.float64)
    if np.isscalar(boxsize):
        xboxsize = boxsize
        yboxsize = boxsize
        zboxsize = boxsize
    else:
        xboxsize = boxsize[0]
        yboxsize = boxsize[1]
        zboxsize = boxsize[2]
    if np.isscalar(periodic):
        xperiodic = periodic
        yperiodic = periodic
        zperiodic = periodic
    else:
        xperiodic = periodic[0]
        yperiodic = periodic[1]
        zperiodic = periodic[2]
    if f is None:
        vgrid = src.get_volume_enclosing_box_3D(
            xboxsize,
            yboxsize,
            zboxsize,
            dgrid,
            idgrid,
            minpart,
            xperiodic=xperiodic,
            yperiodic=yperiodic,
            zperiodic=zperiodic,
        )
        dgridSPH = minpart / vgrid
        return dgridSPH
    else:
        fgridSPH = src.get_volume_enclosing_box_3D(
            xboxsize,
            yboxsize,
            zboxsize,
            dgrid,
            idgrid,
            minpart,
            xperiodic=xperiodic,
            yperiodic=yperiodic,
            zperiodic=zperiodic,
            ifgrid=ifgrid,
        )
        return fgridSPH
