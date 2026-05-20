import numpy as np
import shift

from typing import List, Union, Optional, Tuple

from . import integral_image

from .. import p2g
from .. import src


def mpi_gridSPH2D(
    x: np.ndarray,
    y: np.ndarray,
    boxsize: Union[float, List[float]],
    ngrid: int,
    MPI: object,
    minpart: int = 1,
    w: Optional[np.ndarray] = None,
    f: Optional[np.ndarray] = None,
    periodic: Union[bool, List[bool]] = True,
    buffer_size: int = 2,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    The grid SPH method. This method is similar to a k-Nearest Neighbour method although
    performed on a grid for speed.

    Parameters
    ----------
    x, y : array
        X and Y coordinates of the points.
    boxsize : float or list
        Size of the 2D grid.
    ngrid : int or list
       Grid size.
    MPI : class
        MPIutils MPI class object.
    minpart : int, optional
        Minimum number of particles.
    w : array, optional
        Weights for the points, if None assumed to be unitary for all.
    f : array, optional
        Field values for the points, if None it is assumed the intention is to compute density.
    periodic : bool or list, optional
        Periodic boundary condition.
    buffer_size : float, optional
        Buffer size for MPI communication for padding the grid for the SPH estimation.
        This should be at least the size of the largest SPH kernel in units of the grid cell size.

    Returns
    -------
    dgridSPH : array
        Grid SPH density estimation.
    fgridSPH : array
        If f is not None then the field is estimated via SPH.
    """
    if np.isscalar(boxsize):
        xlength = boxsize
    else:
        xlength = boxsize[0]
    xmin = 0.0
    if np.isscalar(ngrid):
        nxgrid = ngrid
    else:
        nxgrid = ngrid[0]
    xedges, xgrid = shift.cart.mpi_grid1D(xlength, nxgrid, MPI, origin=xmin)
    xmin, xmax = xedges[0], xedges[-1]
    dx = xedges[1] - xedges[0]
    nxgrid = len(xgrid)

    xmin -= buffer_size * dx
    xmax += buffer_size * dx
    nxgrid += 2 * buffer_size

    if w is None:
        w = np.ones(len(x))

    dgrid = p2g.mpi_part2grid2D(
        x, y, w, boxsize, ngrid, MPI, method="NGP", periodic=True, origin=0.0
    )
    if f is not None:
        fgrid = p2g.mpi_part2grid2D(
            x, y, f, boxsize, ngrid, MPI, method="NGP", periodic=True, origin=0.0
        )

    dgrid_send_up = MPI.send_up(dgrid[-buffer_size:])
    dgrid_send_down = MPI.send_down(dgrid[:buffer_size])
    dgrid = np.concatenate([dgrid_send_up, dgrid, dgrid_send_down], axis=0)

    if f is not None:
        fgrid_send_up = MPI.send_up(fgrid[-buffer_size:])
        fgrid_send_down = MPI.send_down(fgrid[:buffer_size])
        fgrid = np.concatenate([fgrid_send_up, fgrid, fgrid_send_down], axis=0)

    idgrid = integral_image.get_integral_image2D(dgrid)
    idgrid = idgrid.astype(np.float64)
    if f is not None:
        ifgrid = integral_image.get_integral_image2D(fgrid)
        ifgrid = ifgrid.astype(np.float64)

    xboxsize = xmax - xmin
    if np.isscalar(boxsize):
        yboxsize = boxsize
    else:
        yboxsize = boxsize[1]

    xperiodic = False
    if np.isscalar(periodic):
        yperiodic = periodic
    else:
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
        dgridSPH = dgridSPH[buffer_size:-buffer_size]
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
        fgridSPH = fgridSPH[buffer_size:-buffer_size]
        return fgridSPH


def mpi_gridSPH3D(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    boxsize: Union[float, List[float]],
    ngrid: int,
    MPI: object,
    minpart: int = 1,
    w: Optional[np.ndarray] = None,
    f: Optional[np.ndarray] = None,
    periodic: Union[bool, List[bool]] = True,
    buffer_size: int = 2,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    The grid SPH method. This method is similar to a k-Nearest Neighbour method although
    performed on a grid for speed.

    Parameters
    ----------
    x, y, z : array
        X, Y, and Z coordinates of the points.
    boxsize : float or list
        Size of the 3D grid.
    ngrid : int or list
       Grid size.
    MPI : class
        MPIutils MPI class object.
    minpart : int, optional
        Minimum number of particles.
    w : array, optional
        Weights for the points, if None assumed to be unitary for all.
    f : array, optional
        Field values for the points, if None it is assumed the intention is to compute density.
    periodic : bool or list, optional
        Periodic boundary condition.
    buffer_size : float, optional
        Buffer size for MPI communication for padding the grid for the SPH estimation.
        This should be at least the size of the largest SPH kernel in units of the grid cell size.

    Returns
    -------
    dgridSPH : array
        Grid SPH density estimation.
    fgridSPH : array
        If f is not None then the field is estimated via SPH.
    """
    if np.isscalar(boxsize):
        xlength = boxsize
    else:
        xlength = boxsize[0]
    xmin = 0.0
    if np.isscalar(ngrid):
        nxgrid = ngrid
    else:
        nxgrid = ngrid[0]
    xedges, xgrid = shift.cart.mpi_grid1D(xlength, nxgrid, MPI, origin=xmin)
    xmin, xmax = xedges[0], xedges[-1]
    dx = xedges[1] - xedges[0]
    nxgrid = len(xgrid)

    xmin -= buffer_size * dx
    xmax += buffer_size * dx
    nxgrid += 2 * buffer_size

    if w is None:
        w = np.ones(len(x))

    dgrid = p2g.mpi_part2grid3D(
        x, y, z, w, boxsize, ngrid, MPI, method="NGP", periodic=True, origin=0.0
    )
    if f is not None:
        fgrid = p2g.mpi_part2grid3D(
            x, y, z, f, boxsize, ngrid, MPI, method="NGP", periodic=True, origin=0.0
        )

    dgrid_send_up = MPI.send_up(dgrid[-buffer_size:])
    dgrid_send_down = MPI.send_down(dgrid[:buffer_size])
    dgrid = np.concatenate([dgrid_send_up, dgrid, dgrid_send_down], axis=0)

    if f is not None:
        fgrid_send_up = MPI.send_up(fgrid[-buffer_size:])
        fgrid_send_down = MPI.send_down(fgrid[:buffer_size])
        fgrid = np.concatenate([fgrid_send_up, fgrid, fgrid_send_down], axis=0)

    idgrid = integral_image.get_integral_image3D(dgrid)
    idgrid = idgrid.astype(np.float64)
    if f is not None:
        ifgrid = integral_image.get_integral_image3D(fgrid)
        ifgrid = ifgrid.astype(np.float64)

    xboxsize = xmax - xmin
    if np.isscalar(boxsize):
        yboxsize = boxsize
        zboxsize = boxsize
    else:
        yboxsize = boxsize[1]
        zboxsize = boxsize[2]

    xperiodic = False
    if np.isscalar(periodic):
        yperiodic = periodic
        zperiodic = periodic
    else:
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
        dgridSPH = dgridSPH[buffer_size:-buffer_size]
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
        fgridSPH = fgridSPH[buffer_size:-buffer_size]
        return fgridSPH
