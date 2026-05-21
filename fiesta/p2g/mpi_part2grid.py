import numpy as np
import shift

from typing import List, Union
from .. import coords
from .. import src


def mpi_part2grid2D(
    x: np.ndarray,
    y: np.ndarray,
    f: np.ndarray,
    boxsize: Union[float, List[float]],
    ngrid: Union[int, List[int]],
    MPI: object,
    method: str = "TSC",
    periodic: Union[bool, List[bool]] = True,
    origin: Union[float, List[float]] = 0.0,
) -> np.ndarray:
    """
    Returns the density contrast for the nearest grid point grid assignment.

    Parameters
    ----------
    x : array
        X coordinates of the particle.
    y : array
        Y coordinates of the particle.
    f : array
        Value of each particle to be assigned to the grid.
    boxsize : float or list
        Box size.
    ngrid : int or list
        Grid divisions across one axis.
    MPI : class
        MPIutils MPI class object.
    method : str, optional
        Grid assignment scheme, either 'NGP', 'CIC', 'TSC', 'PCS'.
    periodic : bool or list, optional
        Assign particles with periodic boundaries.
    origin : float, optional
        Origin.

    Returns
    -------
    fgrid : array
        Grid assigned values.
    """
    if x is None:
        data = None
    else:
        data = coords.coord2points([x, y, f])
    data = coords.distribute_points_by_x(data, boxsize, ngrid, origin, MPI)
    x, y, f = data[:,0], data[:,1], data[:,2]
    if np.isscalar(boxsize):
        xlength, ylength = boxsize, boxsize
    else:
        xlength, ylength = boxsize[0], boxsize[1]
    if np.isscalar(origin):
        xmin = origin
        ymin = origin
    else:
        xmin, ymin = origin[0], origin[1]
    if np.isscalar(ngrid):
        nxgrid, nygrid = ngrid, ngrid
    else:
        nxgrid, nygrid = ngrid[0], ngrid[1]
    if np.isscalar(periodic):
        periodx = periodic
        periody = periodic
    else:
        periodx, periody = periodic[0], periodic[1]
    xedges, xgrid = shift.cart.mpi_grid1D(xlength, nxgrid, MPI, origin=xmin)
    xmin, xmax = xedges[0], xedges[-1]
    dx = xedges[1] - xedges[0]
    nxgrid = len(xgrid)

    if method != "NGP":
        xmin -= dx
        xmax += dx
        nxgrid += 2
    xlength = xmax - xmin
    if method == "NGP":
        fgrid = src.part2grid_ngp_2d(
            x, y, f, xlength, ylength, xmin, ymin, nxgrid, nygrid
        )
    elif method == "CIC":
        fgrid = src.part2grid_cic_2d(
            x, y, f, xlength, ylength, xmin, ymin, nxgrid, nygrid, False, periody
        )
    elif method == "TSC":
        fgrid = src.part2grid_tsc_2d(
            x, y, f, xlength, ylength, xmin, ymin, nxgrid, nygrid, False, periody
        )
    elif method == "PCS":
        fgrid = src.part2grid_pcs_2d(
            x, y, f, xlength, ylength, xmin, ymin, nxgrid, nygrid, False, periody
        )
    fgrid = fgrid.reshape(nxgrid, nygrid)
    if method != "NGP":
        fgrid_send_up = MPI.send_up(fgrid[-1])
        fgrid_send_down = MPI.send_down(fgrid[0])
        fgrid = fgrid[1:-1]
        if periodx is True or MPI.rank > 0:
            fgrid[0] += fgrid_send_up
        if periodx is True or MPI.rank < MPI.size - 1:
            fgrid[-1] += fgrid_send_down
    return fgrid


def mpi_part2grid3D(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    f: np.ndarray,
    boxsize: Union[float, List[float]],
    ngrid: Union[int, List[int]],
    MPI: object,
    method: str = "TSC",
    periodic: Union[bool, List[bool]] = True,
    origin: Union[float, List[float]] = 0.0,
) -> np.ndarray:
    """
    Returns the density contrast for the nearest grid point grid assignment.

    Parameters
    ----------
    x : array
        X coordinates of the particle.
    y : array
        Y coordinates of the particle.
    z : array
        Z coordinates of the particle.
    f : array
        Value of each particle to be assigned to the grid.
    boxsize : float or list
        Box size.
    ngrid : int or list
        Grid divisions across one axis.
    MPI : class
        MPIutils MPI class object.
    method : str, optional
        Grid assignment scheme, either 'NGP', 'CIC' or 'TSC'.
    periodic : bool or list, optional
        Assign particles with periodic boundaries.
    origin : float, optional
        Origin.

    Returns
    -------
    fgrid : array
        Grid assigned values.
    """
    if x is None:
        data = None
    else:
        data = coords.coord2points([x, y, z, f])
    data = coords.distribute_points_by_x(data, boxsize, ngrid, origin, MPI)
    x, y, z, f = data[:,0], data[:,1], data[:,2], data[:,3]
    if np.isscalar(boxsize):
        xlength, ylength, zlength = boxsize, boxsize, boxsize
    else:
        xlength, ylength, zlength = boxsize[0], boxsize[1], boxsize[2]
    if np.isscalar(origin):
        xmin = origin
        ymin = origin
        zmin = origin
    else:
        xmin, ymin, zmin = origin[0], origin[1], origin[2]
    if np.isscalar(ngrid):
        nxgrid, nygrid, nzgrid = ngrid, ngrid, ngrid
    else:
        nxgrid, nygrid, nzgrid = ngrid[0], ngrid[1], ngrid[2]
    if np.isscalar(periodic):
        periodx = periodic
        periody = periodic
        periodz = periodic
    else:
        periodx, periody, periodz = periodic[0], periodic[1], periodic[2]
    xedges, xgrid = shift.cart.mpi_grid1D(xlength, nxgrid, MPI, origin=xmin)
    xmin, xmax = xedges[0], xedges[-1]
    dx = xedges[1] - xedges[0]
    nxgrid = len(xgrid)
    if method != "NGP":
        xmin -= dx
        xmax += dx
        nxgrid += 2
    xlength = xmax - xmin
    if method == "NGP":
        fgrid = src.part2grid_ngp_3d(
            x,
            y,
            z,
            f,
            xlength,
            ylength,
            zlength,
            xmin,
            ymin,
            zmin,
            nxgrid,
            nygrid,
            nzgrid,
        )
    elif method == "CIC":
        fgrid = src.part2grid_cic_3d(
            x,
            y,
            z,
            f,
            xlength,
            ylength,
            zlength,
            xmin,
            ymin,
            zmin,
            nxgrid,
            nygrid,
            nzgrid,
            False,
            periody,
            periodz,
        )
    elif method == "TSC":
        fgrid = src.part2grid_tsc_3d(
            x,
            y,
            z,
            f,
            xlength,
            ylength,
            zlength,
            xmin,
            ymin,
            zmin,
            nxgrid,
            nygrid,
            nzgrid,
            False,
            periody,
            periodz,
        )
    elif method == "PCS":
        fgrid = src.part2grid_pcs_3d(
            x,
            y,
            z,
            f,
            xlength,
            ylength,
            zlength,
            xmin,
            ymin,
            zmin,
            nxgrid,
            nygrid,
            nzgrid,
            False,
            periody,
            periodz,
        )
    fgrid = fgrid.reshape(nxgrid, nygrid, nzgrid)
    if method != "NGP":
        fgrid_send_up = MPI.send_up(fgrid[-1])
        fgrid_send_down = MPI.send_down(fgrid[0])
        fgrid = fgrid[1:-1]
        if periodx is True or MPI.rank > 0:
            fgrid[0] += fgrid_send_up
        if periodx is True or MPI.rank < MPI.size - 1:
            fgrid[-1] += fgrid_send_down
    return fgrid
