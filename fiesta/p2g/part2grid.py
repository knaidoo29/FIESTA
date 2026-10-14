import numpy as np

from typing import Union, List

from . import part2grid_ngp
from . import part2grid_cic
from . import part2grid_tsc
from . import part2grid_pcs


def part2grid2D(
    x: np.ndarray,
    y: np.ndarray,
    boxsize: Union[float, List[float]],
    ngrid: Union[int, List[int]],
    f: np.ndarray = None,
    method: str = "TSC",
    periodic: Union[bool, List[bool]] = True,
    origin: Union[float, List[float]] = 0.0,
    dtype: np.dtype = np.float64
) -> np.ndarray:
    """
    Returns the density contrast for the nearest grid point grid assignment.

    Parameters
    ----------
    x : array
        X coordinates of the particle.
    y : array
        Y coordinates of the particle.
    boxsize : float or list
        Box size.
    ngrid : int
        Grid divisions across one axis.
    f : array
        Value of each particle to be assigned to the grid.
    method : str, optional
        Grid assignment scheme, either 'NGP', 'CIC', 'TSC' or 'PCS'.
    periodic : bool or list, optional
        Assign particles with periodic boundaries.
    origin : float or list, optional
        Origin.
    dtype : numpy dtype, optional
        Data type used for the output grid. Default is ``np.float64``.
    
    Returns
    -------
    fgrid : array
        Grid assigned values.
    """
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
        nxgrid, nygrid = int(ngrid), int(ngrid)
    else:
        nxgrid, nygrid = int(ngrid[0]), int(ngrid[1])
    if np.isscalar(periodic):
        periodx = periodic
        periody = periodic
    else:
        periodx, periody = periodic[0], periodic[1]
    if method == "NGP":
        if f is None:
            fgrid = part2grid_ngp.part2grid_ngp_2d_unit(
                x, y, xlength, ylength, xmin, ymin, nxgrid, nygrid, dtype=dtype
            )
        else:
            fgrid = part2grid_ngp.part2grid_ngp_2d(
                x, y, f, xlength, ylength, xmin, ymin, nxgrid, nygrid, dtype=dtype
            )
    elif method == "CIC":
        if f is None:
            fgrid = part2grid_cic.part2grid_cic_2d_unit(
                x, y, xlength, ylength, xmin, ymin, nxgrid, nygrid, periodx, periody, dtype=dtype
            )
        else:
            fgrid = part2grid_cic.part2grid_cic_2d(
                x, y, f, xlength, ylength, xmin, ymin, nxgrid, nygrid, periodx, periody, dtype=dtype
            )
    elif method == "TSC":
        if f is None:
            fgrid = part2grid_tsc.part2grid_tsc_2d_unit(
                x, y, xlength, ylength, xmin, ymin, nxgrid, nygrid, periodx, periody, dtype=dtype
            )
        else:
            fgrid = part2grid_tsc.part2grid_tsc_2d(
                x, y, f, xlength, ylength, xmin, ymin, nxgrid, nygrid, periodx, periody, dtype=dtype
            )
    elif method == "PCS":
        if f is None:
            fgrid = part2grid_pcs.part2grid_pcs_2d_unit(
                x, y, xlength, ylength, xmin, ymin, nxgrid, nygrid, periodx, periody, dtype=dtype
            )
        else:
            fgrid = part2grid_pcs.part2grid_pcs_2d(
                x, y, f, xlength, ylength, xmin, ymin, nxgrid, nygrid, periodx, periody, dtype=dtype
            )
    return fgrid


def part2grid3D(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    boxsize: Union[float, List[float]],
    ngrid: Union[float, List[float]],
    f: np.ndarray = None,
    method: str = "TSC",
    periodic: Union[bool, List[bool]] = True,
    origin: Union[float, List[float]] = 0.0,
    dtype: np.dtype = np.float64
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
    boxsize : float or list
        Box size.
    ngrid : int or list
        Grid divisions across one axis.
    f : array, optional
        Value of each particle to be assigned to the grid.
    method : str, optional
        Grid assignment scheme, either 'NGP', 'CIC', 'TSC' or 'PCS'.
    periodic : bool or list, optional
        Assign particles with periodic boundaries.
    origin : float or list, optional
        Origin.
    dtype : numpy dtype, optional
        Data type used for the output grid. Default is ``np.float64``.

    Returns
    -------
    fgrid : array
        Grid assigned values.
    """
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
        nxgrid, nygrid, nzgrid = int(ngrid), int(ngrid), int(ngrid)
    else:
        nxgrid, nygrid, nzgrid = int(ngrid[0]), int(ngrid[1]), int(ngrid[2])
    if np.isscalar(periodic):
        periodx = periodic
        periody = periodic
        periodz = periodic
    else:
        periodx, periody, periodz = periodic[0], periodic[1], periodic[2]
    if method == "NGP":
        if f is None:
            fgrid = part2grid_ngp.part2grid_ngp_3d_unit(
                x, y, z, xlength, ylength, zlength, xmin, ymin, zmin, 
                nxgrid, nygrid, nzgrid, dtype=dtype
            )
        else:
            fgrid = part2grid_ngp.part2grid_ngp_3d(
                x, y, z, f, xlength, ylength, zlength, xmin, ymin, zmin, 
                nxgrid, nygrid, nzgrid, dtype=dtype
            )
    elif method == "CIC":
        if f is None:
            fgrid = part2grid_cic.part2grid_cic_3d_unit(
                x, y, z, xlength, ylength, zlength, xmin, ymin, zmin,
                nxgrid, nygrid, nzgrid, periodx, periody, periodz, dtype=dtype
            )
        else:
            fgrid = part2grid_cic.part2grid_cic_3d(
                x, y, z, f, xlength, ylength, zlength, xmin, ymin, zmin,
                nxgrid, nygrid, nzgrid, periodx, periody, periodz, dtype=dtype
            )
    elif method == "TSC":
        if f is None:
            fgrid = part2grid_tsc.part2grid_tsc_3d_unit(
                x, y, z, xlength, ylength, zlength, xmin, ymin, zmin,
                nxgrid, nygrid, nzgrid, periodx, periody, periodz, dtype=dtype
            )
        else:
            fgrid = part2grid_tsc.part2grid_tsc_3d(
                x, y, z, f, xlength, ylength, zlength, xmin, ymin, zmin,
                nxgrid, nygrid, nzgrid, periodx, periody, periodz, dtype=dtype
            )
    elif method == "PCS":
        if f is None:
            fgrid = part2grid_pcs.part2grid_pcs_3d_unit(
                x, y, z, xlength, ylength, zlength, xmin, ymin, zmin,
                nxgrid, nygrid, nzgrid, periodx, periody, periodz, dtype=dtype
            )
        else:
            fgrid = part2grid_pcs.part2grid_pcs_3d(
                x, y, z, f, xlength, ylength, zlength, xmin, ymin, zmin,
                nxgrid, nygrid, nzgrid, periodx, periody, periodz, dtype=dtype
            )
    return fgrid
