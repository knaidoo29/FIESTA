import numpy as np

import shift

from .sph2d import SPH2D
from .sph3d import SPH3D

from typing import Union, Tuple, List, Optional


def sph4grid2D(
    x: np.ndarray,
    y: np.ndarray,
    boxsize: Union[float, List[float]],
    ngrid: Union[int, List[int]],
    f: Optional[np.ndarray] = None,
    origin: Union[float, List[float]] = 0.,
    mass: Optional[np.ndarray] = None,
    k: int = 20,
    periodic: bool = True,
    subsampling: int = 1,
    outputgrid: bool = False,
) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """
    Estimates a field on a regular grid based on SPH k neighbours.
    
    Parameters
    ----------
    x, y : array
        Cartesian coordinates of the input points.
    boxsize : float or list of float
        Size of the box in each dimension.
    ngrid : int or list of int
        Number of grid points in each dimension.
    f : array, optional
        Field values at the input points.
    origin : float or list of float, optional
        Origin of the grid.
    mass: array, optional
        Mass of the input points. If not provided, it is assumed to be 1 for all points.
    k : int, optional
        Number of nearest neighbors to consider.
    periodic : bool or list of bool, optional
        Whether to use periodic boundary conditions.
    subsampling : int, optional
        Subsampling factor for the output grid.
    outputgrid : bool, optional
        Whether to output the grid coordinates.

    Returns
    -------
    dgrid : array
        Density field on the regular grid, if f is None.
    fgrid : array
        Estimated field on the regular grid, if f is provided.
    x2d, y2d : array
        Grid coordinates, if outputgrid is True.
    """

    SPH = SPH2D()
    if periodic == True:
        SPH.build_tree(x, y, boxsize=boxsize)
    else:
        SPH.build_tree(x, y)
    
    if mass is None:
        mass = np.ones(len(x))
    
    SPH.setup(k=k, mass=mass)

    # define boxsize on each axis
    if np.isscalar(boxsize):
        xboxsize, yboxsize = boxsize, boxsize
    else:
        xboxsize, yboxsize = boxsize[0], boxsize[1]

    # define boxsize on each axis
    if np.isscalar(origin):
        xorigin, yorigin = origin, origin
    else:
        xorigin, yorigin = origin[0], origin[1]

    # define grid on each axis
    if np.isscalar(ngrid):
        nxgrid, nygrid = ngrid, ngrid
    else:
        nxgrid, nygrid = ngrid[0], ngrid[1]

    # define grid coordinates
    x2d, y2d = shift.cart.grid2D(boxsize, ngrid, origin=[xorigin, yorigin])
    x2d, y2d = x2d.flatten(), y2d.flatten()
    dgrid = np.zeros(len(x2d))
    if f is not None:
        fgrid = np.zeros(len(x2d))

    # define pixel length across each axis
    dx = xboxsize / nxgrid
    dy = yboxsize / nygrid

    dnxgrid, dnygrid = subsampling, subsampling

    # define subsampling points for each pixel
    dx2d, dy2d = shift.cart.grid2D(
        [dx, dy], [dnxgrid, dnygrid], origin=[-dx / 2.0, -dy / 2.0]
    )
    dx2d, dy2d = dx2d.flatten(), dy2d.flatten()

    for i in range(len(dx2d)):
        _x2d = x2d + dx2d[i]
        _y2d = y2d + dy2d[i]

        if f is None:
            dgrid += SPH.get_density(_x2d, _y2d)
        else:
            fgrid += SPH.sph_estimate(_x2d, _y2d, f=f)

    if f is None:
        dgrid /= subsampling**2
        dgrid = dgrid.reshape((nxgrid, nygrid))
        if outputgrid:
            x2d = x2d.reshape((nxgrid, nygrid))
            y2d = y2d.reshape((nxgrid, nygrid))
            return dgrid, x2d, y2d
        else:
            return dgrid
    else:
        fgrid /= subsampling**2
        fgrid = fgrid.reshape((nxgrid, nygrid))
        if outputgrid:
            x2d = x2d.reshape((nxgrid, nygrid))
            y2d = y2d.reshape((nxgrid, nygrid))
            return fgrid, x2d, y2d
        else:
            return fgrid


def sph4grid3D(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    boxsize: Union[float, List[float]],
    ngrid: Union[int, List[int]],
    f: Optional[np.ndarray] = None,
    origin: Union[float, List[float]] = 0.,
    mass: Optional[np.ndarray] = None,
    k: int = 32,
    periodic: bool = True,
    subsampling: int = 1,
    outputgrid: bool = False,
) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """
    Estimates a field on a regular grid based on SPH k neighbours.
    
    Parameters
    ----------
    x, y, z : array
        Cartesian coordinates of the input points.
    boxsize : float or list of float
        Size of the box in each dimension.
    ngrid : int or list of int
        Number of grid points in each dimension.
    f : array, optional
        Field values at the input points.
    origin : float or list of float, optional
        Origin of the grid.
    mass: array, optional
        Mass of the input points. If not provided, it is assumed to be 1 for all points.
    k : int, optional
        Number of nearest neighbors to consider.
    periodic : bool or list of bool, optional
        Whether to use periodic boundary conditions.
    subsampling : int, optional
        Subsampling factor for the output grid.
    outputgrid : bool, optional
        Whether to output the grid coordinates.

    Returns
    -------
    dgrid : array
        Density field on the regular grid, if f is None.
    fgrid : array
        Estimated field on the regular grid, if f is provided.
    x3d, y3d, z3d : array
        Grid coordinates, if outputgrid is True.
    """

    SPH = SPH3D()
    if periodic == True:
        SPH.build_tree(x, y, z, boxsize=boxsize)
    else:
        SPH.build_tree(x, y, z)
    
    if mass is None:
        mass = np.ones(len(x))
    
    SPH.setup(k=k, mass=mass)

    # define boxsize on each axis
    if np.isscalar(boxsize):
        xboxsize, yboxsize, zboxsize = boxsize, boxsize, boxsize
    else:
        xboxsize, yboxsize, zboxsize = boxsize[0], boxsize[1], boxsize[2]

    # define boxsize on each axis
    if np.isscalar(origin):
        xorigin, yorigin, zorigin = origin, origin, origin
    else:
        xorigin, yorigin, zorigin = origin[0], origin[1], origin[2]

    # define grid on each axis
    if np.isscalar(ngrid):
        nxgrid, nygrid, nzgrid = ngrid, ngrid, ngrid
    else:
        nxgrid, nygrid, nzgrid = ngrid[0], ngrid[1], ngrid[2]

    # define grid coordinates
    x3d, y3d, z3d = shift.cart.grid3D(boxsize, ngrid, origin=[xorigin, yorigin, zorigin])
    x3d, y3d, z3d = x3d.flatten(), y3d.flatten(), z3d.flatten()
    dgrid = np.zeros(len(x3d))
    if f is not None:
        fgrid = np.zeros(len(x3d))

    # define pixel length across each axis
    dx = xboxsize / nxgrid
    dy = yboxsize / nygrid
    dz = zboxsize / nzgrid

    dnxgrid, dnygrid, dnzgrid = subsampling, subsampling, subsampling

    # define subsampling points for each pixel
    dx3d, dy3d, dz3d = shift.cart.grid3D(
        [dx, dy, dz], [dnxgrid, dnygrid, dnzgrid], origin=[-dx / 2.0, -dy / 2.0, -dz / 2.0]
    )
    dx3d, dy3d, dz3d = dx3d.flatten(), dy3d.flatten(), dz3d.flatten()

    for i in range(len(dx3d)):
        _x3d = x3d + dx3d[i]
        _y3d = y3d + dy3d[i]
        _z3d = z3d + dz3d[i]

        if f is None:
            dgrid += SPH.get_density(_x3d, _y3d, _z3d)
        else:
            fgrid += SPH.sph_estimate(_x3d, _y3d, _z3d, f=f)

    if f is None:
        dgrid /= subsampling**3
        dgrid = dgrid.reshape((nxgrid, nygrid, nzgrid))
        if outputgrid:
            x3d = x3d.reshape((nxgrid, nygrid, nzgrid))
            y3d = y3d.reshape((nxgrid, nygrid, nzgrid))
            z3d = z3d.reshape((nxgrid, nygrid, nzgrid))
            return dgrid, x3d, y3d, z3d
        else:
            return dgrid
    else:
        fgrid /= subsampling**3
        fgrid = fgrid.reshape((nxgrid, nygrid, nzgrid))
        if outputgrid:
            x3d = x3d.reshape((nxgrid, nygrid, nzgrid))
            y3d = y3d.reshape((nxgrid, nygrid, nzgrid))
            z3d = z3d.reshape((nxgrid, nygrid, nzgrid))
            return fgrid, x3d, y3d, z3d
        else:
            return fgrid