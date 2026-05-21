import numpy as np

import shift

from .. import coords
from .. import src

from typing import Union, List


def trilinear(
    fgrid: np.ndarray,
    boxsize: Union[float, List[float]],
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    origin: Union[float, List[float]] = 0.0,
    fill_value: float = np.nan,
    periodic: bool = True,
) -> np.ndarray:
    """Trilinear interpolation from a 3D grid defined in box of [0., boxsize].

    Parameter
    ---------
    fgrid : array
        Field values on a 3D grid.
    boxsize : float or list
        Box size in one or all axes.
    x : array
        x coordinate values.
    y : array
        y coordinate values.
    z : array
        z coordinate values.
    origin : float or list, optional
        Origin for the axes.
    fill_value : float, optional
        Fill outside boundary values.
    periodic : bool, optional
        Determines whether to interpolate on a periodic grid.

    Returns
    -------
    f : array
        Field interpolation values.
    """
    # determine ngrid from fgrid
    ngrids = np.shape(fgrid)
    if np.isscalar(boxsize):
        xbox = boxsize
        ybox = boxsize
        zbox = boxsize
    else:
        xbox, ybox, zbox = boxsize[0], boxsize[1], boxsize[2]
    if np.isscalar(origin):
        _x = np.copy(x) - origin
        _y = np.copy(y) - origin
        _z = np.copy(z) - origin
    else:
        _x = np.copy(x) - origin[0]
        _y = np.copy(y) - origin[1]
        _z = np.copy(z) - origin[2]
    # check if particles are inside the box
    cond = np.where(
        (_x >= 0.0)
        & (_x < xbox)
        & (_y >= 0.0)
        & (_y < ybox)
        & (_z >= 0.0)
        & (_z < zbox)
    )[0]
    if len(cond) == len(_x):
        # All particles are within the boundaries so no boundary management is necessary.
        npart = len(_x)
        if np.isscalar(periodic):
            if periodic == True:
                f = src.trilinear_periodic(
                    fgrid.flatten(),
                    _x,
                    _y,
                    _z,
                    xbox,
                    ybox,
                    zbox,
                    ngrids[0],
                    ngrids[1],
                    ngrids[2],
                )
            else:
                f = src.trilinear_nonperiodic(
                    fgrid.flatten(),
                    _x,
                    _y,
                    _z,
                    xbox,
                    ybox,
                    zbox,
                    ngrids[0],
                    ngrids[1],
                    ngrids[2],
                )
        else:
            if periodic[0] is True:
                perix = 1
            else:
                perix = 0
            if periodic[1] is True:
                periy = 1
            else:
                periy = 0
            if periodic[2] is True:
                periz = 1
            else:
                periz = 0
            f = src.trilinear_axisperiodic(
                fgrid.flatten(),
                _x,
                _y,
                _z,
                xbox,
                ybox,
                zbox,
                perix,
                periy,
                periz,
                ngrids[0],
                ngrids[1],
                ngrids[2],
            )
    else:
        # Some particles are outside the boundary.
        # create a mask for in and outside the box
        mask = np.zeros(len(_x))
        # assign particles in the boundary a binary mask of 1.
        mask[cond] = 1.0
        # find trilinear interpolation for points inside the boundary.
        npart = len(x[cond])
        f = np.zeros(len(_x))
        if np.isscalar(periodic):
            if periodic == True:
                f[cond] = src.trilinear_periodic(
                    fgrid.flatten(),
                    _x[cond],
                    _y[cond],
                    _z[cond],
                    xbox,
                    ybox,
                    zbox,
                    ngrids[0],
                    ngrids[1],
                    ngrids[2],
                )
            else:
                f[cond] = src.trilinear_nonperiodic(
                    fgrid.flatten(),
                    _x[cond],
                    _y[cond],
                    _z[cond],
                    xbox,
                    ybox,
                    zbox,
                    ngrids[0],
                    ngrids[1],
                    ngrids[2],
                )
        else:
            if periodic[0] is True:
                perix = 1
            else:
                perix = 0
            if periodic[1] is True:
                periy = 1
            else:
                periy = 0
            if periodic[2] is True:
                periz = 1
            else:
                periz = 0
            f[cond] = src.trilinear_axisperiodic(
                fgrid.flatten(),
                _x[cond],
                _y[cond],
                _z[cond],
                xbox,
                ybox,
                zbox,
                perix,
                periy,
                periz,
                ngrids[0],
                ngrids[1],
                ngrids[2],
            )
        # fill outside boundary with fill values.
        cond = np.where(mask == 0.0)[0]
        f[cond] = fill_value
    return f



def mpi_trilinear(
    fgrid: np.ndarray,
    ngrid: Union[int, List[int]],
    boxsize: Union[float, List[float]],
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    MPI: object,
    origin: Union[float, List[float]] = 0.0,
    fill_value: float = np.nan,
    periodic: bool = True,
) -> np.ndarray:
    """
    Bilinear interpolation from a 2D grid defined in box of [0., boxsize].

    Parameter
    ---------
    fgrid : array
        Field values on a 2D grid
    ngrid : int or int list
        Grid dimensions.
    boxsize : float or list
        Box size.
    x : array
        x coordinate values.
    y : array
        y coordinate values.
    z : array
        z coordinate values.
    MPI : class object
        shift.mpiutils MPI object.
    origin : float or list, optional
        Origin for x and y coordinates.
    fill_value : float, optional
        Fill outside boundary values.
    periodic : bool, optional
        Determines whether to interpolate on a periodic grid.

    Returns
    -------
    x : array
        x coordinate values, redistributed for specific slab.
    y : array
        y coordinate values, redistributed for specific slab.
    z : array
        z coordinate values, redistributed for specific slab.
    f : array
        Field interpolation values.
    """
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

    # define grid on each axis
    if np.isscalar(periodic):
        xperiodic, yperiodic, zperiodic = periodic, periodic, periodic
    else:
        xperiodic, yperiodic, zperiodic = periodic[0], periodic[1], periodic[2]
    
    # define pixel length across each axis
    dx = xboxsize / nxgrid

    xedges, xgrid = shift.cart.mpi_grid1D(xboxsize, nxgrid, MPI, origin=xorigin)

    _xorigin = xedges[0]
    _xboxsize = xedges[-1] - xedges[0]

    if np.shape(fgrid)[0] == len(xgrid):
        correct_shape = True
    else:
        correct_shape = False

    correct_shapes = MPI.collect([correct_shape])

    if MPI.rank == 0:
        if all(correct_shapes):
            correct_shape = True
        else:
            correct_shape = False
    
    correct_shape = MPI.broadcast(correct_shape)
    
    if correct_shape:

        fgrid_sendup = np.array([MPI.send_up(fgrid[-1])])
        fgrid_senddown = np.array([MPI.send_down(fgrid[0])])

        if MPI.rank == 0:
            if xperiodic:
                fgrid = np.concatenate([fgrid_sendup, fgrid, fgrid_senddown])
                _xorigin -= dx
                _xboxsize += 2*dx
            else:
                fgrid = np.concatenate([fgrid, fgrid_senddown])
                _xboxsize += dx
        elif MPI.rank == MPI.size - 1:
            if xperiodic:
                fgrid = np.concatenate([fgrid_sendup, fgrid, fgrid_senddown])
                _xorigin -= dx
                _xboxsize += 2*dx
            else:
                fgrid = np.concatenate([fgrid_sendup, fgrid])
                _xorigin -= dx
                _xboxsize += dx
        else:
            fgrid = np.concatenate([fgrid_sendup, fgrid, fgrid_senddown])
            _xorigin -= dx
            _xboxsize += 2*dx
        
        if x is not None:
            data = coords.xyz2points(x, y, z)
        else: 
            data = None
        data = coords.distribute_points_by_x(data, boxsize, ngrid, origin, MPI)
        x, y, z = coords.points2xyz(data)
        f = trilinear(fgrid, [_xboxsize, yboxsize, zboxsize], x, y, z, [_xorigin, yorigin, zorigin], fill_value=fill_value, periodic=[False, yperiodic, zperiodic])
        return x, y, z, f
    else:
        MPI.mpi_print_zero("ERROR: Shape of fgrid does not match expectation for distributed array")
        return np.nan, np.nan, np.nan, np.nan
        

        


