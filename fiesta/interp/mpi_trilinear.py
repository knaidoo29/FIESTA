import numpy as np

import shift

from . import interp_trilinear
from .. import coords

from typing import Union, List


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
    check_distributed: bool = True,
    dtype: np.dtype = np.float64
) -> np.ndarray:
    """
    Trilinear interpolation from a 3D grid defined in box of [0., boxsize].

    Parameter
    ---------
    fgrid : array
        Field values on a 3D grid
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
    check_distributed : bool, optional
        If True, checks if the coordinates are distributed correctly across MPI ranks.
    dtype : np.dtype
        Data type for the output array.

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

    # Check that the local grid has the expected shape on every rank.
    correct_shape = fgrid.shape[0] == len(xgrid)

    correct_shape = bool(
        MPI.comm.allreduce(correct_shape, op=MPI.mpi.LAND)
    )
    
    if correct_shape:

        fgrid_sendup = MPI.Send_up(np.ascontiguousarray(fgrid[-1:]))
        fgrid_senddown = MPI.Send_down(np.ascontiguousarray(fgrid[:1]))

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
        
        if check_distributed:
            if x is not None:
                data = coords.xyz2points(x, y, z)
            else: 
                data = None
            data = coords.distribute_points_by_x(data, boxsize, ngrid, origin, MPI)
            x, y, z = coords.points2xyz(data)
        
        f = interp_trilinear.trilinear(
            fgrid, [_xboxsize, yboxsize, zboxsize], x, y, z, [_xorigin, yorigin, zorigin], 
            fill_value=fill_value, periodic=[False, yperiodic, zperiodic], dtype=dtype
        )
    
        return x, y, z, f
    
    else:
        MPI.mpi_print_zero("ERROR: Shape of fgrid does not match expectation for distributed array")
        return np.nan, np.nan, np.nan, np.nan
        

        


