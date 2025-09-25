import numpy as np

import shift

from .. import coords
from . import dtfe2d, dtfe3d
from . import dtfe4grid

from typing import Union, Tuple, Optional, List
from copy import deepcopy


def mpi_delaunay_density4grid2D(
    x: np.ndarray, y: np.ndarray, boxsize: Union[float, List[float]], ngrid: Union[int, List[int]], 
    MPI: object, origin: Union[float, List[float]] = 0., mass: Optional[np.ndarray] = None, 
    partition: int = 1, periodic : Union[bool, List[bool]] = True, fbuffer: float = 0.5, subsampling: int = 1, 
    outputgrid: bool = False
) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """
    Returns the Delaunay tesselation density on a grid.

    Parameters
    ----------
    x, y : array
        Coordinates of particles.
    boxsize : float or list
        Dimensions of the grid.
    ngrid : int or int list
        Grid dimensions.
    MPI : class object
        shift.mpiutils MPI object.
    origin : float or list, optional
        Origin for grid.
    mass : array, optional
        Mass or mass weights for each the particles.
    partition : int or list, optional
        The number of internal Delaunay tesselations used to compute the density grid.
    periodic : bool or list, optional
        Determines whether periodic boundaries are applied.
    fbuffer : float, optional
        Buffer factor length, used to trim exterior buffer region.
    subsampling : int, optional
        The pixel subsampling rate. Each pixel is evaluated subsampling^2 points
        on a grid within each pixel. This is to ensure each pixel is assigned a
        mean pixel value and not the value at the center.
    outputgrid : bool, optional
        Outputs coordinate grid.
    outputexterior : bool, optional
        Ouput exterior information, including exterior border dtfe and current unassigned pixels.

    Returns
    -------
    dens : ndarray
        Density field values on a grid.
    x2d, y2d : ndarray, optional
        Pixel coordinate points.
    """
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
    
    # define grid on each axis
    if np.isscalar(periodic):
        xperiodic, yperiodic = periodic, periodic
    else:
        xperiodic, yperiodic = periodic[0], periodic[1]
    
    # define pixel length across each axis
    dx = xboxsize/nxgrid
    dy = yboxsize/nygrid

    if np.isscalar(subsampling):
        dnxgrid, dnygrid = subsampling, subsampling
    else:
        dnxgrid, dnygrid = subsampling[0], subsampling[1]
    
    # define subsampling points for each pixel
    dx2d, dy2d = shift.cart.grid2D([dx, dy], [dnxgrid, dnygrid], origin=[-dx/2., -dy/2.])
    dx2d, dy2d = dx2d.flatten(), dy2d.flatten()

    # collapse data
    if x is not None:
        if mass is None:
            data = coords.coord2points([x, y, np.ones(len(x))])
        else:
            data = coords.coord2points([x, y, mass])
    else:
        data = None
    
    # check size of data
    if data is not None:
        lendata = len(data)
    else:
        lendata = 0
    size = MPI.collect(lendata, outlist=True)
    if MPI.rank == 0:
        npart = np.sum(np.array(size))
        MPI.send(npart, tag=11)
    else:
        npart = MPI.recv(0, tag=11)
    MPI.wait()

    # sort coordinates and distribute by coordinate system
    MPI_SBX = coords.MPI_SortByX(MPI)
    MPI_SBX.settings(xboxsize, nxgrid, origin=xorigin)
    MPI_SBX.input(data)
    MPI_SBX.limits4grid()
    data = MPI_SBX.distribute()
    limits = [MPI_SBX.limits[0], MPI_SBX.limits[1], yorigin, yorigin+yboxsize]
    
    # coordinates only inside subbox
    x, y, mass = data[:,0], data[:,1], data[:,2]

    # create grids
    xedges, xgrid = shift.cart.mpi_grid1D(xboxsize, nxgrid, MPI, origin=xorigin)
    yedges, ygrid = shift.cart.grid1D(yboxsize, nygrid, origin=yorigin)
    x2d, y2d = shift.cart.mpi_grid2D(boxsize, ngrid, MPI, origin=origin)
    x2d = x2d.flatten()
    y2d = y2d.flatten()

    # local grid, we will use an underscore to differentiate from global
    _xorigin = xedges[0]
    _xboxsize = xedges[-1] - xedges[0]
    _yorigin = yedges[0]
    _yboxsize = yedges[-1] - yedges[0]

    _nxgrid = len(xgrid)
    _nygrid = len(ygrid)

    dens, exterior_border, pixID, count = dtfe4grid.delaunay_density4grid2D(
        x, y, [_xboxsize, _yboxsize], [_nxgrid, _nygrid], origin=[_xorigin, _yorigin], mass=mass, 
        partition=partition, periodic=[False, False], fbuffer=fbuffer, subsampling=subsampling, 
        outputgrid=False, outputexterior=True, normalise=False, flatten=True
    )

    MPI.wait()

    # get border particles from adjacent slabs with the slab below

    exterior_border_below = {}
    
    exterior_border_below['x'] = MPI.send_up(exterior_border['x'])
    exterior_border_below['y'] = MPI.send_up(exterior_border['y'])
    exterior_border_below['mass'] = MPI.send_up(exterior_border['mass'])
    exterior_border_below['ptype'] = MPI.send_up(exterior_border['ptype'])
    exterior_border_below['f'] = MPI.send_up(exterior_border['f'])
    exterior_border_below['dtfe_mode'] = exterior_border['dtfe_mode']
    exterior_border_below['simplices'] = MPI.send_up(exterior_border['simplices'])
    exterior_border_below['simptypes'] = MPI.send_up(exterior_border['simptypes'])

    # get border particles from adjacent slabs with the slab above
    
    exterior_border_above = {}

    exterior_border_above['x'] = MPI.send_down(exterior_border['x'])
    exterior_border_above['y'] = MPI.send_down(exterior_border['y'])
    exterior_border_above['mass'] = MPI.send_down(exterior_border['mass'])
    exterior_border_above['ptype'] = MPI.send_down(exterior_border['ptype'])
    exterior_border_above['f'] = MPI.send_down(exterior_border['f'])
    exterior_border_above['dtfe_mode'] = exterior_border['dtfe_mode']
    exterior_border_above['simplices'] = MPI.send_down(exterior_border['simplices'])
    exterior_border_above['simptypes'] = MPI.send_down(exterior_border['simptypes'])
    
    # Construct merged Delaunay tesselation

    dtfe = dtfe2d.DelaunayMerger2D()

    if yperiodic:
        jmin, jmax = -1, 2
    else:
        jmin, jmax = 0, 1

    for j in range(jmin, jmax):
        _exterior_border_below = deepcopy(exterior_border_below)
        _exterior_border_below['y'] += j*yboxsize
        if MPI.rank != 0:
            dtfe.add_border(_exterior_border_below)
        else:
            if xperiodic:
                _exterior_border_below['x'] -= xboxsize
                dtfe.add_border(_exterior_border_below)
        
        _exterior_border = deepcopy(exterior_border)
        _exterior_border['y'] += j*yboxsize
        dtfe.add_border(_exterior_border)

        _exterior_border_above = deepcopy(exterior_border_above)
        _exterior_border_above['y'] += j*yboxsize
        if MPI.rank != MPI.size-1:
            dtfe.add_border(_exterior_border_above)
        else:
            if xperiodic:
                _exterior_border_above['x'] += xboxsize
                dtfe.add_border(_exterior_border_above)

    _xbuffer = fbuffer*_xboxsize
    _ybuffer = fbuffer*_yboxsize
    
    boundary = [_xorigin-_xbuffer, _xorigin+_xboxsize+_xbuffer, _yorigin-_ybuffer, _yorigin+_yboxsize+_ybuffer]

    dtfe.run(boundary, apply_filter=True)

    cond = np.where((x2d[pixID] >= boundary[0]) & (x2d[pixID] < boundary[1]) & (y2d[pixID] >= boundary[2]) & (y2d[pixID] < boundary[3]))[0]

    for _dx2d, _dy2d in zip(dx2d, dy2d):
        dtfe_estimate = dtfe.estimate(x2d[pixID[cond]]+_dx2d, y2d[pixID[cond]]+_dy2d)
        cond_isfinite = np.where(np.isfinite(dtfe_estimate) & (count[cond] != len(dx2d)))[0]
        dens[pixID[cond[cond_isfinite]]] += dtfe_estimate[cond_isfinite]
        count[cond[cond_isfinite]] += 1.

    cond2 = np.where(count[cond] == len(dx2d))[0]
    pixID = np.delete(pixID, cond[cond2])
    count = np.delete(count, cond[cond2])

    exterior_border = dtfe.get_border()

    dens /= len(dx2d)

    cond = np.where(count == 0.)[0]
    dens[pixID[cond]] *= len(dx2d)
    cond = np.where(count != 0.)[0]
    dens[pixID[cond]] *= len(dx2d)/count[cond]

    dens = dens.reshape(_nxgrid, _nygrid)

    if outputgrid:
        return dens, x2d, y2d
    else:
        return dens


def mpi_delaunay_field4grid2D(
    x: np.ndarray, y: np.ndarray, f: np.ndarray, boxsize: Union[float, List[float]], ngrid: Union[int, List[int]], 
    MPI: object, origin: Union[float, List[float]] = 0., mass: Optional[np.ndarray] = None, 
    partition: int = 1, periodic : Union[bool, List[bool]] = True, fbuffer: float = 0.5, subsampling: int = 1, 
    outputgrid: bool = False
):
    """
    Returns the Delaunay tesselation field on a grid.

    Parameters
    ----------
    x, y : array
        Coordinates of particles.
    f : array
        Field values to be interpolated.
    boxsize : float or list
        Dimensions of the grid.
    ngrid : int or int list
        Grid dimensions.
    MPI : class object
        shift.mpiutils MPI object.
    origin : float or list, optional
        Origin for grid.
    mass : array, optional
        Mass or mass weights for each the particles.
    partition : int or list, optional
        The number of internal Delaunay tesselations used to compute the field grid.
    periodic : bool or list, optional
        Determines whether periodic boundaries are applied.
    fbuffer : float, optional
        Buffer factor length, used to trim exterior buffer region.
    subsampling : int, optional
        The pixel subsampling rate. Each pixel is evaluated subsampling^2 points
        on a grid within each pixel. This is to ensure each pixel is assigned a
        mean pixel value and not the value at the center.
    outputgrid : bool, optional
        Outputs coordinate grid.
    outputexterior : bool, optional
        Ouput exterior information, including exterior border dtfe and current unassigned pixels.

    Returns
    -------
    field : ndarray
        Density field values on a grid.
    x2d, y2d : ndarray, optional
        Pixel coordinate points.
    """
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
    
    # define grid on each axis
    if np.isscalar(periodic):
        xperiodic, yperiodic = periodic, periodic
    else:
        xperiodic, yperiodic = periodic[0], periodic[1]
    
    # define pixel length across each axis
    dx = xboxsize/nxgrid
    dy = yboxsize/nygrid

    if np.isscalar(subsampling):
        dnxgrid, dnygrid = subsampling, subsampling
    else:
        dnxgrid, dnygrid = subsampling[0], subsampling[1]
    
    # define subsampling points for each pixel
    dx2d, dy2d = shift.cart.grid2D([dx, dy], [dnxgrid, dnygrid], origin=[-dx/2., -dy/2.])
    dx2d, dy2d = dx2d.flatten(), dy2d.flatten()

    # collapse data
    if x is not None:
        if mass is None:
            data = coords.coord2points([x, y, f, np.ones(len(x))])
        else:
            data = coords.coord2points([x, y, f, mass])
    else:
        data = None
    
    # check size of data
    if data is not None:
        lendata = len(data)
    else:
        lendata = 0
    size = MPI.collect(lendata, outlist=True)
    if MPI.rank == 0:
        npart = np.sum(np.array(size))
        MPI.send(npart, tag=11)
    else:
        npart = MPI.recv(0, tag=11)
    MPI.wait()

    # sort coordinates and distribute by coordinate system
    MPI_SBX = coords.MPI_SortByX(MPI)
    MPI_SBX.settings(xboxsize, nxgrid, origin=xorigin)
    MPI_SBX.input(data)
    MPI_SBX.limits4grid()
    data = MPI_SBX.distribute()
    limits = [MPI_SBX.limits[0], MPI_SBX.limits[1], yorigin, yorigin+yboxsize]
    
    # coordinates only inside subbox
    x, y, f, mass = data[:,0], data[:,1], data[:,2], data[:,3]

    # create grids
    xedges, xgrid = shift.cart.mpi_grid1D(xboxsize, nxgrid, MPI, origin=xorigin)
    yedges, ygrid = shift.cart.grid1D(yboxsize, nygrid, origin=yorigin)
    x2d, y2d = shift.cart.mpi_grid2D(boxsize, ngrid, MPI, origin=origin)
    x2d = x2d.flatten()
    y2d = y2d.flatten()

    # local grid, we will use an underscore to differentiate from global
    _xorigin = xedges[0]
    _xboxsize = xedges[-1] - xedges[0]
    _yorigin = yedges[0]
    _yboxsize = yedges[-1] - yedges[0]

    _nxgrid = len(xgrid)
    _nygrid = len(ygrid)

    field, exterior_border, pixID, count = dtfe4grid.delaunay_field4grid2D(
        x, y, f, [_xboxsize, _yboxsize], [_nxgrid, _nygrid], origin=[_xorigin, _yorigin], mass=mass, 
        partition=partition, periodic=[False, False], fbuffer=fbuffer, subsampling=subsampling, 
        outputgrid=False, outputexterior=True, normalise=False, flatten=True
    )

    MPI.wait()

    # get border particles from adjacent slabs with the slab below

    exterior_border_below = {}
    
    exterior_border_below['x'] = MPI.send_up(exterior_border['x'])
    exterior_border_below['y'] = MPI.send_up(exterior_border['y'])
    exterior_border_below['mass'] = MPI.send_up(exterior_border['mass'])
    exterior_border_below['ptype'] = MPI.send_up(exterior_border['ptype'])
    exterior_border_below['f'] = MPI.send_up(exterior_border['f'])
    exterior_border_below['dtfe_mode'] = exterior_border['dtfe_mode']
    exterior_border_below['simplices'] = MPI.send_up(exterior_border['simplices'])
    exterior_border_below['simptypes'] = MPI.send_up(exterior_border['simptypes'])

    # get border particles from adjacent slabs with the slab above
    
    exterior_border_above = {}

    exterior_border_above['x'] = MPI.send_down(exterior_border['x'])
    exterior_border_above['y'] = MPI.send_down(exterior_border['y'])
    exterior_border_above['mass'] = MPI.send_down(exterior_border['mass'])
    exterior_border_above['ptype'] = MPI.send_down(exterior_border['ptype'])
    exterior_border_above['f'] = MPI.send_down(exterior_border['f'])
    exterior_border_above['dtfe_mode'] = exterior_border['dtfe_mode']
    exterior_border_above['simplices'] = MPI.send_down(exterior_border['simplices'])
    exterior_border_above['simptypes'] = MPI.send_down(exterior_border['simptypes'])
    
    # Construct merged Delaunay tesselation

    dtfe = dtfe2d.DelaunayMerger2D()

    if yperiodic:
        jmin, jmax = -1, 2
    else:
        jmin, jmax = 0, 1

    for j in range(jmin, jmax):
        _exterior_border_below = deepcopy(exterior_border_below)
        _exterior_border_below['y'] += j*yboxsize
        if MPI.rank != 0:
            dtfe.add_border(_exterior_border_below)
        else:
            if xperiodic:
                _exterior_border_below['x'] -= xboxsize
                dtfe.add_border(_exterior_border_below)
        
        _exterior_border = deepcopy(exterior_border)
        _exterior_border['y'] += j*yboxsize
        dtfe.add_border(_exterior_border)

        _exterior_border_above = deepcopy(exterior_border_above)
        _exterior_border_above['y'] += j*yboxsize
        if MPI.rank != MPI.size-1:
            dtfe.add_border(_exterior_border_above)
        else:
            if xperiodic:
                _exterior_border_above['x'] += xboxsize
                dtfe.add_border(_exterior_border_above)

    _xbuffer = fbuffer*_xboxsize
    _ybuffer = fbuffer*_yboxsize
    
    boundary = [_xorigin-_xbuffer, _xorigin+_xboxsize+_xbuffer, _yorigin-_ybuffer, _yorigin+_yboxsize+_ybuffer]

    dtfe.run(boundary, apply_filter=True)

    cond = np.where((x2d[pixID] >= boundary[0]) & (x2d[pixID] < boundary[1]) & (y2d[pixID] >= boundary[2]) & (y2d[pixID] < boundary[3]))[0]

    for _dx2d, _dy2d in zip(dx2d, dy2d):
        dtfe_estimate = dtfe.estimate(x2d[pixID[cond]]+_dx2d, y2d[pixID[cond]]+_dy2d)
        cond_isfinite = np.where(np.isfinite(dtfe_estimate) & (count[cond] != len(dx2d)))[0]
        field[pixID[cond[cond_isfinite]]] += dtfe_estimate[cond_isfinite]
        count[cond[cond_isfinite]] += 1.
    
    cond2 = np.where(count[cond] == len(dx2d))[0]
    pixID = np.delete(pixID, cond[cond2])
    count = np.delete(count, cond[cond2])

    exterior_border = dtfe.get_border()

    field /= len(dx2d)

    cond = np.where(count == 0.)[0]
    field[pixID[cond]] *= len(dx2d)
    cond = np.where(count != 0.)[0]
    field[pixID[cond]] *= len(dx2d)/count[cond]

    field = field.reshape(_nxgrid, _nygrid)

    if outputgrid:
        return field, x2d, y2d
    else:
        return field


def mpi_delaunay_density4grid3D(
    x: np.ndarray, y: np.ndarray, z: np.ndarray, boxsize: Union[float, List[float]], ngrid: Union[int, List[int]], 
    MPI: object, origin: Union[float, List[float]] = 0., mass: Optional[np.ndarray] = None, 
    partition: int = 1, periodic : Union[bool, List[bool]] = True, fbuffer: float = 0.5, subsampling: int = 1, 
    outputgrid: bool = False
) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """
    Returns the Delaunay tesselation density on a grid.

    Parameters
    ----------
    x, y, z: array
        Coordinates of particles.
    boxsize : float or list
        Dimensions of the grid.
    ngrid : int or int list
        Grid dimensions.
    MPI : class object
        shift.mpiutils MPI object.
    origin : float or list, optional
        Origin for grid.
    mass : array, optional
        Mass or mass weights for each the particles.
    partition : int or list, optional
        The number of internal Delaunay tesselations used to compute the density grid.
    periodic : bool or list, optional
        Determines whether periodic boundaries are applied.
    fbuffer : float, optional
        Buffer factor length, used to trim exterior buffer region.
    subsampling : int, optional
        The pixel subsampling rate. Each pixel is evaluated subsampling^2 points
        on a grid within each pixel. This is to ensure each pixel is assigned a
        mean pixel value and not the value at the center.
    outputgrid : bool, optional
        Outputs coordinate grid.
    outputexterior : bool, optional
        Ouput exterior information, including exterior border dtfe and current unassigned pixels.

    Returns
    -------
    dens : ndarray
        Density field values on a grid.
    x3d, y3d, z3d : ndarray, optional
        Pixel coordinate points.
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
    dx = xboxsize/nxgrid
    dy = yboxsize/nygrid
    dz = yboxsize/nzgrid

    if np.isscalar(subsampling):
        dnxgrid, dnygrid, dnzgrid = subsampling, subsampling, subsampling
    else:
        dnxgrid, dnygrid, dnzgrid = subsampling[0], subsampling[1], subsampling[2]
    
    # define subsampling points for each pixel
    dx3d, dy3d, dz3d = shift.cart.grid3D([dx, dy, dz], [dnxgrid, dnygrid, dnzgrid], origin=[-dx/2., -dy/2., -dz/2.])
    dx3d, dy3d, dz3d = dx3d.flatten(), dy3d.flatten(), dz3d.flatten()

    # collapse data
    if x is not None:
        if mass is None:
            data = coords.coord2points([x, y, z, np.ones(len(x))])
        else:
            data = coords.coord2points([x, y, z, mass])
    else:
        data = None
    
    # check size of data
    if data is not None:
        lendata = len(data)
    else:
        lendata = 0
    size = MPI.collect(lendata, outlist=True)
    if MPI.rank == 0:
        npart = np.sum(np.array(size))
        MPI.send(npart, tag=11)
    else:
        npart = MPI.recv(0, tag=11)
    MPI.wait()

    # sort coordinates and distribute by coordinate system
    MPI_SBX = coords.MPI_SortByX(MPI)
    MPI_SBX.settings(xboxsize, nxgrid, origin=xorigin)
    MPI_SBX.input(data)
    MPI_SBX.limits4grid()
    data = MPI_SBX.distribute()

    # TODO: This isn't used for anything but hesitant to remove since it's good for diagnostics...
    # limits = [MPI_SBX.limits[0], MPI_SBX.limits[1], yorigin, yorigin+yboxsize, zorigin, zorigin+zboxsize] 
    
    # coordinates only inside subbox
    x, y, z, mass = data[:,0], data[:,1], data[:,2], data[:,3]

    # create grids
    xedges, xgrid = shift.cart.mpi_grid1D(xboxsize, nxgrid, MPI, origin=xorigin)
    yedges, ygrid = shift.cart.grid1D(yboxsize, nygrid, origin=yorigin)
    zedges, zgrid = shift.cart.grid1D(zboxsize, nzgrid, origin=zorigin)
    x3d, y3d, z3d = shift.cart.mpi_grid3D(boxsize, ngrid, MPI, origin=origin)
    x3d = x3d.flatten()
    y3d = y3d.flatten()
    z3d = z3d.flatten()

    # local grid, we will use an underscore to differentiate from global
    _xorigin = xedges[0]
    _xboxsize = xedges[-1] - xedges[0]
    _yorigin = yedges[0]
    _yboxsize = yedges[-1] - yedges[0]
    _zorigin = zedges[0]
    _zboxsize = zedges[-1] - zedges[0]

    _nxgrid = len(xgrid)
    _nygrid = len(ygrid)
    _nzgrid = len(zgrid)

    dens, exterior_border, pixID, count = dtfe4grid.delaunay_density4grid3D(
        x, y, z, [_xboxsize, _yboxsize, _zboxsize], [_nxgrid, _nygrid, _nzgrid], 
        origin=[_xorigin, _yorigin, _zorigin], mass=mass, partition=partition, 
        periodic=[False, False, False], fbuffer=fbuffer, subsampling=subsampling, 
        outputgrid=False, outputexterior=True, normalise=False, flatten=True
    )

    MPI.wait()

    # get border particles from adjacent slabs with the slab below

    exterior_border_below = {}
    
    exterior_border_below['x'] = MPI.send_up(exterior_border['x'])
    exterior_border_below['y'] = MPI.send_up(exterior_border['y'])
    exterior_border_below['z'] = MPI.send_up(exterior_border['z'])
    exterior_border_below['mass'] = MPI.send_up(exterior_border['mass'])
    exterior_border_below['ptype'] = MPI.send_up(exterior_border['ptype'])
    exterior_border_below['f'] = MPI.send_up(exterior_border['f'])
    exterior_border_below['dtfe_mode'] = exterior_border['dtfe_mode']
    exterior_border_below['simplices'] = MPI.send_up(exterior_border['simplices'])
    exterior_border_below['simptypes'] = MPI.send_up(exterior_border['simptypes'])

    # get border particles from adjacent slabs with the slab above
    
    exterior_border_above = {}

    exterior_border_above['x'] = MPI.send_down(exterior_border['x'])
    exterior_border_above['y'] = MPI.send_down(exterior_border['y'])
    exterior_border_above['z'] = MPI.send_down(exterior_border['z'])
    exterior_border_above['mass'] = MPI.send_down(exterior_border['mass'])
    exterior_border_above['ptype'] = MPI.send_down(exterior_border['ptype'])
    exterior_border_above['f'] = MPI.send_down(exterior_border['f'])
    exterior_border_above['dtfe_mode'] = exterior_border['dtfe_mode']
    exterior_border_above['simplices'] = MPI.send_down(exterior_border['simplices'])
    exterior_border_above['simptypes'] = MPI.send_down(exterior_border['simptypes'])
    
    # Construct merged Delaunay tesselation

    dtfe = dtfe3d.DelaunayMerger3D()

    if yperiodic:
        jmin, jmax = -1, 2
    else:
        jmin, jmax = 0, 1

    if zperiodic:
        kmin, kmax = -1, 2
    else:
        kmin, kmax = 0, 1

    for j in range(jmin, jmax):
        for k in range(kmin, kmax):
            _exterior_border_below = deepcopy(exterior_border_below)
            _exterior_border_below['y'] += j*yboxsize
            _exterior_border_below['z'] += k*zboxsize
            if MPI.rank != 0:
                dtfe.add_border(_exterior_border_below)
            else:
                if xperiodic:
                    _exterior_border_below['x'] -= xboxsize
                    dtfe.add_border(_exterior_border_below)
            
            _exterior_border = deepcopy(exterior_border)
            _exterior_border['y'] += j*yboxsize
            _exterior_border['z'] += k*zboxsize
            dtfe.add_border(_exterior_border)

            _exterior_border_above = deepcopy(exterior_border_above)
            _exterior_border_above['y'] += j*yboxsize
            _exterior_border_above['z'] += k*zboxsize
            if MPI.rank != MPI.size-1:
                dtfe.add_border(_exterior_border_above)
            else:
                if xperiodic:
                    _exterior_border_above['x'] += xboxsize
                    dtfe.add_border(_exterior_border_above)

    _xbuffer = fbuffer*_xboxsize
    _ybuffer = fbuffer*_yboxsize
    _zbuffer = fbuffer*_zboxsize
    
    boundary = [_xorigin-_xbuffer, _xorigin+_xboxsize+_xbuffer, _yorigin-_ybuffer, _yorigin+_yboxsize+_ybuffer, _zorigin-_zbuffer, _zorigin+_zboxsize+_zbuffer]

    dtfe.run(boundary, apply_filter=True)

    cond = np.where(
        (x3d[pixID] >= boundary[0]) & (x3d[pixID] < boundary[1]) & 
        (y3d[pixID] >= boundary[2]) & (y3d[pixID] < boundary[3]) & 
        (z3d[pixID] >= boundary[4]) & (z3d[pixID] < boundary[5]))[0]

    for _dx3d, _dy3d, _dz3d in zip(dx3d, dy3d, dz3d):
        dtfe_estimate = dtfe.estimate(x3d[pixID[cond]]+_dx3d, y3d[pixID[cond]]+_dy3d, z3d[pixID[cond]]+_dz3d)
        cond_isfinite = np.where(np.isfinite(dtfe_estimate) & (count[cond] != len(dx3d)))[0]
        dens[pixID[cond[cond_isfinite]]] += dtfe_estimate[cond_isfinite]
        count[cond[cond_isfinite]] += 1.

    cond2 = np.where(count[cond] == len(dx3d))[0]
    pixID = np.delete(pixID, cond[cond2])
    count = np.delete(count, cond[cond2])

    exterior_border = dtfe.get_border()

    dens /= len(dx3d)

    cond = np.where(count == 0.)[0]
    dens[pixID[cond]] *= len(dx3d)
    cond = np.where(count != 0.)[0]
    dens[pixID[cond]] *= len(dx3d)/count[cond]

    dens = dens.reshape(_nxgrid, _nygrid, _nzgrid)

    if outputgrid:
        return dens, x3d, y3d, z3d
    else:
        return dens


def mpi_delaunay_field4grid3D(
    x: np.ndarray, y: np.ndarray, z: np.ndarray, f: np.ndarray, boxsize: Union[float, List[float]], 
    ngrid: Union[int, List[int]], MPI: object, origin: Union[float, List[float]] = 0., mass: Optional[np.ndarray] = None, 
    partition: int = 1, periodic : Union[bool, List[bool]] = True, fbuffer: float = 0.5, subsampling: int = 1, 
    outputgrid: bool = False
) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """
    Returns the Delaunay tesselation field a grid.

    Parameters
    ----------
    x, y, z: array
        Coordinates of particles.
    f : array
        Field values to be interpolated.
    boxsize : float or list
        Dimensions of the grid.
    ngrid : int or int list
        Grid dimensions.
    MPI : class object
        shift.mpiutils MPI object.
    origin : float or list, optional
        Origin for grid.
    mass : array, optional
        Mass or mass weights for each the particles.
    partition : int or list, optional
        The number of internal Delaunay tesselations used to compute the field grid.
    periodic : bool or list, optional
        Determines whether periodic boundaries are applied.
    fbuffer : float, optional
        Buffer factor length, used to trim exterior buffer region.
    subsampling : int, optional
        The pixel subsampling rate. Each pixel is evaluated subsampling^2 points
        on a grid within each pixel. This is to ensure each pixel is assigned a
        mean pixel value and not the value at the center.
    outputgrid : bool, optional
        Outputs coordinate grid.
    outputexterior : bool, optional
        Ouput exterior information, including exterior border dtfe and current unassigned pixels.

    Returns
    -------
    field : ndarray
        Density field values on a grid.
    x3d, y3d, z3d : ndarray, optional
        Pixel coordinate points.
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
    dx = xboxsize/nxgrid
    dy = yboxsize/nygrid
    dz = yboxsize/nzgrid

    if np.isscalar(subsampling):
        dnxgrid, dnygrid, dnzgrid = subsampling, subsampling, subsampling
    else:
        dnxgrid, dnygrid, dnzgrid = subsampling[0], subsampling[1], subsampling[2]
    
    # define subsampling points for each pixel
    dx3d, dy3d, dz3d = shift.cart.grid3D([dx, dy, dz], [dnxgrid, dnygrid, dnzgrid], origin=[-dx/2., -dy/2., -dz/2.])
    dx3d, dy3d, dz3d = dx3d.flatten(), dy3d.flatten(), dz3d.flatten()

    # collapse data
    if x is not None:
        if mass is None:
            data = coords.coord2points([x, y, z, f, np.ones(len(x))])
        else:
            data = coords.coord2points([x, y, z, f, mass])
    else:
        data = None
    
    # check size of data
    if data is not None:
        lendata = len(data)
    else:
        lendata = 0
    size = MPI.collect(lendata, outlist=True)
    if MPI.rank == 0:
        npart = np.sum(np.array(size))
        MPI.send(npart, tag=11)
    else:
        npart = MPI.recv(0, tag=11)
    MPI.wait()

    # sort coordinates and distribute by coordinate system
    MPI_SBX = coords.MPI_SortByX(MPI)
    MPI_SBX.settings(xboxsize, nxgrid, origin=xorigin)
    MPI_SBX.input(data)
    MPI_SBX.limits4grid()
    data = MPI_SBX.distribute()

    # TODO: This isn't used for anything but hesitant to remove since it's good for diagnostics...
    # limits = [MPI_SBX.limits[0], MPI_SBX.limits[1], yorigin, yorigin+yboxsize, zorigin, zorigin+zboxsize] 
    
    # coordinates only inside subbox
    x, y, z, f, mass = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4]

    # create grids
    xedges, xgrid = shift.cart.mpi_grid1D(xboxsize, nxgrid, MPI, origin=xorigin)
    yedges, ygrid = shift.cart.grid1D(yboxsize, nygrid, origin=yorigin)
    zedges, zgrid = shift.cart.grid1D(zboxsize, nzgrid, origin=zorigin)
    x3d, y3d, z3d = shift.cart.mpi_grid3D(boxsize, ngrid, MPI, origin=origin)
    x3d = x3d.flatten()
    y3d = y3d.flatten()
    z3d = z3d.flatten()

    # local grid, we will use an underscore to differentiate from global
    _xorigin = xedges[0]
    _xboxsize = xedges[-1] - xedges[0]
    _yorigin = yedges[0]
    _yboxsize = yedges[-1] - yedges[0]
    _zorigin = zedges[0]
    _zboxsize = zedges[-1] - zedges[0]

    _nxgrid = len(xgrid)
    _nygrid = len(ygrid)
    _nzgrid = len(zgrid)

    field, exterior_border, pixID, count = dtfe4grid.delaunay_field4grid3D(
        x, y, z, f, [_xboxsize, _yboxsize, _zboxsize], [_nxgrid, _nygrid, _nzgrid], 
        origin=[_xorigin, _yorigin, _zorigin], mass=mass, partition=partition, 
        periodic=[False, False, False], fbuffer=fbuffer, subsampling=subsampling, 
        outputgrid=False, outputexterior=True, normalise=False, flatten=True
    )

    MPI.wait()

    # get border particles from adjacent slabs with the slab below

    exterior_border_below = {}
    
    exterior_border_below['x'] = MPI.send_up(exterior_border['x'])
    exterior_border_below['y'] = MPI.send_up(exterior_border['y'])
    exterior_border_below['z'] = MPI.send_up(exterior_border['z'])
    exterior_border_below['mass'] = MPI.send_up(exterior_border['mass'])
    exterior_border_below['ptype'] = MPI.send_up(exterior_border['ptype'])
    exterior_border_below['f'] = MPI.send_up(exterior_border['f'])
    exterior_border_below['dtfe_mode'] = exterior_border['dtfe_mode']
    exterior_border_below['simplices'] = MPI.send_up(exterior_border['simplices'])
    exterior_border_below['simptypes'] = MPI.send_up(exterior_border['simptypes'])

    # get border particles from adjacent slabs with the slab above
    
    exterior_border_above = {}

    exterior_border_above['x'] = MPI.send_down(exterior_border['x'])
    exterior_border_above['y'] = MPI.send_down(exterior_border['y'])
    exterior_border_above['z'] = MPI.send_down(exterior_border['z'])
    exterior_border_above['mass'] = MPI.send_down(exterior_border['mass'])
    exterior_border_above['ptype'] = MPI.send_down(exterior_border['ptype'])
    exterior_border_above['f'] = MPI.send_down(exterior_border['f'])
    exterior_border_above['dtfe_mode'] = exterior_border['dtfe_mode']
    exterior_border_above['simplices'] = MPI.send_down(exterior_border['simplices'])
    exterior_border_above['simptypes'] = MPI.send_down(exterior_border['simptypes'])
    
    # Construct merged Delaunay tesselation

    dtfe = dtfe3d.DelaunayMerger3D()

    if yperiodic:
        jmin, jmax = -1, 2
    else:
        jmin, jmax = 0, 1

    if zperiodic:
        kmin, kmax = -1, 2
    else:
        kmin, kmax = 0, 1

    for j in range(jmin, jmax):
        for k in range(kmin, kmax):
            _exterior_border_below = deepcopy(exterior_border_below)
            _exterior_border_below['y'] += j*yboxsize
            _exterior_border_below['z'] += k*zboxsize
            if MPI.rank != 0:
                dtfe.add_border(_exterior_border_below)
            else:
                if xperiodic:
                    _exterior_border_below['x'] -= xboxsize
                    dtfe.add_border(_exterior_border_below)
            
            _exterior_border = deepcopy(exterior_border)
            _exterior_border['y'] += j*yboxsize
            _exterior_border['z'] += k*zboxsize
            dtfe.add_border(_exterior_border)

            _exterior_border_above = deepcopy(exterior_border_above)
            _exterior_border_above['y'] += j*yboxsize
            _exterior_border_above['z'] += k*zboxsize
            if MPI.rank != MPI.size-1:
                dtfe.add_border(_exterior_border_above)
            else:
                if xperiodic:
                    _exterior_border_above['x'] += xboxsize
                    dtfe.add_border(_exterior_border_above)

    _xbuffer = fbuffer*_xboxsize
    _ybuffer = fbuffer*_yboxsize
    _zbuffer = fbuffer*_zboxsize
    
    boundary = [_xorigin-_xbuffer, _xorigin+_xboxsize+_xbuffer, _yorigin-_ybuffer, _yorigin+_yboxsize+_ybuffer, _zorigin-_zbuffer, _zorigin+_zboxsize+_zbuffer]

    dtfe.run(boundary, apply_filter=True)

    cond = np.where(
        (x3d[pixID] >= boundary[0]) & (x3d[pixID] < boundary[1]) & 
        (y3d[pixID] >= boundary[2]) & (y3d[pixID] < boundary[3]) & 
        (z3d[pixID] >= boundary[4]) & (z3d[pixID] < boundary[5]))[0]

    for _dx3d, _dy3d, _dz3d in zip(dx3d, dy3d, dz3d):
        dtfe_estimate = dtfe.estimate(x3d[pixID[cond]]+_dx3d, y3d[pixID[cond]]+_dy3d, z3d[pixID[cond]]+_dz3d)
        cond_isfinite = np.where(np.isfinite(dtfe_estimate) & (count[cond] != len(dx3d)))[0]
        field[pixID[cond[cond_isfinite]]] += dtfe_estimate[cond_isfinite]
        count[cond[cond_isfinite]] += 1.

    cond2 = np.where(count[cond] == len(dx3d))[0]
    pixID = np.delete(pixID, cond[cond2])
    count = np.delete(count, cond[cond2])

    exterior_border = dtfe.get_border()

    field /= len(dx3d)

    cond = np.where(count == 0.)[0]
    field[pixID[cond]] *= len(dx3d)
    cond = np.where(count != 0.)[0]
    field[pixID[cond]] *= len(dx3d)/count[cond]

    field = field.reshape(_nxgrid, _nygrid, _nzgrid)

    if outputgrid:
        return field, x3d, y3d, z3d
    else:
        return field

