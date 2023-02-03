import numpy as np
import matplotlib.pylab as plt

import shift

from .. import coords
from .. import boundary

from . import dtfe4grid


def mpi_dtfe4grid2D(x, y, ngrid, boxsize, MPI, MPI_split, f=None,
                    buffer_type=None, buffer_length=0., buffer_val=0.,
                    origin=0., subsampling=4, outputgrid=False):
    """Returns the Delaunay tesselation density or field on a grid.

    Parameters
    ----------
    x, y : array
        Coordinates of particles.
    ngrid : int or int list
        Grid dimensions.
    boxsize : float or list
        Dimensions of the grid.
    MPI : class object
        MPIutils MPI class object.
    MPI_split : int or int list
        Determines how to split each axis for serial DTFE calculations.
    f : array, optional
        Field values, if None assumed output is density.
    buffer_type : str, optional
        Buffer particle type, either:
            - 'random' for random buffer particles.
            - 'periodic' for periodic buffer particles.
            - None for no buffer particles.
    buffer_length : float, optional
        Buffer length.
    buffer_val : float, optional
        Value given to random buffer particles.
    subsampling : int, optional
        The pixel subsampling rate. Each pixel is evaluated subsampling^2 points
        on a grid within each pixel. This is to ensure each pixel is assigned a
        mean pixel value and not the value at the center.
    outputgrid : bool, optional
        Outputs coordinate grid.

    Returns
    -------
    f2D : ndarray
        Field values on a grid.
    x2D, y2D : ndarray, optional
        Pixel coordinate points.
    """
    # check basic inputs
    if np.isscalar(ngrid):
        nxgrid, nygrid = ngrid, ngrid
    else:
        nxgrid, nygrid = ngrid[0], ngrid[1]
    if np.isscalar(boxsize):
        xboxsize, yboxsize = boxsize, boxsize
    else:
        xboxsize, yboxsize = boxsize[0], boxsize[1]
    if np.isscalar(origin):
        xorigin, yorigin = origin, origin
    else:
        xorigin, yorigin = origin[0], origin[1]
    # check buffer length
    if buffer_type is not None:
        buffer_length = buffer_length
    else:
        buffer_length = 0.
    # collapse data
    if x is not None:
        if f is None:
            data = coords.coord2points([x, y])
        else:
            data = coords.coord2points([x, y, f])
    else:
        data = None
    # check size of data and calculate densities
    size = MPI.collect(len(data), outlist=True)
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
    # add buffer particles.
    if buffer_type == 'random':
        _data = boundary.mpi_buffer_random_utils(data, limits, buffer_length, MPI)
        xr, yr = boundary.mpi_buffer_random_2D(npart, boxsize, limits, buffer_length, MPI)
        _x, _y = _data[:,0], _data[:,1]
        _x = np.concatenate([_x, xr])
        _y = np.concatenate([_y, yr])
        if len(_data) == 3:
            _f = np.concatenate([_data[:,2], np.ones(npart)*buffer_val])
    elif buffer_type == 'periodic':
        _data = boundary.mpi_buffer_periodic_2D(data, xboxsize, limits, buffer_length, MPI)
        _x, _y = _data[:,0], _data[:,1]
        if len(data) == 3:
            _f = _data[:,2]
    else:
        _x, _y = data[:,0], data[:,1]
        if len(data) == 3:
            _f = data[:,2]
    # create grids
    xedges, xgrid = shift.cart.mpi_grid1D(xboxsize, nxgrid, MPI, origin=xorigin)
    yedges, ygrid = shift.cart.grid1D(yboxsize, nygrid, origin=yorigin)
    x2D, y2D = shift.cart.mpi_grid2D(boxsize, ngrid, MPI, origin=origin)
    shape = np.shape(x2D)
    nxgrid, nygrid = shape[0], shape[1]
    f2D = np.zeros((nxgrid, nygrid))
    if np.isscalar(MPI_split):
        xs1, xs2 = 0, len(xgrid)
        ys1, ys2 = MPI.split(len(ygrid), size=MPI_split)
    else:
        xs1, xs2 = MPI.split(len(xgrid), size=MPI_split[0])
        ys1, ys2 = MPI.split(len(ygrid), size=MPI_split[1])
    xs1, ys1 = np.meshgrid(xs1, ys1, indexing='ij')
    xs1, ys1 = xs1.flatten(), ys1.flatten()
    xs2, ys2 = np.meshgrid(xs2, ys2, indexing='ij')
    xs2, ys2 = xs2.flatten(), ys2.flatten()
    # calculate dtfe
    for i in range(0, len(xs1)):
        _xmin, _xmax = xedges[xs1[i]], xedges[xs2[i]]
        _ymin, _ymax = yedges[ys1[i]], yedges[ys2[i]]
        _nxgrid, _nygrid = xs2[i]-xs1[i], ys2[i]-ys1[i]
        cond = np.where((_x >= _xmin-buffer_length) &
                        (_x < _xmax+buffer_length) &
                        (_y >= _ymin-buffer_length) &
                        (_y < _ymax+buffer_length))[0]
        _xx, _yy = _x[cond], _y[cond]
        if f is None:
            _ff = None
        else:
            _ff = _f[cond]
        _f2d = dtfe4grid.dtfe4grid2D(_xx, _yy, [_nxgrid, _nygrid], [_xmax-_xmin, _ymax-_ymin],
                                     f=_ff, origin=[_xmin, _ymin], buffer_type=None,
                                     subsampling=subsampling, outputgrid=False)
        f2D[xs1[i]:xs2[i],ys1[i]:ys2[i]] = _f2d
    if outputgrid:
        return x2D, y2D, f2D
    else:
        return f2D


# def mpi_dtfe4grid3D(x, y, z, ngrid, boxsize, MPI, MPI_split, f=None, origin=0.,
#                     buffer_type=None, buffer_length=0., buffer_val=0., subsampling=4,
#                     outputgrid=False):
#     """Returns the Delaunay tesselation density or field on a grid.
#
#     Parameters
#     ----------
#     x, y, z : array
#         Coordinates of particles.
#     ngrid : int or int list
#         Grid dimensions.
#     boxsize : float or list
#         Dimensions of the grid.
#     MPI : class object
#         MPIutils MPI class object.
#     MPI_split : int or int list
#         Determines how to split each axis for serial DTFE calculations.
#     f : array, optional
#         Field values, if None assumed output is density.
#     origin : float or list, optional
#         Origin for grid.
#     buffer_type : str, optional
#         Buffer particle type, either:
#             - 'random' for random buffer particles.
#             - 'periodic' for periodic buffer particles.
#             - None for no buffer particles.
#     buffer_factor : float, optional
#         Buffer length given as a multiple of the interparticle separation.
#     buffer_val : float, optional
#         Value given to random buffer particles.
#     subsampling : int, optional
#         The pixel subsampling rate. Each pixel is evaluated subsampling^2 points
#         on a grid within each pixel. This is to ensure each pixel is assigned a
#         mean pixel value and not the value at the center.
#     outputgrid : bool, optional
#         Outputs coordinate grid.
#
#     Returns
#     -------
#     f3D : ndarray
#         Field values on a grid.
#     x3D, y3D, z3D : ndarray, optional
#         Pixel coordinate points.
#     """
#     if np.isscalar(ngrid):
#         nxgrid, nygrid, nzgrid = ngrid, ngrid, ngrid
#     else:
#         nxgrid, nygrid, nzgrid = ngrid[0], ngrid[1], ngrid[2]
#     if np.isscalar(boxsize):
#         xboxsize, yboxsize, zboxsize = boxsize, boxsize, boxsize
#     else:
#         xboxsize, yboxsize, zboxsize = boxsize[0], boxsize[1], boxsize[2]
#     if np.isscalar(boxsize):
#         xmin, ymin, zmin = origin, origin, origin
#     else:
#         xmin, ymin, zmin = origin[0], origin[1], origin[2]
#     if MPI.rank == 0:
#         if buffer_type is not None:
#             buffer_length = buffer_factor*dtfe4grid.mean_separation_3D(len(x), boxsize)
#         else:
#             buffer_length = 0.
#         MPI.send(buffer_length, tag=11)
#     else:
#         buffer_length = MPI.recv(0, tag=11)
#     MPI.wait()
#     # if MPI.rank == 0:
#     #     if f is None:
#     #         data = coords.coord2points([x, y, z])
#     #     else:
#     #         data = coords.coord2points([x, y, z, f])
#     #     if buffer_type == 'random':
#     #         xb, yb, zb = boundary.buffer_random_3D(len(x), boxsize, buffer_length, origin=origin)
#     #         _x = np.concatenate([data[:,0], xb])
#     #         _y = np.concatenate([data[:,1], yb])
#     #         _z = np.concatenate([data[:,2], zb])
#     #         if f is not None:
#     #             _f = np.concatenate([data[:,3], buffer_val*np.ones(len(xb))])
#     #     elif buffer_type == 'periodic':
#     #         if f is None:
#     #             data = coords.coord2points([x, y, z])
#     #         else:
#     #             data = coords.coord2points([x, y, z, f])
#     #         datab = boundary.buffer_periodic_3D(data, boxsize, buffer_length, origin=origin)
#     #         _data = np.vstack([data, datab])
#     #         _x, _y, _z = _data[:,0], _data[:,1], _data[:,2]
#     #         if f is not None:
#     #             _f = _data[:,3]
#     #     else:
#     #         _x, _y, _z = x, y, z
#     #         if f is not None:
#     #             _f = f
#
#     #     xedges, xgrid = shift.cart.grid1D(xboxsize, nxgrid, origin=ymin)
#     #     xs1, xs2 = MPI.split(len(xgrid))
#     #     for i in range(1, MPI.size):
#     #         _xmin = xedges[xs1[i]]-buffer_length
#     #         _xmax = xedges[xs2[i]]+buffer_length
#     #         cond = np.where((_x >= _xmin) & (_x < _xmax))[0]
#     #         if f is None:
#     #             _data = coords.coord2points([_x[cond], _y[cond], _z[cond]])
#     #         else:
#     #             _data = coords.coord2points([_x[cond], _y[cond], _z[cond], _f[cond]])
#     #         MPI.send(_data, to_rank=i, tag=10+i)
#     #     _xmin = xedges[xs1[0]]-buffer_length
#     #     _xmax = xedges[xs2[0]]+buffer_length
#     #     cond = np.where((_x >= _xmin) & (_x < _xmax))[0]
#     #     _x, _y, _z = _x[cond], _y[cond], _z[cond]
#     #     if f is not None:
#     #         _f = _f[cond]
#     # else:
#     #     _data = MPI.recv(0, tag=10+MPI.rank)
#     #     _x, _y, _z = _data[:,0], _data[:,1], _data[:,2]
#     #     if len(_data[0]) == 4:
#     #         _f = _data[:,3]
#     if f is None:
#         _f = None
#     MPI.wait()
#     dx = xboxsize/nxgrid
#     dy = yboxsize/nygrid
#     dz = zboxsize/nzgrid
#     xedges, xgrid = shift.cart.mpi_grid1D(xboxsize, nxgrid, MPI, origin=xmin)
#     yedges, ygrid = shift.cart.grid1D(yboxsize, nygrid, origin=ymin)
#     zedges, zgrid = shift.cart.grid1D(zboxsize, nzgrid, origin=zmin)
#     x3D, y3D, z3D = shift.cart.mpi_grid3D(boxsize, ngrid, MPI, origin=origin)
#     shape = np.shape(x3D)
#     nxgrid, nygrid, nzgrid = shape[0], shape[1], shape[2]
#     f3D = np.zeros((nxgrid, nygrid, nzgrid))
#     if np.isscalar(MPI_split):
#         xs1, xs2 = 0, len(xgrid)
#         ys1, ys2 = MPI.split(len(ygrid), size=MPI_split)
#         zs1, zs2 = MPI.split(len(zgrid), size=MPI_split)
#     else:
#         xs1, xs2 = MPI.split(len(xgrid), size=MPI_split[0])
#         ys1, ys2 = MPI.split(len(ygrid), size=MPI_split[1])
#         zs1, zs2 = MPI.split(len(zgrid), size=MPI_split[2])
#     xs1, ys1, zs1 = np.meshgrid(xs1, ys1, zs1, indexing='ij')
#     xs1, ys1, zs1 = xs1.flatten(), ys1.flatten(), zs1.flatten()
#     xs2, ys2, zs2 = np.meshgrid(xs2, ys2, zs2, indexing='ij')
#     xs2, ys2, zs2 = xs2.flatten(), ys2.flatten(), zs2.flatten()
#     for i in range(0, len(xs1)):
#         _xmin, _xmax = xedges[xs1[i]], xedges[xs2[i]]
#         _ymin, _ymax = yedges[ys1[i]], yedges[ys2[i]]
#         _zmin, _zmax = zedges[zs1[i]], zedges[zs2[i]]
#         _nxgrid, _nygrid, _nzgrid = xs2[i]-xs1[i], ys2[i]-ys1[i], zs2[i]-zs1[i]
#         cond = np.where((_x >= _xmin-buffer_length) &
#                         (_x < _xmax+buffer_length) &
#                         (_y >= _ymin-buffer_length) &
#                         (_y < _ymax+buffer_length) &
#                         (_z >= _zmin-buffer_length) &
#                         (_z < _zmax+buffer_length))[0]
#         _xx, _yy, _zz = _x[cond], _y[cond], _z[cond]
#         if f is None:
#             _ff = None
#         else:
#             _ff = _f[cond]
#         if len(_xx) > 0:
#             _f3d = dtfe4grid.dtfe4grid3D(_xx, _yy, _zz, [_nxgrid, _nygrid, _nzgrid],
#                                          [_xmax-_xmin, _ymax-_ymin, _zmax-_zmin],
#                                          f=_ff, origin=[_xmin, _ymin, _zmin],
#                                          buffer_type=None, subsampling=subsampling,
#                                          outputgrid=False)
#             f3D[xs1[i]:xs2[i],ys1[i]:ys2[i],zs1[i]:zs2[i]] = _f3d
#     if outputgrid:
#         return x3D, y3D, f3D
#     else:
#         return f3D


def mpi_dtfe4grid3D(x, y, z, ngrid, boxsize, MPI, MPI_split, f=None,
                    buffer_type=None, buffer_length=0., buffer_val=0.,
                    origin=0., subsampling=4, outputgrid=False):
    """Returns the Delaunay tesselation density or field on a grid.

    Parameters
    ----------
    x, y, z : array
        Coordinates of particles.
    ngrid : int or int list
        Grid dimensions.
    boxsize : float or list
        Dimensions of the grid.
    MPI : class object
        MPIutils MPI class object.
    MPI_split : int or int list
        Determines how to split each axis for serial DTFE calculations.
    f : array, optional
        Field values, if None assumed output is density.
    buffer_type : str, optional
        Buffer particle type, either:
            - 'random' for random buffer particles.
            - 'periodic' for periodic buffer particles.
            - None for no buffer particles.
    buffer_length : float, optional
        Buffer length.
    buffer_val : float, optional
        Value given to random buffer particles.
    subsampling : int, optional
        The pixel subsampling rate. Each pixel is evaluated subsampling^2 points
        on a grid within each pixel. This is to ensure each pixel is assigned a
        mean pixel value and not the value at the center.
    outputgrid : bool, optional
        Outputs coordinate grid.

    Returns
    -------
    f3D : ndarray
        Field values on a grid.
    x3D, y3D, z3D : ndarray, optional
        Pixel coordinate points.
    """
    # check basic inputs
    if np.isscalar(ngrid):
        nxgrid, nygrid, nzgrid = ngrid, ngrid, ngrid
    else:
        nxgrid, nygrid, nzgrid = ngrid[0], ngrid[1], ngrid[2]
    if np.isscalar(boxsize):
        xboxsize, yboxsize, zboxsize = boxsize, boxsize, boxsize
    else:
        xboxsize, yboxsize, zboxsize = boxsize[0], boxsize[1], boxsize[2]
    if np.isscalar(origin):
        xorigin, yorigin, zorigin = origin, origin, origin
    else:
        xorigin, yorigin, zorigin = origin[0], origin[1], origin[2]
    # check buffer length
    if buffer_type is not None:
        buffer_length = buffer_length
    else:
        buffer_length = 0.
    # collapse data
    if x is not None:
        if f is None:
            data = coords.coord2points([x, y, z])
        else:
            data = coords.coord2points([x, y, z, f])
    else:
        data = None
    # check size of data and calculate densities
    size = MPI.collect(len(data), outlist=True)
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
    limits = [MPI_SBX.limits[0], MPI_SBX.limits[1], yorigin, yorigin+yboxsize, zorigin, zorigin+zboxsize]
    # add buffer particles.
    if buffer_type == 'random':
        _data = boundary.mpi_buffer_random_utils(data, limits, buffer_length, MPI)
        xr, yr, zr = boundary.mpi_buffer_random_3D(npart, boxsize, limits, buffer_length, MPI)
        _x, _y, _z = _data[:,0], _data[:,1], _data[:,2]
        _x = np.concatenate([_x, xr])
        _y = np.concatenate([_y, yr])
        _z = np.concatenate([_z, zr])
        if len(_data) == 4:
            _f = np.concatenate([_data[:,3], np.ones(len(xr))*buffer_val])
    elif buffer_type == 'periodic':
        _data = boundary.mpi_buffer_periodic_3D(data, xboxsize, limits, buffer_length, MPI)
        _x, _y, _z = _data[:,0], _data[:,1], _data[:,2]
        if len(data) == 4:
            _f = _data[:,3]
    else:
        _x, _y, _z = data[:,0], data[:,1], data[:,2]
        if len(data) == 4:
            _f = data[:,3]
    # create grids
    xedges, xgrid = shift.cart.mpi_grid1D(xboxsize, nxgrid, MPI, origin=xorigin)
    yedges, ygrid = shift.cart.grid1D(yboxsize, nygrid, origin=yorigin)
    zedges, zgrid = shift.cart.grid1D(zboxsize, nzgrid, origin=zorigin)
    x3D, y3D, z3D = shift.cart.mpi_grid3D(boxsize, ngrid, MPI, origin=origin)
    shape = np.shape(x3D)
    nxgrid, nygrid, nzgrid = shape[0], shape[1], shape[2]
    f3D = np.zeros((nxgrid, nygrid, nzgrid))
    if np.isscalar(MPI_split):
        xs1, xs2 = 0, len(xgrid)
        ys1, ys2 = MPI.split(len(ygrid), size=MPI_split)
        zs1, zs2 = MPI.split(len(zgrid), size=MPI_split)
    else:
        xs1, xs2 = MPI.split(len(xgrid), size=MPI_split[0])
        ys1, ys2 = MPI.split(len(ygrid), size=MPI_split[1])
        zs1, zs2 = MPI.split(len(zgrid), size=MPI_split[2])
    xs1, ys1, zs1 = np.meshgrid(xs1, ys1, zs1, indexing='ij')
    xs1, ys1, zs1 = xs1.flatten(), ys1.flatten(), zs1.flatten()
    xs2, ys2, zs2 = np.meshgrid(xs2, ys2, zs2, indexing='ij')
    xs2, ys2, zs2 = xs2.flatten(), ys2.flatten(), zs2.flatten()
    # calculate dtfe
    for i in range(0, len(xs1)):
        _xmin, _xmax = xedges[xs1[i]], xedges[xs2[i]]
        _ymin, _ymax = yedges[ys1[i]], yedges[ys2[i]]
        _zmin, _zmax = zedges[zs1[i]], zedges[zs2[i]]
        _nxgrid, _nygrid, _nzgrid = xs2[i]-xs1[i], ys2[i]-ys1[i], zs2[i]-zs1[i]
        cond = np.where((_x >= _xmin-buffer_length) &
                        (_x < _xmax+buffer_length) &
                        (_y >= _ymin-buffer_length) &
                        (_y < _ymax+buffer_length) &
                        (_z >= _zmin-buffer_length) &
                        (_z < _zmax+buffer_length))[0]
        _xx, _yy, _zz = _x[cond], _y[cond], _z[cond]
        if f is None:
            _ff = None
        else:
            _ff = _f[cond]
        _f3d = dtfe4grid.dtfe4grid3D(_xx, _yy, _zz, [_nxgrid, _nygrid, _nzgrid],
                                     [_xmax-_xmin, _ymax-_ymin, _zmax-_zmin],
                                     f=_ff, origin=[_xmin, _ymin, _zmin],
                                     buffer_type=None, subsampling=subsampling,
                                     outputgrid=False)
        f3D[xs1[i]:xs2[i],ys1[i]:ys2[i],zs1[i]:zs2[i]] = _f3d
    if outputgrid:
        return x3D, y3D, z3D, f3D
    else:
        return f3D
