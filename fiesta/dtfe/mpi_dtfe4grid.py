import numpy as np

import shift

from .. import coords
from .. import boundary

from . import dtfe4grid


def mpi_dtfe4grid2D(x, y, ngrid, boxsize, MPI, MPI_ngrid, f=None, origin=0.,
                    buffer_type=None, buffer_factor=2, buffer_val=0., subsampling=4,
                    outputgrid=False):
    if np.isscalar(ngrid):
        nxgrid, nygrid = ngrid, ngrid
    else:
        nxgrid, nygrid = ngrid[0], ngrid[1]
    if np.isscalar(boxsize):
        xboxsize, yboxsize = boxsize, boxsize
    else:
        xboxsize, yboxsize = boxsize[0], boxsize[1]
    if np.isscalar(boxsize):
        xmin, ymin = origin, origin
    else:
        xmin, ymin = origin[0], origin[1]
    if MPI.rank == 0:
        if buffer_type is not None:
            buffer_length = buffer_factor*dtfe4grid.mean_separation_2D(len(x), boxsize)
        else:
            buffer_length = 0.
        MPI.send(buffer_length, tag=11)
    else:
        buffer_length = MPI.recv(0, tag=11)
    MPI.wait()
    if MPI.rank == 0:
        if f is None:
            data = coords.coord2points([x, y])
        else:
            data = coords.coord2points([x, y, f])
        if buffer_type == 'random':
            xb, yb = boundary.buffer_random_2D(len(x), boxsize, buffer_length, origin=origin)
            _x = np.concatenate([data[:,0], xb])
            _y = np.concatenate([data[:,1], yb])
            if f is not None:
                _f = np.concatenate([data[:,2], buffer_val*np.ones(len(xb))])
        elif buffer_type == 'periodic':
            if f is None:
                data = coords.coord2points([x, y])
            else:
                data = coords.coord2points([x, y, f])
            datab = boundary.buffer_periodic_2D(data, boxsize, buffer_length, origin=origin)
            _data = np.vstack([data, datab])
            _x, _y = _data[:,0], _data[:,1]
            if f is not None:
                _f = _data[:,2]
        else:
            _x, _y = x, y
            if f is not None:
                _f = f
        xedges, xgrid = shift.cart.grid1D(xboxsize, nxgrid, origin=ymin)
        xs1, xs2 = MPI.split(len(xgrid))
        for i in range(1, MPI.size):
            _xmin = xedges[xs1[i]]-buffer_length
            _xmax = xedges[xs2[i]]+buffer_length
            cond = np.where((_x >= _xmin) & (_x < _xmax))[0]
            if f is None:
                _data = coords.coord2points([_x[cond], _y[cond]])
            else:
                _data = coords.coord2points([_x[cond], _y[cond], _f[cond]])
            MPI.send(_data, to_rank=i, tag=10+i)
        _xmin = xedges[xs1[0]]-buffer_length
        _xmax = xedges[xs2[0]]+buffer_length
        cond = np.where((_x >= _xmin) & (_x < _xmax))[0]
        _x, _y = _x[cond], _y[cond]
        if f is not None:
            _f = _f[cond]
    else:
        _data = MPI.recv(0, tag=10+MPI.rank)
        _x, _y = _data[:,0], _data[:,1]
        if len(_data[0]) == 3:
            _f = data[:,2]
    if f is None:
        _f = None
    MPI.wait()
    dx = xboxsize/nxgrid
    dy = yboxsize/nygrid
    xedges, xgrid = shift.cart.mpi_grid1D(xboxsize, nxgrid, MPI, origin=xmin)
    yedges, ygrid = shift.cart.grid1D(yboxsize, nygrid, origin=ymin)
    x2D, y2D = shift.cart.mpi_grid2D(boxsize, ngrid, MPI, origin=origin)
    shape = np.shape(x2D)
    nxgrid, nygrid = shape[0], shape[1]
    f2D = np.zeros((nxgrid, nygrid))
    if np.isscalar(MPI_ngrid):
        xs1, xs2 = 0, len(xgrid)
        ys1, ys2 = MPI.split(len(ygrid), size=MPI_ngrid)
    else:
        xs1, xs2 = MPI.split(len(xgrid), size=MPI_ngrid[0])
        ys1, ys2 = MPI.split(len(ygrid), size=MPI_ngrid[1])
    xs1, ys1 = np.meshgrid(xs1, ys1, indexing='ij')
    xs1, ys1 = xs1.flatten(), ys1.flatten()
    xs2, ys2 = np.meshgrid(xs2, ys2, indexing='ij')
    xs2, ys2 = xs2.flatten(), ys2.flatten()
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


def mpi_dtfe4grid3D(x, y, z, ngrid, boxsize, MPI, MPI_ngrid, f=None, origin=0.,
                    buffer_type=None, buffer_factor=6, buffer_val=0., subsampling=4,
                    outputgrid=False):
    if np.isscalar(ngrid):
        nxgrid, nygrid, nzgrid = ngrid, ngrid, ngrid
    else:
        nxgrid, nygrid, nzgrid = ngrid[0], ngrid[1], ngrid[2]
    if np.isscalar(boxsize):
        xboxsize, yboxsize, zboxsize = boxsize, boxsize, boxsize
    else:
        xboxsize, yboxsize, zboxsize = boxsize[0], boxsize[1], boxsize[2]
    if np.isscalar(boxsize):
        xmin, ymin, zmin = origin, origin, origin
    else:
        xmin, ymin, zmin = origin[0], origin[1], origin[2]
    if MPI.rank == 0:
        if buffer_type is not None:
            buffer_length = buffer_factor*dtfe4grid.mean_separation_3D(len(x), boxsize)
        else:
            buffer_length = 0.
        MPI.send(buffer_length, tag=11)
    else:
        buffer_length = MPI.recv(0, tag=11)
    MPI.wait()
    if MPI.rank == 0:
        if f is None:
            data = coords.coord2points([x, y, z])
        else:
            data = coords.coord2points([x, y, z, f])
        if buffer_type == 'random':
            xb, yb, zb = boundary.buffer_random_3D(len(x), boxsize, buffer_length, origin=origin)
            _x = np.concatenate([data[:,0], xb])
            _y = np.concatenate([data[:,1], yb])
            _z = np.concatenate([data[:,2], zb])
            if f is not None:
                _f = np.concatenate([data[:,3], buffer_val*np.ones(len(xb))])
        elif buffer_type == 'periodic':
            if f is None:
                data = coords.coord2points([x, y, z])
            else:
                data = coords.coord2points([x, y, z, f])
            datab = boundary.buffer_periodic_3D(data, boxsize, buffer_length, origin=origin)
            _data = np.vstack([data, datab])
            _x, _y, _z = _data[:,0], _data[:,1], _data[:,2]
            if f is not None:
                _f = _data[:,3]
        else:
            _x, _y, _z = x, y, z
            if f is not None:
                _f = f
        xedges, xgrid = shift.cart.grid1D(xboxsize, nxgrid, origin=ymin)
        xs1, xs2 = MPI.split(len(xgrid))
        for i in range(1, MPI.size):
            _xmin = xedges[xs1[i]]-buffer_length
            _xmax = xedges[xs2[i]]+buffer_length
            cond = np.where((_x >= _xmin) & (_x < _xmax))[0]
            if f is None:
                _data = coords.coord2points([_x[cond], _y[cond], _z[cond]])
            else:
                _data = coords.coord2points([_x[cond], _y[cond], _z[cond], _f[cond]])
            MPI.send(_data, to_rank=i, tag=10+i)
        _xmin = xedges[xs1[0]]-buffer_length
        _xmax = xedges[xs2[0]]+buffer_length
        cond = np.where((_x >= _xmin) & (_x < _xmax))[0]
        _x, _y, _z = _x[cond], _y[cond], _z[cond]
        if f is not None:
            _f = _f[cond]
    else:
        _data = MPI.recv(0, tag=10+MPI.rank)
        _x, _y, _z = _data[:,0], _data[:,1], _data[:,2]
        if len(_data[0]) == 4:
            _f = _data[:,3]
    if f is None:
        _f = None
    MPI.wait()
    dx = xboxsize/nxgrid
    dy = yboxsize/nygrid
    dz = zboxsize/nzgrid
    xedges, xgrid = shift.cart.mpi_grid1D(xboxsize, nxgrid, MPI, origin=xmin)
    yedges, ygrid = shift.cart.grid1D(yboxsize, nygrid, origin=ymin)
    zedges, zgrid = shift.cart.grid1D(zboxsize, nzgrid, origin=zmin)
    x3D, y3D, z3D = shift.cart.mpi_grid3D(boxsize, ngrid, MPI, origin=origin)
    shape = np.shape(x3D)
    nxgrid, nygrid, nzgrid = shape[0], shape[1], shape[2]
    f3D = np.zeros((nxgrid, nygrid, nzgrid))
    if np.isscalar(MPI_ngrid):
        xs1, xs2 = 0, len(xgrid)
        ys1, ys2 = MPI.split(len(ygrid), size=MPI_ngrid)
        zs1, zs2 = MPI.split(len(zgrid), size=MPI_ngrid)
    else:
        xs1, xs2 = MPI.split(len(xgrid), size=MPI_ngrid[0])
        ys1, ys2 = MPI.split(len(ygrid), size=MPI_ngrid[1])
        zs1, zs2 = MPI.split(len(zgrid), size=MPI_ngrid[2])
    xs1, ys1, zs1 = np.meshgrid(xs1, ys1, zs1, indexing='ij')
    xs1, ys1, zs1 = xs1.flatten(), ys1.flatten(), zs1.flatten()
    xs2, ys2, zs2 = np.meshgrid(xs2, ys2, zs2, indexing='ij')
    xs2, ys2, zs2 = xs2.flatten(), ys2.flatten(), zs2.flatten()
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
        if len(_xx) > 0:
            _f3d = dtfe4grid.dtfe4grid3D(_xx, _yy, _zz, [_nxgrid, _nygrid, _nzgrid],
                                         [_xmax-_xmin, _ymax-_ymin, _zmax-_zmin],
                                         f=_ff, origin=[_xmin, _ymin, _zmin],
                                         buffer_type=None, subsampling=subsampling,
                                         outputgrid=False)
            f3D[xs1[i]:xs2[i],ys1[i]:ys2[i],zs1[i]:zs2[i]] = _f3d
    if outputgrid:
        return x3D, y3D, f3D
    else:
        return f3D
