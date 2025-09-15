import numpy as np

import shift

from . import dtfe2d
from . import dtfe2dv2
from . import dtfe3d


def mean_separation_2D(npart, boxsize):
    """Calculates the mean separation for particles inside a 2D box.

    Parameters
    ----------
    npart : int
        Number of particles.
    boxsize : float or list
        Dimensions of the box.

    Returns
    -------
    mean_sep : float
        Mean separation of particles.
    """
    if np.isscalar(boxsize):
        area = boxsize**2.
    else:
        area = boxsize[0]*boxsize[1]
    mean_sep = np.sqrt(area/npart)
    return mean_sep


def mean_separation_3D(npart, boxsize):
    """Calculates the mean separation for particles inside a 3D box.

    Parameters
    ----------
    npart : int
        Number of particles.
    boxsize : float or list
        Dimensions of the box.

    Returns
    -------
    mean_sep : float
        Mean separation of particles.
    """
    if np.isscalar(boxsize):
        vol = boxsize**3.
    else:
        vol = boxsize[0]*boxsize[1]*boxsize[2]
    mean_sep = (vol/npart)**(1./3.)
    return mean_sep


def dtfe4grid2D(x, y, ngrid, boxsize, f=None, mass=None, origin=0., buffer_type=None,
                buffer_length=0., buffer_val=0., buffer_mass=None, subsampling=4, 
                outputgrid=False, calcdens=True):
    """Returns the Delaunay tesselation density or field on a grid.

    Parameters
    ----------
    x, y : array
        Coordinates of particles.
    ngrid : int or int list
        Grid dimensions.
    boxsize : float or list
        Dimensions of the grid.
    f : array, optional
        Field values, if None assumed output is density.
    mass : array, optional
        Mass of the particles.
    origin : float or list, optional
        Origin for grid.
    buffer_type : str, optional
        Buffer particle type, either:
            - 'random' for random buffer particles.
            - 'periodic' for periodic buffer particles.
            - None for no buffer particles.
    buffer_length : float, optional
        Buffer length.
    buffer_val : float, optional
        Value given to random buffer particles.
    buffer_mass : float, optional
        Must be provided if mass is provided and buffer particles are provided.
    subsampling : int, optional
        The pixel subsampling rate. Each pixel is evaluated subsampling^2 points
        on a grid within each pixel. This is to ensure each pixel is assigned a
        mean pixel value and not the value at the center.
    outputgrid : bool, optional
        Outputs coordinate grid.
    calcdens : bool, optional
        Calculates density.

    Returns
    -------
    f2d : ndarray
        Field values on a grid.
    x2d, y2d : ndarray, optional
        Pixel coordinate points.
    """
    # define boxsize on each axis
    if np.isscalar(boxsize):
        xboxsize, yboxsize = boxsize, boxsize
    else:
        xboxsize, yboxsize = boxsize[0], boxsize[1]
    # define grid on each axis
    if np.isscalar(ngrid):
        nxgrid, nygrid = ngrid, ngrid
    else:
        nxgrid, nygrid = ngrid[0], ngrid[1]
    # define subsampling rate for each pixel across each axis
    if np.isscalar(subsampling):
        dnxgrid, dnygrid = subsampling, subsampling
    else:
        dnxgrid, dnygrid = subsampling[0], subsampling[1]
    # define grid coordinates
    x2d, y2d = shift.cart.grid2D(boxsize, ngrid, origin=origin)
    x2d, y2d = x2d.flatten(), y2d.flatten()
    # define pixel length across each axis
    dx = xboxsize/nxgrid
    dy = yboxsize/nygrid
    # define subsampling points for each pixel
    dx2d, dy2d = shift.cart.grid2D([dx, dy], [dnxgrid, dnygrid], origin=[-dx/2., -dy/2.])
    dx2d, dy2d = dx2d.flatten(), dy2d.flatten()
    # initialise Delaunay tesselation
    D2D = dtfe2d.Delaunay2D()
    # add points
    if mass is None:
        if f is None:
            D2D.set_points(x, y, np.ones(len(x)))
        else:
            D2D.set_points(x, y, f)
    else:
        if f is None:
            D2D.set_points(x, y, np.ones(len(x)), mass=mass)
        else:
            D2D.set_points(x, y, f, mass=mass)
    # set boundary buffer points, either periodic or random buffer points
    if buffer_type == 'periodic':
        D2D.set_periodic(boxsize, buffer_length)
    elif buffer_type == 'random':
        if mass is None:
            D2D.set_buffer(boxsize, buffer_length, buffer_val=buffer_val)
        else:
            if buffer_mass is not None:
                D2D.set_buffer(boxsize, buffer_length, buffer_val=buffer_val, buffer_mass=buffer_mass)
            else:
                assert False, "buffer_mass must be defined."
    # construct delaunay tesselation triangles
    D2D.construct()
    # calculate delaunay tesselation field, if f is None we compute the density
    if calcdens:
        D2D.get_dens()
        D2D.set_field(f=D2D.points_dens)
    else:
        D2D.set_field()
    # calculate the field at the grid points
    for i in range(0, len(dx2d)):
        _f2d = (1./len(dx2d))*D2D.estimate(x2d+dx2d[i], y2d+dy2d[i]).reshape(nxgrid, nygrid)
        if i == 0:
            f2d = _f2d
        else:
            f2d += _f2d
    # ouput the field on a grid
    if outputgrid:
        return x2d.reshape(nxgrid, nygrid), y2d.reshape(nxgrid, nygrid), f2d.reshape(nxgrid, nygrid)
    else:
        return f2d.reshape(nxgrid, nygrid)


def dtfe4grid3D(x, y, z, ngrid, boxsize, f=None, mass=None, origin=0., buffer_length=0.,
                buffer_val=0., buffer_mass=None, buffer_type=None, subsampling=4,
                outputgrid=False, calcdens=True):
    """Returns the Delaunay tesselation density or field on a grid.

    Parameters
    ----------
    x, y, z : array
        Coordinates of particles.
    ngrid : int or int list
        Grid dimensions.
    boxsize : float or list
        Dimensions of the grid.
    f : array, optional
        Field values, if None assumed output is density.
    mass : array, optional
        Mass of the particles.
    origin : float or list, optional
        Origin for grid.
    buffer_type : str, optional
        Buffer particle type, either:
            - 'random' for random buffer particles.
            - 'periodic' for periodic buffer particles.
            - None for no buffer particles.
    buffer_length : float, optional
        Buffer length.
    buffer_val : float, optional
        Value given to random buffer particles.
    buffer_mass : float, optional
        Must be provided if mass is provided and buffer particles are provided.
    subsampling : int, optional
        The pixel subsampling rate. Each pixel is evaluated subsampling^2 points
        on a grid within each pixel. This is to ensure each pixel is assigned a
        mean pixel value and not the value at the center.
    outputgrid : bool, optional
        Outputs coordinate grid.
    calcdens : bool, optional
        Calculates density.

    Returns
    -------
    f3d : ndarray
        Field values on a grid.
    x3d, y3d, z3d : ndarray, optional
        Pixel coordinate points.
    """
    # define boxsize on each axis
    if np.isscalar(boxsize):
        xboxsize, yboxsize, zboxsize = boxsize, boxsize, boxsize
    else:
        xboxsize, yboxsize, zboxsize = boxsize[0], boxsize[1], boxsize[2]
    # define grid on each axis
    if np.isscalar(ngrid):
        nxgrid, nygrid, nzgrid = ngrid, ngrid, ngrid
    else:
        nxgrid, nygrid, nzgrid = ngrid[0], ngrid[1], ngrid[2]
    # define subsampling rate for each pixel across each axis
    if np.isscalar(subsampling):
        dnxgrid, dnygrid, dnzgrid = subsampling, subsampling, subsampling
    else:
        dnxgrid, dnygrid, dnzgrid = subsampling[0], subsampling[1], subsampling[2]
    # define grid coordinates
    x3d, y3d, z3d = shift.cart.grid3D(boxsize, ngrid, origin=origin)
    x3d, y3d, z3d = x3d.flatten(), y3d.flatten(), z3d.flatten()
    # define pixel length across each axis
    dx = xboxsize/nxgrid
    dy = yboxsize/nygrid
    dz = zboxsize/nzgrid
    # define subsampling points for each pixel
    dx3d, dy3d, dz3d = shift.cart.grid3D([dx, dy, dz], [dnxgrid, dnygrid, dnzgrid], origin=[-dx/2., -dy/2., -dz/2.])
    dx3d, dy3d, dz3d = dx3d.flatten(), dy3d.flatten(), dz3d.flatten()
    # initialise Delaunay tesselation
    D3D = dtfe3d.Delaunay3D()
    # add points
    if mass is None:
        if f is None:
            D3D.set_points(x, y, z, np.ones(len(x)))
        else:
            D3D.set_points(x, y, z, f)
    else:
        if f is None:
            D3D.set_points(x, y, z, np.ones(len(x)), mass=mass)
        else:
            D3D.set_points(x, y, z, f, mass=mass)
    # set boundary buffer points, either periodic or random buffer points
    if buffer_type == 'periodic':
        D3D.set_periodic(boxsize, buffer_length)
    elif buffer_type == 'random':
        if mass is None:
            D3D.set_buffer(boxsize, buffer_length, buffer_val=buffer_val)
        else:
            if buffer_mass is not None:
                D3D.set_buffer(boxsize, buffer_length, buffer_val=buffer_val, buffer_mass=buffer_mass)
            else:
                assert False, "buffer_mass must be defined."
    # construct delaunay tesselation triangles
    D3D.construct()
    # calculate delaunay tesselation field, if f is None we compute the density
    # if f is None:
    #     D3D.get_dens()
    #     D3D.set_field(f=D3D.points_dens)
    # else:
    #     D3D.set_field()
    if calcdens:
        D3D.get_dens()
        D3D.set_field(f=D3D.points_dens)
    else:
        D3D.set_field()
    # calculate the field at the grid points
    for i in range(0, len(dx3d)):
        _f3d = (1./len(dx3d))*D3D.estimate(x3d+dx3d[i], y3d+dy3d[i], z3d+dz3d[i])
        if i == 0:
            f3d = _f3d
        else:
            f3d += _f3d
    # ouput the field on a grid
    if outputgrid:
        return x3d.reshape(nxgrid, nygrid, nzgrid), y3d.reshape(nxgrid, nygrid, nzgrid), z3d.reshape(nxgrid, nygrid, nzgrid), f3d.reshape(nxgrid, nygrid, nzgrid)
    else:
        return f3d.reshape(nxgrid, nygrid, nzgrid)


def get_density_dtfe4grid2D(x, y, boxsize, ngrid, origin=0., mass=None, partition=1, periodic=True, fbuffer=0.2, subsampling=1, outputgrid=False):
    """Returns the Delaunay tesselation density or field on a grid.

    Parameters
    ----------
    x, y : array
        Coordinates of particles.
    boxsize : float or list
        Dimensions of the grid.
    ngrid : int or int list
        Grid dimensions.
    origin : float or list, optional
        Origin for grid.
    mass : array, optional
        Mass or mass weights for each the particles.
    partition : int or list, optional
        The number of internal Delaunay tesselations used to compute the density grid.
    periodic : bool or list, optional
        Determines whether periodic boundaries are applied.
    fbuffer : float, optional
        Buffer factor length.
    subsampling : int, optional
        The pixel subsampling rate. Each pixel is evaluated subsampling^2 points
        on a grid within each pixel. This is to ensure each pixel is assigned a
        mean pixel value and not the value at the center.
    outputgrid : bool, optional
        Outputs coordinate grid.

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
    
    # define grid coordinates
    _x2d, _y2d = shift.cart.grid2D(boxsize, ngrid, origin=[xorigin, yorigin])
    _x2d, _y2d = _x2d.flatten(), _y2d.flatten()
    _dens = np.zeros(len(_x2d))

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

    # partition data
    if np.isscalar(partition):
        xpartition, ypartition = partition, partition
    else:
        xpartition, ypartition = partition[0], partition[1]
    
    # construct partition divisions
    xedges, _ = shift.cart.grid1D(xboxsize, xpartition, origin=xorigin)
    yedges, _ = shift.cart.grid1D(yboxsize, ypartition, origin=yorigin)

    xedges[:-1] = np.floor(xedges[:-1] / dx) * dx
    yedges[:-1] = np.floor(yedges[:-1] / dy) * dy
    
    # construct partition indexes
    idx = np.arange(xpartition)
    idy = np.arange(ypartition)

    # assign mass weights
    if mass is None:
        mass = np.ones(len(x))

    for _dx2d, _dy2d in zip(dx2d, dy2d):

        x2d, y2d = _x2d + _dx2d, _y2d + _dy2d

        pixID = np.arange(len(x2d))
        dens = np.zeros(len(x2d))

        print(_dx2d, _dy2d, x2d, y2d)

        # lists for concatenating boundary particles
        x_boundary = []
        y_boundary = []
        mass_boundary = []
        type_boundary = []
        dens_boundary = []
        
        dtfe = dtfe2dv2.Delaunay2Dv2()

        # cycle through partition regions
        for i in range(0, len(idx)):
            
            xcond = np.where((x >= xedges[idx[i]]) & (x < xedges[idx[i]+1]))[0]

            pix_xcond = np.where((x2d[pixID] >= xedges[idx[i]]) & (x2d[pixID] < xedges[idx[i]+1]))[0]
            
            for j in range(0, len(idy)):
                
                # local dtfe
                dtfe.clean()
                
                ycond = np.where((y[xcond] >= yedges[idy[j]]) & (y[xcond] < yedges[idy[j]+1]))[0]
                
                # partition boundary
                boundary = [xedges[idx[i]], xedges[idx[i]+1], yedges[idy[j]], yedges[idy[j]+1]]
                    
                dtfe.setup(x[xcond[ycond]], y[xcond[ycond]], boundary, mass=mass[xcond[ycond]])

                dtfe.run()
                dtfe.get_dens()
                dtfe.set_field2dens()

                cond = np.where(dtfe.point_type != 1.)[0]
                
                # store boundary particles
                x_boundary.append(x[xcond[ycond[cond]]])
                y_boundary.append(y[xcond[ycond[cond]]])
                mass_boundary.append(mass[xcond[ycond[cond]]])
                type_boundary.append(dtfe.point_type[cond])
                dens_boundary.append(dtfe.point_dens[cond])

                pix_ycond = np.where((y2d[pixID[pix_xcond]] >= yedges[idy[j]]) & (y2d[pixID[pix_xcond]] < yedges[idy[j]+1]))[0]

                dens[pixID[pix_xcond[pix_ycond]]] = dtfe.estimate(x2d[pixID[pix_xcond[pix_ycond]]], y2d[pixID[pix_xcond[pix_ycond]]])

            cond = np.where(np.isfinite(dens[pixID[pix_xcond]]))[0]
            pixID = np.delete(pixID, pix_xcond[cond])

        # concatenate boundary particles
        x_boundary = np.concatenate(x_boundary)
        y_boundary = np.concatenate(y_boundary)
        mass_boundary = np.concatenate(mass_boundary)
        type_boundary = np.concatenate(type_boundary)
        dens_boundary = np.concatenate(dens_boundary)

        if len(idx) == 1 and len(idy) == 1:

            _x_border = x_boundary
            _y_border = y_boundary
            _mass_border = mass_boundary
            _type_border = type_boundary
            _dens_border = dens_boundary

        else:

            # lists for concatenating boundary particles globally
            x_gboundary = []
            y_gboundary = []
            mass_gboundary = []
            type_gboundary = []
            dens_gboundary = []

            for i in range(0, len(idx)):

                # partition boundary
                boundary = [xedges[idx[i]], xedges[idx[i]+1], yorigin, yorigin+yboxsize]

                # local dtfe
                dtfe.clean()
                
                cond = np.where((x_boundary >= boundary[0]) & (x_boundary < boundary[1]) & (y_boundary >= boundary[2]) & (y_boundary < boundary[3]))[0]

                dtfe.setup(x_boundary[cond], y_boundary[cond], boundary, mass=mass_boundary[cond])

                dtfe.run()
                dtfe.get_dens()

                cond2 = np.where(type_boundary[cond] == 0.)[0]
                cond3 = np.where(dtfe.point_type[cond2] != 1.)[0]

                dtfe.point_dens[cond2] = dens_boundary[cond[cond2]]
                dtfe.point_type[cond2[cond3]] = type_boundary[cond[cond2[cond3]]]

                dtfe.set_field2dens()
                
                cond2 = np.where(dtfe.point_type != 1.)[0]
                
                # store boundary particles
                x_gboundary.append(x_boundary[cond[cond2]])
                y_gboundary.append(y_boundary[cond[cond2]])
                mass_gboundary.append(mass_boundary[cond[cond2]])
                type_gboundary.append(dtfe.point_type[cond2])
                dens_gboundary.append(dtfe.point_dens[cond2])

                cond = np.where((x2d[pixID] >= boundary[0]) & (x2d[pixID] < boundary[1]) & (y2d[pixID] >= boundary[2]) & (y2d[pixID] < boundary[3]))[0]

                dens[pixID[cond]] = dtfe.estimate(x2d[pixID[cond]], y2d[pixID[cond]])

                cond2 = np.where(np.isfinite(dens[pixID[cond]]))[0]
                
                pixID = np.delete(pixID, cond[cond2])

            x_gboundary = np.concatenate(x_gboundary)
            y_gboundary = np.concatenate(y_gboundary)
            mass_gboundary = np.concatenate(mass_gboundary)
            type_gboundary = np.concatenate(type_gboundary)
            dens_gboundary = np.concatenate(dens_gboundary)

            # partition boundary
            boundary = [xorigin, xorigin+xboxsize, yorigin, yorigin+yboxsize]
            
            dtfe.clean()
            
            dtfe.setup(x_gboundary, y_gboundary, boundary, mass=mass_gboundary)

            dtfe.run()
            dtfe.get_dens()

            cond = np.where(type_gboundary[:] == 0.)[0]
            cond2 = np.where(dtfe.point_type[cond] != 1.)[0]

            dtfe.point_dens[cond] = dens_gboundary[cond]
            dtfe.point_type[cond[cond2]] = type_gboundary[cond[cond2]]

            dtfe.set_field2dens()
            
            cond = np.where(dtfe.point_type != 1.)[0]
            
            # store border particles
            _x_border = x_gboundary[cond]
            _y_border = y_gboundary[cond]
            _mass_border = mass_gboundary[cond]
            _type_border = dtfe.point_type[cond]
            _dens_border = dtfe.point_dens[cond]

            cond = np.where((x2d[pixID] >= boundary[0]) & (x2d[pixID] < boundary[1]) & (y2d[pixID] >= boundary[2]) & (y2d[pixID] < boundary[3]))[0]

            dens[pixID[cond]] = dtfe.estimate(x2d[pixID[cond]], y2d[pixID[cond]])

            cond2 = np.where(np.isfinite(dens[pixID[cond]]))[0]
            pixID = np.delete(pixID, cond[cond2])

        boundary = [xorigin, xorigin+xboxsize, yorigin, yorigin+yboxsize]

        if xperiodic == True or yperiodic == True:

            xbuffer = fbuffer*xboxsize
            ybuffer = fbuffer*yboxsize

            x_border = []
            y_border = []
            mass_border = []
            type_border = []
            dens_border = []

            if xperiodic == True and yperiodic == False:
                for i in range(-1, 2):
                    __x_border = _x_border + i*xboxsize
                    cond = np.where((__x_border >= boundary[0] - xbuffer) & (__x_border <= boundary[1] + xbuffer))[0]
                    x_border.append(__x_border[cond])
                    y_border.append(_y_border[cond])
                    mass_border.append(_mass_border[cond])
                    type_border.append(_type_border[cond])
                    dens_border.append(_dens_border[cond])
            elif xperiodic == False and yperiodic == True:
                for i in range(-1, 2):
                    __y_border = _y_border + i*yboxsize
                    cond = np.where((__y_border >= boundary[2] - ybuffer) & (__y_border <= boundary[3] + ybuffer))[0]
                    x_border.append(_x_border[cond])
                    y_border.append(__y_border[cond])
                    mass_border.append(_mass_border[cond])
                    type_border.append(_type_border[cond])
                    dens_border.append(_dens_border[cond])
            else:
                for i in range(-1, 2):
                    __x_border = _x_border + i*xboxsize
                    for j in range(-1, 2):
                        __y_border = _y_border + j*yboxsize
                        cond = np.where((__x_border >= boundary[0] - xbuffer) & (__x_border <= boundary[1] + xbuffer) & 
                                        (__y_border >= boundary[2] - ybuffer) & (__y_border <= boundary[3] + ybuffer))[0]
                        x_border.append(__x_border[cond])
                        y_border.append(__y_border[cond])
                        mass_border.append(_mass_border[cond])
                        type_border.append(_type_border[cond])
                        dens_border.append(_dens_border[cond])

            x_border = np.concatenate(x_border)
            y_border = np.concatenate(y_border)
            mass_border = np.concatenate(mass_border)
            type_border = np.concatenate(type_border)
            dens_border = np.concatenate(dens_border)

            dtfe.clean()
            
            boundary_border = np.copy(boundary)
            boundary_border[0] -= xbuffer
            boundary_border[1] += xbuffer
            boundary_border[2] -= ybuffer
            boundary_border[3] += ybuffer

            dtfe.setup(x_border, y_border, boundary_border, mass=mass_border)

            dtfe.run()
            dtfe.get_dens()
            
            cond = np.where(type_border == 0.)[0]
            cond2 = np.where(dtfe.point_type[cond] != 1.)[0]

            dtfe.point_dens[cond] = dens_border[cond]
            dtfe.point_type[cond[cond2]] = type_border[cond[cond2]]

            dtfe.set_field2dens()

            cond = np.where((x2d[pixID] >= boundary[0]) & (x2d[pixID] < boundary[1]) & (y2d[pixID] >= boundary[2]) & (y2d[pixID] < boundary[3]))[0]

            dens[pixID[cond]] = dtfe.estimate(x2d[pixID[cond]], y2d[pixID[cond]])

            cond2 = np.where(np.isfinite(dens[pixID[cond]]))[0]

            pixID = np.delete(pixID, cond[cond2])
        
        _dens += dens

    _dens /= len(dx2d)

    dens = _dens.reshape(nxgrid,nygrid)

    if outputgrid:
        return dens, x2d, y2d
    else:
        return dens
    
# TODO add field version of above function


def get_density_dtfe4grid2D_v2(x, y, boxsize, ngrid, origin=0., mass=None, partition=1, periodic=True, fbuffer=0.2, subsampling=1, outputgrid=False):
    """Returns the Delaunay tesselation density or field on a grid.

    Parameters
    ----------
    x, y : array
        Coordinates of particles.
    boxsize : float or list
        Dimensions of the grid.
    ngrid : int or int list
        Grid dimensions.
    origin : float or list, optional
        Origin for grid.
    mass : array, optional
        Mass or mass weights for each the particles.
    partition : int or list, optional
        The number of internal Delaunay tesselations used to compute the density grid.
    periodic : bool or list, optional
        Determines whether periodic boundaries are applied.
    fbuffer : float, optional
        Buffer factor length.
    subsampling : int, optional
        The pixel subsampling rate. Each pixel is evaluated subsampling^2 points
        on a grid within each pixel. This is to ensure each pixel is assigned a
        mean pixel value and not the value at the center.
    outputgrid : bool, optional
        Outputs coordinate grid.

    Returns
    -------
    dens : ndarray
        Density field values on a grid.
    x2d, y2d : ndarray, optional
        Pixel coordinate points.
    """

    # TODO figure out how to stitch Delaunay tesselations, this will require going back to the class, doing this
    # here seems problematic.

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
    
    # define grid coordinates
    x2d, y2d = shift.cart.grid2D(boxsize, ngrid, origin=[xorigin, yorigin])
    x2d, y2d = x2d.flatten(), y2d.flatten()
    dens = np.zeros(len(x2d))

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

    # partition data
    if np.isscalar(partition):
        xpartition, ypartition = partition, partition
    else:
        xpartition, ypartition = partition[0], partition[1]
    
    # construct partition divisions
    xedges, _ = shift.cart.grid1D(xboxsize, xpartition, origin=xorigin)
    yedges, _ = shift.cart.grid1D(yboxsize, ypartition, origin=yorigin)

    xedges[:-1] = np.floor(xedges[:-1] / dx) * dx
    yedges[:-1] = np.floor(yedges[:-1] / dy) * dy
    
    # construct partition indexes
    idx = np.arange(xpartition)
    idy = np.arange(ypartition)

    # assign mass weights
    if mass is None:
        mass = np.ones(len(x))

    pixID = np.arange(len(x2d))
    count = np.zeros(len(x2d))

    # lists for concatenating boundary particles
    x_boundary = []
    y_boundary = []
    mass_boundary = []
    type_boundary = []
    dens_boundary = []
    
    dtfe = dtfe2dv2.Delaunay2Dv2()

    # cycle through partition regions
    for i in range(0, len(idx)):
        
        pix_xcond = np.where((x2d[pixID] >= xedges[idx[i]]) & (x2d[pixID] < xedges[idx[i]+1]))[0]
        
        xcond = np.where((x >= xedges[idx[i]]) & (x < xedges[idx[i]+1]))[0]
        
        for j in range(0, len(idy)):
            
            # local dtfe
            dtfe.clean()
            
            ycond = np.where((y[xcond] >= yedges[idy[j]]) & (y[xcond] < yedges[idy[j]+1]))[0]
            
            # partition boundary
            boundary = [xedges[idx[i]], xedges[idx[i]+1], yedges[idy[j]], yedges[idy[j]+1]]
                
            dtfe.setup(x[xcond[ycond]], y[xcond[ycond]], boundary, mass=mass[xcond[ycond]])

            dtfe.run()
            dtfe.get_dens()
            dtfe.set_field2dens()

            cond = np.where(dtfe.point_type != 1.)[0]
            
            # store boundary particles
            x_boundary.append(x[xcond[ycond[cond]]])
            y_boundary.append(y[xcond[ycond[cond]]])
            mass_boundary.append(mass[xcond[ycond[cond]]])
            type_boundary.append(dtfe.point_type[cond])
            dens_boundary.append(dtfe.point_dens[cond])

            pix_ycond = np.where((y2d[pixID[pix_xcond]] >= yedges[idy[j]]) & (y2d[pixID[pix_xcond]] < yedges[idy[j]+1]))[0]

            for _dx2d, _dy2d in zip(dx2d, dy2d):
                dtfe_estimate = dtfe.estimate(x2d[pixID[pix_xcond[pix_ycond]]]+_dx2d, y2d[pixID[pix_xcond[pix_ycond]]]+_dy2d)
                cond = np.where(np.isfinite(dtfe_estimate) & (count[pix_xcond[pix_ycond]] != len(dx2d)))[0]
                dens[pixID[pix_xcond[pix_ycond[cond]]]] += dtfe_estimate[cond]
                count[pix_xcond[pix_ycond[cond]]] += 1.

        cond = np.where(count[pix_xcond] == len(dx2d))[0]
        pixID = np.delete(pixID, pix_xcond[cond])
        count = np.delete(count, pix_xcond[cond])

    # concatenate boundary particles
    x_boundary = np.concatenate(x_boundary)
    y_boundary = np.concatenate(y_boundary)
    mass_boundary = np.concatenate(mass_boundary)
    type_boundary = np.concatenate(type_boundary)
    dens_boundary = np.concatenate(dens_boundary)
    
    import matplotlib.pylab as plt

    plt.figure(figsize=(10,10))
    plt.imshow(np.log10(dens.reshape(nxgrid,nygrid)).T, origin='lower')
    plt.show()

    if len(idx) == 1 and len(idy) == 1:

        _x_border = x_boundary
        _y_border = y_boundary
        _mass_border = mass_boundary
        _type_border = type_boundary
        _dens_border = dens_boundary

    else:

        partID = np.searchsorted(yedges, y_boundary, side="right") - 1

        # lists for concatenating boundary particles globally
        x_gboundary = []
        y_gboundary = []
        mass_gboundary = []
        type_gboundary = []
        dens_gboundary = []

        for i in range(0, len(idx)):

            # partition boundary
            boundary = [xedges[idx[i]], xedges[idx[i]+1], yorigin, yorigin+yboxsize]

            # local dtfe
            dtfe.clean()
            
            cond = np.where((x_boundary >= boundary[0]) & (x_boundary < boundary[1]) & (y_boundary >= boundary[2]) & (y_boundary < boundary[3]))[0]

            dtfe.setup(x_boundary[cond], y_boundary[cond], boundary, mass=mass_boundary[cond])

            dtfe.run()
            dtfe.get_dens()

            cond2 = np.where(type_boundary[cond] == 0.)[0]

            dtfe.point_dens[cond2] = dens_boundary[cond[cond2]]
            dtfe.point_type[cond2] = type_boundary[cond[cond2]]

            dtfe.set_field2dens()

            # simplices_cond = np.where((dtfe.point_type[dtfe.simplices[:,0]] != 1.) & (dtfe.point_type[dtfe.simplices[:,1]] != 1.) & (dtfe.point_type[dtfe.simplices[:,2]] != 1))[0]
            # dtfe.simplices_type[simplices_cond] = -1.
            # simplices_cond = np.where(
            #     (dtfe.simplices_type == 1) &
            #     ((dtfe.point_type[dtfe.simplices[:,0]] == 0.) | (dtfe.point_type[dtfe.simplices[:,1]] == 0.) | (dtfe.point_type[dtfe.simplices[:,2]] == 0)) &
            #     (partID[dtfe.simplices[:,0]] == partID[dtfe.simplices[:,1]]) & 
            #     (partID[dtfe.simplices[:,1]] == partID[dtfe.simplices[:,2]])
            #     )[0]
            # dtfe.simplices_type[simplices_cond] = 1.
            # simplices_cond = np.where(
            #     (dtfe.simplices_type == 1) &
            #     ((dtfe.point_type[dtfe.simplices[:,0]] == 0.) & (dtfe.point_type[dtfe.simplices[:,1]] == 0.) & (dtfe.point_type[dtfe.simplices[:,2]] == 0)) &
            #     (partID[dtfe.simplices[:,0]] == partID[dtfe.simplices[:,1]]) & 
            #     (partID[dtfe.simplices[:,1]] == partID[dtfe.simplices[:,2]])
            #     )[0]
            # dtfe.simplices_type[simplices_cond] = 0.

            simplices_cond = np.where(
                ((partID[dtfe.simplices[:,0]] == partID[dtfe.simplices[:,1]]) & (partID[dtfe.simplices[:,1]] == partID[dtfe.simplices[:,2]])) & 
                (dtfe.simplices_type == 1.)
            )[0]
            dtfe.simplices_type[simplices_cond] = 1.

            # _cond = np.where(dtfe.simplices_type == 2.)[0]
            # _cond = np.unique(dtfe.simplices[_cond])
            # dtfe.point_type[_cond] = 2.
            _cond = np.where(dtfe.simplices_type == 1.)[0]
            _cond = np.unique(dtfe.simplices[_cond])
            dtfe.point_type[_cond] = 1.
            _cond = np.where(dtfe.simplices_type == 0.)[0]
            _cond = np.unique(dtfe.simplices[_cond])
            dtfe.point_type[_cond] = 0.
            _cond = np.where(dtfe.simplices_type == -1.)[0]
            _cond = np.unique(dtfe.simplices[_cond])
            dtfe.point_type[_cond] = -1.

            # cond_simplices = np.where(dtfe.simplices_type == -1)[0]
            # dtfe.point_type[np.unique(dtfe.simplices[cond_simplices])] = -1.
            
            cond2 = np.where(dtfe.point_type < 1.)[0]
            
            # store boundary particles
            x_gboundary.append(x_boundary[cond[cond2]])
            y_gboundary.append(y_boundary[cond[cond2]])
            mass_gboundary.append(mass_boundary[cond[cond2]])
            type_gboundary.append(dtfe.point_type[cond2])
            dens_gboundary.append(dtfe.point_dens[cond2])

            cond = np.where((x2d[pixID] >= boundary[0]) & (x2d[pixID] < boundary[1]) & (y2d[pixID] >= boundary[2]) & (y2d[pixID] < boundary[3]))[0]

            for _dx2d, _dy2d in zip(dx2d, dy2d):
                dtfe_estimate = dtfe.estimate(x2d[pixID[cond]]+_dx2d, y2d[pixID[cond]]+_dy2d)
                cond_isfinite = np.where(np.isfinite(dtfe_estimate))[0]# & (count[cond] != len(dx2d)))[0]
                dens[pixID[cond[cond_isfinite]]] += dtfe_estimate[cond_isfinite]
                count[cond[cond_isfinite]] += 1.

            cond2 = np.where(count[cond] == len(dx2d))[0]

            pixID = np.delete(pixID, cond[cond2])
            count = np.delete(count, cond[cond2])

            plt.figure(figsize=(10,10))
            plt.imshow(np.log10(dens.reshape(nxgrid,nygrid)).T, origin='lower')
            plt.show()

            # plt.scatter(x_boundary, y_boundary, c=type_boundary, s=1)
            # plt.show()

        x_gboundary = np.concatenate(x_gboundary)
        y_gboundary = np.concatenate(y_gboundary)
        mass_gboundary = np.concatenate(mass_gboundary)
        type_gboundary = np.concatenate(type_gboundary)
        dens_gboundary = np.concatenate(dens_gboundary)

        plt.scatter(x_gboundary, y_gboundary, c=type_gboundary, s=1)
        plt.show()

        partID = np.searchsorted(xedges, x_gboundary, side="right") - 1

        # partition boundary
        boundary = [xorigin, xorigin+xboxsize, yorigin, yorigin+yboxsize]
        
        dtfe.clean()
        
        dtfe.setup(x_gboundary, y_gboundary, boundary, mass=mass_gboundary)

        dtfe.run()
        dtfe.get_dens()

        cond = np.where(type_gboundary == 0.)[0]
        # cond2 = np.where(dtfe.point_type[cond] != 1.)[0]

        dtfe.point_dens[cond] = dens_gboundary[cond]
        # dtfe.point_type[cond[cond2]] = type_gboundary[cond[cond2]]
        dtfe.point_type[cond] = type_gboundary[cond]

        dtfe.set_field2dens()

        # simplices_cond = np.where((dtfe.point_type[dtfe.simplices[:,0]] != 1.) & (dtfe.point_type[dtfe.simplices[:,1]] != 1.) & (dtfe.point_type[dtfe.simplices[:,2]] != 1))[0]
        # dtfe.simplices_type[simplices_cond] = -1.
        # simplices_cond = np.where(
        #     (dtfe.simplices_type == 1) &
        #     ((dtfe.point_type[dtfe.simplices[:,0]] == 0.) | (dtfe.point_type[dtfe.simplices[:,1]] == 0.) | (dtfe.point_type[dtfe.simplices[:,2]] == 0)) &
        #     (partID[dtfe.simplices[:,0]] == partID[dtfe.simplices[:,1]]) & 
        #     (partID[dtfe.simplices[:,1]] == partID[dtfe.simplices[:,2]])
        #     )[0]
        # dtfe.simplices_type[simplices_cond] = 1.
        # simplices_cond = np.where(
        #     (dtfe.simplices_type == 1) &
        #     ((dtfe.point_type[dtfe.simplices[:,0]] == 0.) | (dtfe.point_type[dtfe.simplices[:,1]] == 0.) | (dtfe.point_type[dtfe.simplices[:,2]] == 0))
        #     )[0]
        # dtfe.simplices_type[simplices_cond] = 1.
        # simplices_cond = np.where(
        #     (dtfe.simplices_type == 1) &
        #     ((dtfe.point_type[dtfe.simplices[:,0]] == 0.) & (dtfe.point_type[dtfe.simplices[:,1]] == 0.) & (dtfe.point_type[dtfe.simplices[:,2]] == 0)) &
        #     (partID[dtfe.simplices[:,0]] == partID[dtfe.simplices[:,1]]) & 
        #     (partID[dtfe.simplices[:,1]] == partID[dtfe.simplices[:,2]])
        #     )[0]
        # dtfe.simplices_type[simplices_cond] = 0.

        simplices_cond = np.where(
            (dtfe.point_type[dtfe.simplices[:,0]] == 0.) &
            (dtfe.point_type[dtfe.simplices[:,1]] == 0.) &
            (dtfe.point_type[dtfe.simplices[:,2]] == 0.) &
            ((partID[dtfe.simplices[:,0]] == partID[dtfe.simplices[:,1]]) & 
             (partID[dtfe.simplices[:,1]] == partID[dtfe.simplices[:,2]])) & 
            (dtfe.simplices_type == 1.)
        )[0]
        dtfe.simplices_type[simplices_cond] = -1.

        # _cond = np.where(dtfe.simplices_type == 2.)[0]
        # _cond = np.unique(dtfe.simplices[_cond])
        # dtfe.point_type[_cond] = 2.
        _cond = np.where(dtfe.simplices_type == 1.)[0]
        _cond = np.unique(dtfe.simplices[_cond])
        dtfe.point_type[_cond] = 1.
        _cond = np.where(dtfe.simplices_type == 0.)[0]
        _cond = np.unique(dtfe.simplices[_cond])
        dtfe.point_type[_cond] = 0.
        _cond = np.where(dtfe.simplices_type == -1.)[0]
        _cond = np.unique(dtfe.simplices[_cond])
        dtfe.point_type[_cond] = -1.

        # cond_simplices = np.where(dtfe.simplices_type == -1)[0]
        # dtfe.point_type[np.unique(dtfe.simplices[cond_simplices])] = -1.

        cond = np.where(dtfe.point_type != 1.)[0]
        
        # store border particles
        _x_border = x_gboundary[cond]
        _y_border = y_gboundary[cond]
        _mass_border = mass_gboundary[cond]
        _type_border = dtfe.point_type[cond]
        _dens_border = dtfe.point_dens[cond]

        cond = np.where((x2d[pixID] >= boundary[0]) & (x2d[pixID] < boundary[1]) & (y2d[pixID] >= boundary[2]) & (y2d[pixID] < boundary[3]))[0]

        for _dx2d, _dy2d in zip(dx2d, dy2d):
            dtfe_estimate = dtfe.estimate(x2d[pixID[cond]]+_dx2d, y2d[pixID[cond]]+_dy2d)
            cond_isfinite = np.where(np.isfinite(dtfe_estimate))[0]# & (count[cond] != len(dx2d)))[0]
            dens[pixID[cond[cond_isfinite]]] += dtfe_estimate[cond_isfinite]
            count[cond[cond_isfinite]] += 1.

            # __dens = np.zeros(len(x2d))
            # __dens[pixID[cond]] = dtfe_estimate

        plt.figure(figsize=(10,10))
        plt.imshow(np.log10(dens.reshape(nxgrid,nygrid)).T, origin='lower')
        plt.show()
        
        cond2 = np.where(count[cond] == len(dx2d))[0]
        pixID = np.delete(pixID, cond[cond2])
        count = np.delete(count, cond[cond2])

    plt.figure(figsize=(10,10))
    plt.imshow(np.log10(dens.reshape(nxgrid,nygrid)).T, origin='lower')
    plt.show()

    boundary = [xorigin, xorigin+xboxsize, yorigin, yorigin+yboxsize]

    if xperiodic == True or yperiodic == True:

        xbuffer = fbuffer*xboxsize
        ybuffer = fbuffer*yboxsize

        x_border = []
        y_border = []
        mass_border = []
        type_border = []
        dens_border = []
        pID = 0
        partID = []

        if xperiodic == True and yperiodic == False:
            for i in range(-1, 2):
                __x_border = _x_border + i*xboxsize
                cond = np.where((__x_border >= boundary[0] - xbuffer) & (__x_border <= boundary[1] + xbuffer))[0]
                x_border.append(__x_border[cond])
                y_border.append(_y_border[cond])
                mass_border.append(_mass_border[cond])
                type_border.append(_type_border[cond])
                dens_border.append(_dens_border[cond])
                partID.append(np.ones(len(cond))*pID)
                pID += 1
        elif xperiodic == False and yperiodic == True:
            for i in range(-1, 2):
                __y_border = _y_border + i*yboxsize
                cond = np.where((__y_border >= boundary[2] - ybuffer) & (__y_border <= boundary[3] + ybuffer))[0]
                x_border.append(_x_border[cond])
                y_border.append(__y_border[cond])
                mass_border.append(_mass_border[cond])
                type_border.append(_type_border[cond])
                dens_border.append(_dens_border[cond])
                partID.append(np.ones(len(cond))*pID)
                pID += 1
        else:
            for i in range(-1, 2):
                __x_border = _x_border + i*xboxsize
                for j in range(-1, 2):
                    __y_border = _y_border + j*yboxsize
                    cond = np.where((__x_border >= boundary[0] - xbuffer) & (__x_border <= boundary[1] + xbuffer) & 
                                    (__y_border >= boundary[2] - ybuffer) & (__y_border <= boundary[3] + ybuffer))[0]
                    x_border.append(__x_border[cond])
                    y_border.append(__y_border[cond])
                    mass_border.append(_mass_border[cond])
                    type_border.append(_type_border[cond])
                    dens_border.append(_dens_border[cond])
                    partID.append(np.ones(len(cond))*pID)
                    pID += 1

        x_border = np.concatenate(x_border)
        y_border = np.concatenate(y_border)
        mass_border = np.concatenate(mass_border)
        type_border = np.concatenate(type_border)
        dens_border = np.concatenate(dens_border)
        partID = np.concatenate(partID)

        dtfe.clean()
        
        boundary_border = np.copy(boundary)
        boundary_border[0] -= xbuffer
        boundary_border[1] += xbuffer
        boundary_border[2] -= ybuffer
        boundary_border[3] += ybuffer

        dtfe.setup(x_border, y_border, boundary_border, mass=mass_border)

        dtfe.run()
        dtfe.get_dens()
        
        cond = np.where(type_border == 0.)[0]
        # cond2 = np.where(dtfe.point_type[cond] != 1.)[0]

        dtfe.point_dens[cond] = dens_border[cond]
        dtfe.point_type[cond] = type_border[cond]
        # dtfe.point_type[cond[cond2]] = type_border[cond[cond2]]

        dtfe.set_field2dens()

        # simplices_cond = np.where((dtfe.point_type[dtfe.simplices[:,0]] != 1.) & (dtfe.point_type[dtfe.simplices[:,1]] != 1.) & (dtfe.point_type[dtfe.simplices[:,2]] != 1))[0]
        # dtfe.simplices_type[simplices_cond] = -1.
        # simplices_cond = np.where(
        #     (dtfe.simplices_type == 1) &
        #     ((dtfe.point_type[dtfe.simplices[:,0]] == 0.) | (dtfe.point_type[dtfe.simplices[:,1]] == 0.) | (dtfe.point_type[dtfe.simplices[:,2]] == 0)) &
        #     (partID[dtfe.simplices[:,0]] == partID[dtfe.simplices[:,1]]) & 
        #     (partID[dtfe.simplices[:,1]] == partID[dtfe.simplices[:,2]])
        #     )[0]
        # dtfe.simplices_type[simplices_cond] = 1.
        # simplices_cond = np.where(
        #     (dtfe.simplices_type == 1) &
        #     ((dtfe.point_type[dtfe.simplices[:,0]] == 0.) & (dtfe.point_type[dtfe.simplices[:,1]] == 0.) & (dtfe.point_type[dtfe.simplices[:,2]] == 0)) &
        #     (partID[dtfe.simplices[:,0]] == partID[dtfe.simplices[:,1]]) & 
        #     (partID[dtfe.simplices[:,1]] == partID[dtfe.simplices[:,2]])
        #     )[0]
        # dtfe.simplices_type[simplices_cond] = 0.

        simplices_cond = np.where(
            ((partID[dtfe.simplices[:,0]] == partID[dtfe.simplices[:,1]]) & (partID[dtfe.simplices[:,1]] == partID[dtfe.simplices[:,2]])) & 
            (dtfe.simplices_type == 1.)
        )[0]
        dtfe.simplices_type[simplices_cond] = 1.

        # _cond = np.where(dtfe.simplices_type == 2.)[0]
        # _cond = np.unique(dtfe.simplices[_cond])
        # dtfe.point_type[_cond] = 2.
        _cond = np.where(dtfe.simplices_type == 1.)[0]
        _cond = np.unique(dtfe.simplices[_cond])
        dtfe.point_type[_cond] = 1.
        _cond = np.where(dtfe.simplices_type == 0.)[0]
        _cond = np.unique(dtfe.simplices[_cond])
        dtfe.point_type[_cond] = 0.
        _cond = np.where(dtfe.simplices_type == -1.)[0]
        _cond = np.unique(dtfe.simplices[_cond])
        dtfe.point_type[_cond] = -1.

        # cond_simplices = np.where(dtfe.simplices_type == -1)[0]
        # dtfe.point_type[np.unique(dtfe.simplices[cond_simplices])] = -1.

        cond = np.where((x2d[pixID] >= boundary[0]) & (x2d[pixID] < boundary[1]) & (y2d[pixID] >= boundary[2]) & (y2d[pixID] < boundary[3]))[0]

        # dens[pixID[cond]] = dtfe.estimate(x2d[pixID[cond]], y2d[pixID[cond]])

        # cond2 = np.where(np.isfinite(dens[pixID[cond]]))[0]

        for _dx2d, _dy2d in zip(dx2d, dy2d):
            dtfe_estimate = dtfe.estimate(x2d[pixID[cond]]+_dx2d, y2d[pixID[cond]]+_dy2d)
            cond_isfinite = np.where(np.isfinite(dtfe_estimate) & (count[cond] != len(dx2d)))[0]
            dens[pixID[cond[cond_isfinite]]] += dtfe_estimate[cond_isfinite]
            count[cond[cond_isfinite]] += 1.

        cond2 = np.where(count[cond] == len(dx2d))[0]
        pixID = np.delete(pixID, cond[cond2])
        count = np.delete(count, cond[cond2])

        plt.figure(figsize=(10,10))
        plt.imshow(np.log10(dens.reshape(nxgrid,nygrid)).T, origin='lower')
        plt.show()

    dens /= len(dx2d)

    dens = dens.reshape(nxgrid,nygrid)

    if outputgrid:
        return dens, x2d, y2d
    else:
        return dens