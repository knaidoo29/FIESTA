import numpy as np

import shift

from . import dtfe2d
from . import dtfe3d

from typing import Union, Tuple, Optional, List
from copy import deepcopy


# TODO: Get rid of this? Is this used for anything?

# def mean_separation_2D(npart, boxsize):
#     """
#     Calculates the mean separation for particles inside a 2D box.

#     Parameters
#     ----------
#     npart : int
#         Number of particles.
#     boxsize : float or list
#         Dimensions of the box.

#     Returns
#     -------
#     mean_sep : float
#         Mean separation of particles.
#     """
#     if np.isscalar(boxsize):
#         area = boxsize**2.
#     else:
#         area = boxsize[0]*boxsize[1]
#     mean_sep = np.sqrt(area/npart)
#     return mean_sep


# def mean_separation_3D(npart, boxsize):
#     """
#     Calculates the mean separation for particles inside a 3D box.

#     Parameters
#     ----------
#     npart : int
#         Number of particles.
#     boxsize : float or list
#         Dimensions of the box.

#     Returns
#     -------
#     mean_sep : float
#         Mean separation of particles.
#     """
#     if np.isscalar(boxsize):
#         vol = boxsize**3.
#     else:
#         vol = boxsize[0]*boxsize[1]*boxsize[2]
#     mean_sep = (vol/npart)**(1./3.)
#     return mean_sep


def delaunay_density4grid2D(
    x: np.ndarray, y: np.ndarray, boxsize: Union[float, List[float]], ngrid: Union[int, List[int]], 
    origin: Union[float, List[float]] = 0., mass: Optional[np.ndarray] = None, partition: int = 1, 
    periodic : Union[bool, List[bool]] = True, fbuffer: float = 0.2, subsampling: int = 1, 
    outputgrid: bool = False, outputexterior: bool = False, normalise : bool = True, flatten: bool = False
) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray, np.ndarray], 
           Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray], 
           Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
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
        Output exterior information, including exterior border dtfe and current unassigned pixels.
    normalise : bool, optional
        Whether to normalise the outputing density field.
    flatten : bool, optional
        Flatten output density array.

    Returns
    -------
    dens : ndarray
        Density field values on a grid.
    x2d, y2d : ndarray, optional
        Pixel coordinate points.
    exterior_border : dict, optional
        If outputexterior = True, the exterior_border dictionary is outputted.
    pixID : ndarray, optional
        If outputexterior = True, pixels that are currently undersampled.
    count: ndarray, optional
        If outputexterior = True, the counts for each undersampled pixel.
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

    # list of dictionary containing border points and simplices
    dict_boundary = []
    
    dtfe = dtfe2d.Delaunay2D()

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
            
            pix_ycond = np.where((y2d[pixID[pix_xcond]] >= yedges[idy[j]]) & (y2d[pixID[pix_xcond]] < yedges[idy[j]+1]))[0]

            for _dx2d, _dy2d in zip(dx2d, dy2d):
                dtfe_estimate = dtfe.estimate(x2d[pixID[pix_xcond[pix_ycond]]]+_dx2d, y2d[pixID[pix_xcond[pix_ycond]]]+_dy2d)
                cond = np.where(np.isfinite(dtfe_estimate) & (count[pix_xcond[pix_ycond]] != len(dx2d)))[0]
                dens[pixID[pix_xcond[pix_ycond[cond]]]] += dtfe_estimate[cond]
                count[pix_xcond[pix_ycond[cond]]] += 1.
            
            dict_boundary.append(dtfe.get_border())

        cond = np.where(count[pix_xcond] == len(dx2d))[0]
        pixID = np.delete(pixID, pix_xcond[cond])
        count = np.delete(count, pix_xcond[cond])

    dtfe = dtfe2d.DelaunayMerger2D()
    
    if len(idx) != 1 or len(idy) != 1:

        dict_boundary2 = []

        for i in range(0, len(idx)):

            # local dtfe
            dtfe.clean()

            for j in range(0, len(idy)):

                dtfe.add_border(dict_boundary[i*len(idy) + j])
            
            # partition boundary
            boundary = [xedges[idx[i]], xedges[idx[i]+1], yorigin, yorigin+yboxsize]

            dtfe.run(boundary)

            cond = np.where((x2d[pixID] >= boundary[0]) & (x2d[pixID] < boundary[1]) & (y2d[pixID] >= boundary[2]) & (y2d[pixID] < boundary[3]))[0]

            for _dx2d, _dy2d in zip(dx2d, dy2d):
                dtfe_estimate = dtfe.estimate(x2d[pixID[cond]]+_dx2d, y2d[pixID[cond]]+_dy2d)
                cond_isfinite = np.where(np.isfinite(dtfe_estimate) & (count[cond] != len(dx2d)))[0]
                dens[pixID[cond[cond_isfinite]]] += dtfe_estimate[cond_isfinite]
                count[cond[cond_isfinite]] += 1.

            cond2 = np.where(count[cond] == len(dx2d))[0]

            pixID = np.delete(pixID, cond[cond2])
            count = np.delete(count, cond[cond2])
            
            dict_boundary2.append(dtfe.get_border())

        dtfe.clean()
        
        for boundary in dict_boundary2:
            dtfe.add_border(boundary)

        # partition boundary
        boundary = [xorigin, xorigin+xboxsize, yorigin, yorigin+yboxsize]
        
        dtfe.run(boundary)

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
    else:
        exterior_border = dict_boundary[0]

    boundary = [xorigin, xorigin+xboxsize, yorigin, yorigin+yboxsize]

    _exterior_border = []

    if xperiodic == True or yperiodic == True:

        xbuffer = fbuffer*xboxsize
        ybuffer = fbuffer*yboxsize
        
        if xperiodic == True:
            imin, imax = -1, 2
        else:
            imin, imax = 0, 1
        if yperiodic == True:
            jmin, jmax = -1, 2
        else:
            jmin, jmax = 0, 1
        
        for i in range(imin, imax):
            for j in range(jmin, jmax):
                _border = deepcopy(exterior_border)
                _border['x'] += i*xboxsize
                _border['y'] += j*yboxsize
                _exterior_border.append(_border)
                    
        dtfe.clean()

        for border in _exterior_border:
            dtfe.add_border(border)

        boundary = [xorigin-xbuffer, xorigin+xboxsize+xbuffer, yorigin-ybuffer, yorigin+yboxsize+ybuffer]

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

    if normalise:
        dens /= len(dx2d)
        cond = np.where(count == 0.)[0]
        dens[pixID[cond]] *= len(dx2d)
        cond = np.where(count != 0.)[0]
        dens[pixID[cond]] *= len(dx2d)/count[cond]

    if flatten == False:
        dens = dens.reshape(nxgrid,nygrid)

    if outputexterior:
        if outputgrid:
            return dens, x2d, y2d, exterior_border, pixID, count
        else:
            return dens, exterior_border, pixID, count
    else:
        if outputgrid:
            return dens, x2d, y2d
        else:
            return dens


def delaunay_field4grid2D(
    x: np.ndarray, y: np.ndarray, f: np.ndarray, boxsize: Union[float, List[float]], ngrid: Union[int, List[int]], 
    origin: Union[float, List[float]] = 0., mass: Optional[np.ndarray] = None, partition: int = 1, 
    periodic : Union[bool, List[bool]] = True, fbuffer: float = 0.5, subsampling: int = 1, outputgrid: bool = False,
    outputexterior: bool = False, normalise: bool = True, flatten: bool = False
) -> Union[
    np.ndarray,
    Tuple[np.ndarray, np.ndarray, np.ndarray],
    Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]
]:
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
    normalise : bool, optional
        Whether to normalise the outputing field.
    flatten : bool, optional
        Flatten output field array.

    Returns
    -------
    field : ndarray
        Field values on a grid.
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
    x2d, y2d = shift.cart.grid2D(boxsize, ngrid, origin=[xorigin, yorigin])
    x2d, y2d = x2d.flatten(), y2d.flatten()
    field = np.zeros(len(x2d))

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

    # list of dictionary containing border points and simplices
    dict_boundary = []
    
    dtfe = dtfe2d.Delaunay2D()

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
            dtfe.set_field(f[xcond[ycond]])
            
            pix_ycond = np.where((y2d[pixID[pix_xcond]] >= yedges[idy[j]]) & (y2d[pixID[pix_xcond]] < yedges[idy[j]+1]))[0]

            for _dx2d, _dy2d in zip(dx2d, dy2d):
                dtfe_estimate = dtfe.estimate(x2d[pixID[pix_xcond[pix_ycond]]]+_dx2d, y2d[pixID[pix_xcond[pix_ycond]]]+_dy2d)
                cond = np.where(np.isfinite(dtfe_estimate) & (count[pix_xcond[pix_ycond]] != len(dx2d)))[0]
                field[pixID[pix_xcond[pix_ycond[cond]]]] += dtfe_estimate[cond]
                count[pix_xcond[pix_ycond[cond]]] += 1.
            
            dict_boundary.append(dtfe.get_border())

        cond = np.where(count[pix_xcond] == len(dx2d))[0]
        pixID = np.delete(pixID, pix_xcond[cond])
        count = np.delete(count, pix_xcond[cond])

    dtfe = dtfe2d.DelaunayMerger2D()

    if len(idx) != 1 or len(idy) != 1:

        dict_boundary2 = []

        for i in range(0, len(idx)):

            # local dtfe
            dtfe.clean()

            for j in range(0, len(idy)):

                dtfe.add_border(dict_boundary[i*len(idy) + j])
            
            # partition boundary
            boundary = [xedges[idx[i]], xedges[idx[i]+1], yorigin, yorigin+yboxsize]

            dtfe.run(boundary)

            cond = np.where((x2d[pixID] >= boundary[0]) & (x2d[pixID] < boundary[1]) & (y2d[pixID] >= boundary[2]) & (y2d[pixID] < boundary[3]))[0]

            for _dx2d, _dy2d in zip(dx2d, dy2d):
                dtfe_estimate = dtfe.estimate(x2d[pixID[cond]]+_dx2d, y2d[pixID[cond]]+_dy2d)
                cond_isfinite = np.where(np.isfinite(dtfe_estimate) & (count[cond] != len(dx2d)))[0]
                field[pixID[cond[cond_isfinite]]] += dtfe_estimate[cond_isfinite]
                count[cond[cond_isfinite]] += 1.

            cond2 = np.where(count[cond] == len(dx2d))[0]

            pixID = np.delete(pixID, cond[cond2])
            count = np.delete(count, cond[cond2])

            dict_boundary2.append(dtfe.get_border())

        dtfe.clean()
        
        for boundary in dict_boundary2:
            dtfe.add_border(boundary)

        # partition boundary
        boundary = [xorigin, xorigin+xboxsize, yorigin, yorigin+yboxsize]
        
        dtfe.run(boundary)

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
    else:
        exterior_border = dict_boundary[0]

    boundary = [xorigin, xorigin+xboxsize, yorigin, yorigin+yboxsize]

    _exterior_border = []

    if xperiodic == True or yperiodic == True:

        xbuffer = fbuffer*xboxsize
        ybuffer = fbuffer*yboxsize

        if xperiodic == True:
            imin, imax = -1, 2
        else:
            imin, imax = 0, 1
        if yperiodic == True:
            jmin, jmax = -1, 2
        else:
            jmin, jmax = 0, 1
        
        for i in range(imin, imax):
            for j in range(jmin, jmax):
                _border = deepcopy(exterior_border)
                _border['x'] += i*xboxsize
                _border['y'] += j*yboxsize
                _exterior_border.append(_border)
                    
        dtfe.clean()

        for border in _exterior_border:
            dtfe.add_border(border)

        boundary = [xorigin-xbuffer, xorigin+xboxsize+xbuffer, yorigin-ybuffer, yorigin+yboxsize+ybuffer]

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

    if normalise:
        field /= len(dx2d)

        cond = np.where(count == 0.)[0]
        field[pixID[cond]] *= len(dx2d)
        cond = np.where(count != 0.)[0]
        field[pixID[cond]] *= len(dx2d)/count[cond]

    if flatten == False:
        field = field.reshape(nxgrid,nygrid)

    if outputexterior:
        if outputgrid:
            return field, x2d, y2d, exterior_border, pixID, count
        else:
            return field, exterior_border, pixID, count
    else:
        if outputgrid:
            return field, x2d, y2d
        else:
            return field


def delaunay_density4grid3D(
    x: np.ndarray, y: np.ndarray, z: np.ndarray, boxsize: Union[float, List[float]], ngrid: Union[int, List[int]], 
    origin: Union[float, List[float]] = 0., mass: Optional[np.ndarray] = None, partition: int = 1, 
    periodic : Union[bool, List[bool]] = True, fbuffer: float = 0.5, subsampling: int = 1, 
    outputgrid: bool = False, outputexterior: bool = False, normalise: bool = True, flatten: bool = False
) -> Union[
    np.ndarray, 
    Tuple[np.ndarray, np.ndarray, np.ndarray],
    Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]    
]:
    """
    Returns the Delaunay tesselation density on a grid.

    Parameters
    ----------
    x, y, z : array
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
        Buffer factor length, used to trim exterior buffer region.
    subsampling : int, optional
        The pixel subsampling rate. Each pixel is evaluated subsampling^2 points
        on a grid within each pixel. This is to ensure each pixel is assigned a
        mean pixel value and not the value at the center.
    outputgrid : bool, optional
        Outputs coordinate grid.
    normalise : bool, optional
        Whether to normalise the outputing density field.
    flatten : bool, optional
        Flatten output density array.
    
    Returns
    -------
    dens : ndarray
        Density field values on a grid.
    x2d, y2d, z3d : ndarray, optional
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
    
    # define grid coordinates
    x3d, y3d, z3d = shift.cart.grid3D(boxsize, ngrid, origin=[xorigin, yorigin, zorigin])
    x3d, y3d, z3d = x3d.flatten(), y3d.flatten(), z3d.flatten()
    dens = np.zeros(len(x3d))

    # define pixel length across each axis
    dx = xboxsize/nxgrid
    dy = yboxsize/nygrid
    dz = zboxsize/nzgrid

    if np.isscalar(subsampling):
        dnxgrid, dnygrid, dnzgrid = subsampling, subsampling, subsampling
    else:
        dnxgrid, dnygrid, dnzgrid = subsampling[0], subsampling[1], subsampling[2]

    # define subsampling points for each pixel
    dx3d, dy3d, dz3d = shift.cart.grid3D([dx, dy, dz], [dnxgrid, dnygrid, dnzgrid], origin=[-dx/2., -dy/2., -dz/2.])
    dx3d, dy3d, dz3d = dx3d.flatten(), dy3d.flatten(), dz3d.flatten()

    # partition data
    if np.isscalar(partition):
        xpartition, ypartition, zpartition = partition, partition, partition
    else:
        xpartition, ypartition, zpartition = partition[0], partition[1], partition[2]
    
    # construct partition divisions
    xedges, _ = shift.cart.grid1D(xboxsize, xpartition, origin=xorigin)
    yedges, _ = shift.cart.grid1D(yboxsize, ypartition, origin=yorigin)
    zedges, _ = shift.cart.grid1D(zboxsize, zpartition, origin=zorigin)

    xedges[:-1] = np.floor(xedges[:-1] / dx) * dx
    yedges[:-1] = np.floor(yedges[:-1] / dy) * dy
    zedges[:-1] = np.floor(zedges[:-1] / dz) * dz
    
    # construct partition indexes
    idx = np.arange(xpartition)
    idy = np.arange(ypartition)
    idz = np.arange(zpartition)

    # assign mass weights
    if mass is None:
        mass = np.ones(len(x))

    pixID = np.arange(len(x3d))
    count = np.zeros(len(x3d))

    # list of dictionary containing border points and simplices
    dict_boundary = []
    
    dtfe = dtfe3d.Delaunay3D()

    # cycle through partition regions
    for i in range(0, len(idx)):
        
        pix_xcond = np.where((x3d[pixID] >= xedges[idx[i]]) & (x3d[pixID] < xedges[idx[i]+1]))[0]
        
        xcond = np.where((x >= xedges[idx[i]]) & (x < xedges[idx[i]+1]))[0]
        
        for j in range(0, len(idy)):
            
            pix_ycond = np.where((y3d[pixID[pix_xcond]] >= yedges[idy[j]]) & (y3d[pixID[pix_xcond]] < yedges[idy[j]+1]))[0]

            ycond = np.where((y[xcond] >= yedges[idy[j]]) & (y[xcond] < yedges[idy[j]+1]))[0]

            for k in range(0, len(idz)):
            
                pix_zcond = np.where((z3d[pixID[pix_xcond[pix_ycond]]] >= zedges[idz[k]]) & (z3d[pixID[pix_xcond[pix_ycond]]] < zedges[idz[k]+1]))[0]

                zcond = np.where((z[xcond[ycond]] >= zedges[idz[k]]) & (z[xcond[ycond]] < zedges[idz[k]+1]))[0]
            
                # local dtfe
                dtfe.clean()

                # partition boundary
                boundary = [xedges[idx[i]], xedges[idx[i]+1], yedges[idy[j]], yedges[idy[j]+1], zedges[idz[k]], zedges[idz[k]+1]]
                    
                dtfe.setup(x[xcond[ycond[zcond]]], y[xcond[ycond[zcond]]], z[xcond[ycond[zcond]]], boundary, mass=mass[xcond[ycond[zcond]]])

                dtfe.run()
                dtfe.get_dens()
                dtfe.set_field2dens()
                
                for _dx3d, _dy3d, _dz3d in zip(dx3d, dy3d, dz3d):
                    dtfe_estimate = dtfe.estimate(
                        x3d[pixID[pix_xcond[pix_ycond[pix_zcond]]]]+_dx3d, 
                        y3d[pixID[pix_xcond[pix_ycond[pix_zcond]]]]+_dy3d, 
                        z3d[pixID[pix_xcond[pix_ycond[pix_zcond]]]]+_dz3d)
                    cond = np.where(np.isfinite(dtfe_estimate) & (count[pix_xcond[pix_ycond[pix_zcond]]] != len(dx3d)))[0]
                    dens[pixID[pix_xcond[pix_ycond[pix_zcond[cond]]]]] += dtfe_estimate[cond]
                    count[pix_xcond[pix_ycond[pix_zcond[cond]]]] += 1.

                dict_boundary.append(dtfe.get_border())

        cond = np.where(count[pix_xcond] == len(dx3d))[0]
        pixID = np.delete(pixID, pix_xcond[cond])
        count = np.delete(count, pix_xcond[cond])

    dtfe = dtfe3d.DelaunayMerger3D()

    if len(idx) != 1 or len(idy) != 1 or len(idz) != 1:
        
        dict_boundary3 = []

        for i in range(0, len(idx)):
            
            dict_boundary2 = []

            for j in range(0, len(idy)):

                # local dtfe
                dtfe.clean()
                
                for k in range(0, len(idz)):
                    dtfe.add_border(dict_boundary[(i*len(idy) + j)*len(idz) + k])
            
                # partition boundary
                boundary = [xedges[idx[i]], xedges[idx[i]+1], yedges[idy[j]], yedges[idy[j]+1], zorigin, zorigin+zboxsize]

                dtfe.run(boundary)

                cond = np.where(
                    (x3d[pixID] >= boundary[0]) & (x3d[pixID] < boundary[1]) & 
                    (y3d[pixID] >= boundary[2]) & (y3d[pixID] < boundary[3]) & 
                    (z3d[pixID] >= boundary[4]) & (z3d[pixID] < boundary[5])
                )[0]

                for _dx3d, _dy3d, _dz3d in zip(dx3d, dy3d, dz3d):
                    dtfe_estimate = dtfe.estimate(x3d[pixID[cond]]+_dx3d, y3d[pixID[cond]]+_dy3d, z3d[pixID[cond]]+_dz3d)
                    cond_isfinite = np.where(np.isfinite(dtfe_estimate) & (count[cond] != len(dx3d)))[0]
                    dens[pixID[cond[cond_isfinite]]] += dtfe_estimate[cond_isfinite]
                    count[cond[cond_isfinite]] += 1.

                cond2 = np.where(count[cond] == len(dx3d))[0]
                pixID = np.delete(pixID, cond[cond2])
                count = np.delete(count, cond[cond2])

                dict_boundary2.append(dtfe.get_border())

            dtfe.clean()
            
            for boundary in dict_boundary2:
                dtfe.add_border(boundary)

            # partition boundary
            boundary = [xedges[idx[i]], xedges[idx[i]+1], yorigin, yorigin+yboxsize, zorigin, zorigin+zboxsize]

            dtfe.run(boundary)

            cond = np.where(
                (x3d[pixID] >= boundary[0]) & (x3d[pixID] < boundary[1]) & 
                (y3d[pixID] >= boundary[2]) & (y3d[pixID] < boundary[3]) & 
                (z3d[pixID] >= boundary[4]) & (z3d[pixID] < boundary[5])
            )[0]

            for _dx3d, _dy3d, _dz3d in zip(dx3d, dy3d, dz3d):
                dtfe_estimate = dtfe.estimate(x3d[pixID[cond]]+_dx3d, y3d[pixID[cond]]+_dy3d, z3d[pixID[cond]]+_dz3d)
                cond_isfinite = np.where(np.isfinite(dtfe_estimate) & (count[cond] != len(dx3d)))[0]
                dens[pixID[cond[cond_isfinite]]] += dtfe_estimate[cond_isfinite]
                count[cond[cond_isfinite]] += 1.

            cond2 = np.where(count[cond] == len(dx3d))[0]
            pixID = np.delete(pixID, cond[cond2])
            count = np.delete(count, cond[cond2])

            dict_boundary3.append(dtfe.get_border())
        
        dtfe.clean()

        for boundary in dict_boundary3:
            dtfe.add_border(boundary)

        # partition boundary
        boundary = [xorigin, xorigin+xboxsize, yorigin, yorigin+yboxsize, zorigin, zorigin+zboxsize]
        
        dtfe.run(boundary)

        cond = np.where(
            (x3d[pixID] >= boundary[0]) & (x3d[pixID] < boundary[1]) & 
            (y3d[pixID] >= boundary[2]) & (y3d[pixID] < boundary[3]) & 
            (z3d[pixID] >= boundary[4]) & (z3d[pixID] < boundary[5])
        )[0]

        for _dx3d, _dy3d, _dz3d in zip(dx3d, dy3d, dz3d):
            dtfe_estimate = dtfe.estimate(x3d[pixID[cond]]+_dx3d, y3d[pixID[cond]]+_dy3d, z3d[pixID[cond]]+_dz3d)
            cond_isfinite = np.where(np.isfinite(dtfe_estimate) & (count[cond] != len(dx3d)))[0]
            dens[pixID[cond[cond_isfinite]]] += dtfe_estimate[cond_isfinite]
            count[cond[cond_isfinite]] += 1.

        cond2 = np.where(count[cond] == len(dx3d))[0]
        pixID = np.delete(pixID, cond[cond2])
        count = np.delete(count, cond[cond2])

        exterior_border = dtfe.get_border()

    else:
        exterior_border = dict_boundary[0]
    
    _exterior_border = []

    if xperiodic == True or yperiodic == True or zperiodic == True:

        xbuffer = fbuffer*xboxsize
        ybuffer = fbuffer*yboxsize
        zbuffer = fbuffer*zboxsize

        if xperiodic == True:
            imin, imax = -1, 2
        else:
            imin, imax = 0, 1
        if yperiodic == True:
            jmin, jmax = -1, 2
        else:
            jmin, jmax = 0, 1
        if zperiodic == True:
            kmin, kmax = -1, 2
        else:
            kmin, kmax = 0, 1
        
        for i in range(imin, imax):
            for j in range(jmin, jmax):
                for k in range(kmin, kmax):
                    _border = deepcopy(exterior_border)
                    _border['x'] += i*xboxsize
                    _border['y'] += j*yboxsize
                    _border['z'] += k*zboxsize
                    _exterior_border.append(_border)
    
        dtfe.clean()

        for border in _exterior_border:
            dtfe.add_border(border)

        boundary = [xorigin-xbuffer, xorigin+xboxsize+xbuffer, yorigin-ybuffer, yorigin+yboxsize+ybuffer, zorigin-zbuffer, zorigin+zboxsize+zbuffer]

        dtfe.run(boundary, apply_filter=True)

        cond = np.where(
            (x3d[pixID] >= boundary[0]) & (x3d[pixID] < boundary[1]) & 
            (y3d[pixID] >= boundary[2]) & (y3d[pixID] < boundary[3]) & 
            (z3d[pixID] >= boundary[4]) & (z3d[pixID] < boundary[5])
        )[0]

        for _dx3d, _dy3d, _dz3d in zip(dx3d, dy3d, dz3d):
            dtfe_estimate = dtfe.estimate(x3d[pixID[cond]]+_dx3d, y3d[pixID[cond]]+_dy3d, z3d[pixID[cond]]+_dz3d)
            cond_isfinite = np.where(np.isfinite(dtfe_estimate) & (count[cond] != len(dx3d)))[0]
            dens[pixID[cond[cond_isfinite]]] += dtfe_estimate[cond_isfinite]
            count[cond[cond_isfinite]] += 1.

        cond2 = np.where(count[cond] == len(dx3d))[0]
        pixID = np.delete(pixID, cond[cond2])
        count = np.delete(count, cond[cond2])

        exterior_border = dtfe.get_border()

    if normalise:
        dens /= len(dx3d)
        cond = np.where(count == 0.)[0]
        dens[pixID[cond]] *= len(dx3d)
        cond = np.where(count != 0.)[0]
        dens[pixID[cond]] *= len(dx3d)/count[cond]

    if flatten == False:
        dens = dens.reshape(nxgrid,nygrid,nzgrid)
    
    if outputexterior:
        if outputgrid:
            return dens, x3d, y3d, z3d, exterior_border, pixID, count
        else:
            return dens, exterior_border, pixID, count
    else:
        if outputgrid:
            return dens, x3d, y3d, z3d
        else:
            return dens



def delaunay_field4grid3D(
    x: np.ndarray, y: np.ndarray, z: np.ndarray, f: np.ndarray, boxsize: Union[float, List[float]], ngrid: Union[int, List[int]], 
    origin: Union[float, List[float]] = 0., mass: Optional[np.ndarray] = None, partition: int = 1, 
    periodic : Union[bool, List[bool]] = True, fbuffer: float = 0.2, subsampling: int = 1, 
    outputgrid: bool = False, outputexterior: bool = False, normalise: bool = True, flatten: bool = False
) -> Union[
    np.ndarray, 
    Tuple[np.ndarray, np.ndarray, np.ndarray],
    Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray],
]:
    """
    Returns the Delaunay tesselation field on a grid.

    Parameters
    ----------
    x, y, z : array
        Coordinates of particles.
    f : array
        Field values to be interpolated.
    boxsize : float or list
        Dimensions of the grid.
    ngrid : int or int list
        Grid dimensions.
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
    normalise : bool, optional
        Whether to normalise the outputing field.
    flatten : bool, optional
        Flatten output field array.

    Returns
    -------
    field : ndarray
        Field values on a grid.
    x2d, y2d, z3d : ndarray, optional
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
    
    # define grid coordinates
    x3d, y3d, z3d = shift.cart.grid3D(boxsize, ngrid, origin=[xorigin, yorigin, zorigin])
    x3d, y3d, z3d = x3d.flatten(), y3d.flatten(), z3d.flatten()
    field = np.zeros(len(x3d))

    # define pixel length across each axis
    dx = xboxsize/nxgrid
    dy = yboxsize/nygrid
    dz = zboxsize/nzgrid

    if np.isscalar(subsampling):
        dnxgrid, dnygrid, dnzgrid = subsampling, subsampling, subsampling
    else:
        dnxgrid, dnygrid, dnzgrid = subsampling[0], subsampling[1], subsampling[2]

    # define subsampling points for each pixel
    dx3d, dy3d, dz3d = shift.cart.grid3D([dx, dy, dz], [dnxgrid, dnygrid, dnzgrid], origin=[-dx/2., -dy/2., -dz/2.])
    dx3d, dy3d, dz3d = dx3d.flatten(), dy3d.flatten(), dz3d.flatten()

    # partition data
    if np.isscalar(partition):
        xpartition, ypartition, zpartition = partition, partition, partition
    else:
        xpartition, ypartition, zpartition = partition[0], partition[1], partition[2]
    
    # construct partition divisions
    xedges, _ = shift.cart.grid1D(xboxsize, xpartition, origin=xorigin)
    yedges, _ = shift.cart.grid1D(yboxsize, ypartition, origin=yorigin)
    zedges, _ = shift.cart.grid1D(zboxsize, zpartition, origin=zorigin)

    xedges[:-1] = np.floor(xedges[:-1] / dx) * dx
    yedges[:-1] = np.floor(yedges[:-1] / dy) * dy
    zedges[:-1] = np.floor(zedges[:-1] / dz) * dz
    
    # construct partition indexes
    idx = np.arange(xpartition)
    idy = np.arange(ypartition)
    idz = np.arange(zpartition)

    # assign mass weights
    if mass is None:
        mass = np.ones(len(x))

    pixID = np.arange(len(x3d))
    count = np.zeros(len(x3d))

    # list of dictionary containing border points and simplices
    dict_boundary = []
    
    dtfe = dtfe3d.Delaunay3D()

    # cycle through partition regions
    for i in range(0, len(idx)):
        
        pix_xcond = np.where((x3d[pixID] >= xedges[idx[i]]) & (x3d[pixID] < xedges[idx[i]+1]))[0]
        
        xcond = np.where((x >= xedges[idx[i]]) & (x < xedges[idx[i]+1]))[0]
        
        for j in range(0, len(idy)):
            
            pix_ycond = np.where((y3d[pixID[pix_xcond]] >= yedges[idy[j]]) & (y3d[pixID[pix_xcond]] < yedges[idy[j]+1]))[0]

            ycond = np.where((y[xcond] >= yedges[idy[j]]) & (y[xcond] < yedges[idy[j]+1]))[0]

            for k in range(0, len(idz)):
            
                pix_zcond = np.where((z3d[pixID[pix_xcond[pix_ycond]]] >= zedges[idz[k]]) & (z3d[pixID[pix_xcond[pix_ycond]]] < zedges[idz[k]+1]))[0]

                zcond = np.where((z[xcond[ycond]] >= zedges[idz[k]]) & (z[xcond[ycond]] < zedges[idz[k]+1]))[0]
            
                # local dtfe
                dtfe.clean()

                # partition boundary
                boundary = [xedges[idx[i]], xedges[idx[i]+1], yedges[idy[j]], yedges[idy[j]+1], zedges[idz[k]], zedges[idz[k]+1]]
                    
                dtfe.setup(x[xcond[ycond[zcond]]], y[xcond[ycond[zcond]]], z[xcond[ycond[zcond]]], boundary, mass=mass[xcond[ycond[zcond]]])

                dtfe.run()
                dtfe.set_field(f[xcond[ycond[zcond]]])
                
                for _dx3d, _dy3d, _dz3d in zip(dx3d, dy3d, dz3d):
                    dtfe_estimate = dtfe.estimate(
                        x3d[pixID[pix_xcond[pix_ycond[pix_zcond]]]]+_dx3d, 
                        y3d[pixID[pix_xcond[pix_ycond[pix_zcond]]]]+_dy3d, 
                        z3d[pixID[pix_xcond[pix_ycond[pix_zcond]]]]+_dz3d)
                    cond = np.where(np.isfinite(dtfe_estimate) & (count[pix_xcond[pix_ycond[pix_zcond]]] != len(dx3d)))[0]
                    field[pixID[pix_xcond[pix_ycond[pix_zcond[cond]]]]] += dtfe_estimate[cond]
                    count[pix_xcond[pix_ycond[pix_zcond[cond]]]] += 1.
                
                dict_boundary.append(dtfe.get_border())

        cond = np.where(count[pix_xcond] == len(dx3d))[0]
        pixID = np.delete(pixID, pix_xcond[cond])
        count = np.delete(count, pix_xcond[cond])

    dtfe = dtfe3d.DelaunayMerger3D()

    if len(idx) != 1 or len(idy) != 1 or len(idz) != 1:
        
        dict_boundary3 = []

        for i in range(0, len(idx)):
            
            dict_boundary2 = []

            for j in range(0, len(idy)):

                # local dtfe
                dtfe.clean()
                
                for k in range(0, len(idz)):
                    dtfe.add_border(dict_boundary[(i*len(idy) + j)*len(idz) + k])
            
                # partition boundary
                boundary = [xedges[idx[i]], xedges[idx[i]+1], yedges[idy[j]], yedges[idy[j]+1], zorigin, zorigin+zboxsize]

                dtfe.run(boundary)

                cond = np.where(
                    (x3d[pixID] >= boundary[0]) & (x3d[pixID] < boundary[1]) & 
                    (y3d[pixID] >= boundary[2]) & (y3d[pixID] < boundary[3]) & 
                    (z3d[pixID] >= boundary[4]) & (z3d[pixID] < boundary[5])
                )[0]

                for _dx3d, _dy3d, _dz3d in zip(dx3d, dy3d, dz3d):
                    dtfe_estimate = dtfe.estimate(x3d[pixID[cond]]+_dx3d, y3d[pixID[cond]]+_dy3d, z3d[pixID[cond]]+_dz3d)
                    cond_isfinite = np.where(np.isfinite(dtfe_estimate))[0]
                    field[pixID[cond[cond_isfinite]]] += dtfe_estimate[cond_isfinite]
                    count[cond[cond_isfinite]] += 1.

                cond2 = np.where(count[cond] == len(dx3d))[0]

                pixID = np.delete(pixID, cond[cond2])
                count = np.delete(count, cond[cond2])

                dict_boundary2.append(dtfe.get_border())

            dtfe.clean()
            
            for boundary in dict_boundary2:
                dtfe.add_border(boundary)

            # partition boundary
            boundary = [xedges[idx[i]], xedges[idx[i]+1], yorigin, yorigin+yboxsize, zorigin, zorigin+zboxsize]

            dtfe.run(boundary)

            cond = np.where(
                (x3d[pixID] >= boundary[0]) & (x3d[pixID] < boundary[1]) & 
                (y3d[pixID] >= boundary[2]) & (y3d[pixID] < boundary[3]) & 
                (z3d[pixID] >= boundary[4]) & (z3d[pixID] < boundary[5])
            )[0]

            for _dx3d, _dy3d, _dz3d in zip(dx3d, dy3d, dz3d):
                dtfe_estimate = dtfe.estimate(x3d[pixID[cond]]+_dx3d, y3d[pixID[cond]]+_dy3d, z3d[pixID[cond]]+_dz3d)
                cond_isfinite = np.where(np.isfinite(dtfe_estimate) & (count[cond] != len(dx3d)))[0]
                field[pixID[cond[cond_isfinite]]] += dtfe_estimate[cond_isfinite]
                count[cond[cond_isfinite]] += 1.

            cond2 = np.where(count[cond] == len(dx3d))[0]

            pixID = np.delete(pixID, cond[cond2])
            count = np.delete(count, cond[cond2])

            dict_boundary3.append(dtfe.get_border())
        
        dtfe.clean()
        
        for boundary in dict_boundary3:
            dtfe.add_border(boundary)

        # partition boundary
        boundary = [xorigin, xorigin+xboxsize, yorigin, yorigin+yboxsize, zorigin, zorigin+zboxsize]
        
        dtfe.run(boundary)

        cond = np.where(
            (x3d[pixID] >= boundary[0]) & (x3d[pixID] < boundary[1]) & 
            (y3d[pixID] >= boundary[2]) & (y3d[pixID] < boundary[3]) & 
            (z3d[pixID] >= boundary[4]) & (z3d[pixID] < boundary[5])
        )[0]

        for _dx3d, _dy3d, _dz3d in zip(dx3d, dy3d, dz3d):
            dtfe_estimate = dtfe.estimate(x3d[pixID[cond]]+_dx3d, y3d[pixID[cond]]+_dy3d, z3d[pixID[cond]]+_dz3d)
            cond_isfinite = np.where(np.isfinite(dtfe_estimate) & (count[cond] != len(dx3d)))[0]
            field[pixID[cond[cond_isfinite]]] += dtfe_estimate[cond_isfinite]
            count[cond[cond_isfinite]] += 1.

        cond2 = np.where(count[cond] == len(dx3d))[0]

        pixID = np.delete(pixID, cond[cond2])
        count = np.delete(count, cond[cond2])

        exterior_border = dtfe.get_border()
    else:
        exterior_border = dict_boundary[0]

    _exterior_border = []

    if xperiodic == True or yperiodic == True or zperiodic == True:

        xbuffer = fbuffer*xboxsize
        ybuffer = fbuffer*yboxsize
        zbuffer = fbuffer*zboxsize

        if xperiodic == True:
            imin, imax = -1, 2
        else:
            imin, imax = 0, 1
        if yperiodic == True:
            jmin, jmax = -1, 2
        else:
            jmin, jmax = 0, 1
        if zperiodic == True:
            kmin, kmax = -1, 2
        else:
            kmin, kmax = 0, 1
        
        for i in range(imin, imax):
            for j in range(jmin, jmax):
                for k in range(kmin, kmax):
                    _border = deepcopy(exterior_border)
                    _border['x'] += i*xboxsize
                    _border['y'] += j*yboxsize
                    _border['z'] += k*zboxsize
                    _exterior_border.append(_border)
                    
        dtfe.clean()

        for border in _exterior_border:
            dtfe.add_border(border)

        boundary = [xorigin-xbuffer, xorigin+xboxsize+xbuffer, yorigin-ybuffer, yorigin+yboxsize+ybuffer, zorigin-zbuffer, zorigin+zboxsize+zbuffer]

        dtfe.run(boundary, apply_filter=True)

        cond = np.where(
            (x3d[pixID] >= boundary[0]) & (x3d[pixID] < boundary[1]) & 
            (y3d[pixID] >= boundary[2]) & (y3d[pixID] < boundary[3]) & 
            (z3d[pixID] >= boundary[4]) & (z3d[pixID] < boundary[5])
        )[0]

        for _dx3d, _dy3d, _dz3d in zip(dx3d, dy3d, dz3d):
            dtfe_estimate = dtfe.estimate(x3d[pixID[cond]]+_dx3d, y3d[pixID[cond]]+_dy3d, z3d[pixID[cond]]+_dz3d)
            cond_isfinite = np.where(np.isfinite(dtfe_estimate) & (count[cond] != len(dx3d)))[0]
            field[pixID[cond[cond_isfinite]]] += dtfe_estimate[cond_isfinite]
            count[cond[cond_isfinite]] += 1.

        cond2 = np.where(count[cond] == len(dx3d))[0]
        pixID = np.delete(pixID, cond[cond2])
        count = np.delete(count, cond[cond2])

        exterior_border = dtfe.get_border()

    if normalise:
        field /= len(dx3d)
        cond = np.where(count == 0.)[0]
        field[pixID[cond]] *= len(dx3d)
        cond = np.where(count != 0.)[0]
        field[pixID[cond]] *= len(dx3d)/count[cond]

    if flatten == False:
        field = field.reshape(nxgrid,nygrid,nzgrid)

    if outputexterior:
        if outputgrid:
            return field, x3d, y3d, z3d, exterior_border, pixID, count
        else:
            return field, exterior_border, pixID, count
    else:
        if outputgrid:
            return field, x3d, y3d, z3d
        else:
            return field