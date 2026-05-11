import numpy as np
import shift


from typing import List, Union


def get_deconvol_p(method: str) -> int:
    """
    Returns the deconvolution power.

    Parameters
    ----------
    method : str
        grid assignment scheme, either NGP, CIC, TSC or PCS.

    Returns
    -------
    p : float
        Deconvolution power.
    """
    if method == 'NGP':
        p = 1.
    elif method == 'CIC':
        p = 2.
    elif method == 'TSC':
        p = 3.
    elif method == 'PCS':
        p = 4.
    return p

# TODO: removed this or incorporate it for what it is...

# Note: this is not deconvolution but rather the removal of the window function C(k)
# from https://arxiv.org/pdf/2403.13561v1, these functions should be revised to adjust
# for this.

# def deconvolve_part2grid_2D(field, boxsize, method='TSC'):
#     """Deconvolve the grid assignment scheme in Fourier space.
#
#     Parameters
#     ----------
#     field : ndarray
#         Grid assigned field.
#     boxsize : float
#         Box size.
#     method : str
#         grid assignment scheme, either NGP, CIC, TSC or PCS.
#     """
#     fieldk = shift.cart.fft2D(field, boxsize)
#     ngrid = len(fieldk)
#     kx3d, ky3d = shift.cart.kgrid2D(boxsize, ngrid)
#     kn = shift.cart.get_kn(boxsize, ngrid)
#     if method == 'NGP':
#         C = 1.
#     elif method == 'CIC':
#         C  = 1.-(2./3.)*(np.sin(np.pi*kx3d/(2*kn)))**2.
#         C *= 1.-(2./3.)*(np.sin(np.pi*ky3d/(2*kn)))**2.
#     elif method == 'TSC':
#         C  = 1.-(np.sin(np.pi*kx3d/(2*kn)))**2. + (2./15.)*(np.sin(np.pi*kx3d/(2*kn)))**4.
#         C *= 1.-(np.sin(np.pi*ky3d/(2*kn)))**2. + (2./15.)*(np.sin(np.pi*ky3d/(2*kn)))**4.
#     elif method == 'PCS':
#         C  = 1.-(4./3.)*(np.sin(np.pi*kx3d/(2*kn)))**2. + (2./5.)*(np.sin(np.pi*kx3d/(2*kn)))**4. - (4./315.)*(np.sin(np.pi*kx3d/(2*kn)))**6.
#         C *= 1.-(4./3.)*(np.sin(np.pi*ky3d/(2*kn)))**2. + (2./5.)*(np.sin(np.pi*ky3d/(2*kn)))**4. - (4./315.)*(np.sin(np.pi*ky3d/(2*kn)))**6.
#     fieldk /= np.sqrt(C)
#     field = shift.cart.ifft2D(fieldk, boxsize)
#     return field


# def deconvolve_part2grid_3D(field, boxsize, method='TSC'):
#     """Deconvolve the grid assignment scheme in Fourier space.
#
#     Parameters
#     ----------
#     field : ndarray
#         Grid assigned field.
#     boxsize : float
#         Box size.
#     method : str
#         grid assignment scheme, either NGP, CIC, TSC or PCS.
#     """
#     fieldk = shift.cart.fft3D(field, boxsize)
#     ngrid = len(fieldk)
#     kx3d, ky3d, kz3d = shift.cart.kgrid3D(boxsize, ngrid)
#     kn = shift.cart.get_kn(boxsize, ngrid)
#     if method == 'NGP':
#         C = 1.
#     elif method == 'CIC':
#         C  = 1.-(2./3.)*(np.sin(np.pi*kx3d/(2*kn)))**2.
#         C *= 1.-(2./3.)*(np.sin(np.pi*ky3d/(2*kn)))**2.
#         C *= 1.-(2./3.)*(np.sin(np.pi*kz3d/(2*kn)))**2.
#     elif method == 'TSC':
#         C  = 1.-(np.sin(np.pi*kx3d/(2*kn)))**2. + (2./15.)*(np.sin(np.pi*kx3d/(2*kn)))**4.
#         C *= 1.-(np.sin(np.pi*ky3d/(2*kn)))**2. + (2./15.)*(np.sin(np.pi*ky3d/(2*kn)))**4.
#         C *= 1.-(np.sin(np.pi*kz3d/(2*kn)))**2. + (2./15.)*(np.sin(np.pi*kz3d/(2*kn)))**4.
#     elif method == 'PCS':
#         C  = 1.-(4./3.)*(np.sin(np.pi*kx3d/(2*kn)))**2. + (2./5.)*(np.sin(np.pi*kx3d/(2*kn)))**4. - (4./315.)*(np.sin(np.pi*kx3d/(2*kn)))**6.
#         C *= 1.-(4./3.)*(np.sin(np.pi*ky3d/(2*kn)))**2. + (2./5.)*(np.sin(np.pi*ky3d/(2*kn)))**4. - (4./315.)*(np.sin(np.pi*ky3d/(2*kn)))**6.
#         C *= 1.-(4./3.)*(np.sin(np.pi*kz3d/(2*kn)))**2. + (2./5.)*(np.sin(np.pi*kz3d/(2*kn)))**4. - (4./315.)*(np.sin(np.pi*kz3d/(2*kn)))**6.
#     fieldk /= np.sqrt(C)
#     field = shift.cart.ifft3D(fieldk, boxsize)
#     return field


def get_part2grid2D_kernel(
        kx2d: Union[float, np.ndarray], ky2d: Union[float, np.ndarray], ngrid: Union[int, List[int]], boxsize: Union[float, List[float]], 
        method: str = 'TSC'
    ) -> np.ndarray:
    """
    Returns the convolution kernel for particle-to-grid assignment scheme on a Fourier space grid.

    Parameters
    ----------
    kx2d, ky2d : float or array
        Fourier coefficients along the x and y axis.
    ngrid: int or list
        Grid dimensions.
    boxsize : float or list
        Box size. 
    method : str
        grid assignment scheme, either NGP, CIC, TSC or PCS.
    
    Returns
    -------
    Wk : array
        The part-to-grid mass assignment kernel in Fourier space.
    """
    if np.isscalar(ngrid):
        nxgrid, nygrid = ngrid, ngrid
    else:
        nxgrid, nygrid = ngrid[0], ngrid[1]
    if np.isscalar(boxsize):
        xboxsize, yboxsize = boxsize, boxsize
    else:
        xboxsize, yboxsize = boxsize[0], boxsize[1]
    # knx = shift.cart.get_kn(xboxsize, nxgrid)
    # kny = shift.cart.get_kn(yboxsize, nygrid)
    # Wk =  np.sinc(kx2d/(np.pi*2*knx))
    # Wk *= np.sinc(ky2d/(np.pi*2*kny))
    # Wk =  np.sinc(kx2d/(2*knx))
    # Wk *= np.sinc(ky2d/(2*kny))
    dx = xboxsize/nxgrid
    dy = yboxsize/nygrid
    # dz = zboxsize/nzgrid
    Wk = np.sinc(kx2d*dx/(2*np.pi))**get_deconvol_p(method)
    Wk *= np.sinc(ky2d*dy/(2*np.pi))**get_deconvol_p(method)
    # Wk *= np.sinc(kz3d*dz/(2*np.pi))
    # Wk = Wk**get_deconvol_p(method)
    return Wk


def get_part2grid3D_kernel(
        kx3d: Union[float, np.ndarray], ky3d: Union[float, np.ndarray], kz3d: Union[float, np.ndarray], ngrid: Union[int, List[int]], 
        boxsize: Union[float, List[float]], method: str = 'TSC'
    ) -> np.ndarray:
    """
    Returns the convolution kernel for particle-to-grid assignment scheme on a Fourier space grid.

    Parameters
    ----------
    kx3d, ky3d, kz3d : float or array
        Fourier coefficients along the x, y and z axis.
    ngrid: int or list
        Grid dimensions.
    boxsize : float or list
        Box size. 
    method : str
        grid assignment scheme, either NGP, CIC, TSC or PCS.
    
    Returns
    -------
    Wk : array
        The part-to-grid mass assignment kernel in Fourier space.
    """
    if np.isscalar(ngrid):
        nxgrid, nygrid, nzgrid = ngrid, ngrid, ngrid
    else:
        nxgrid, nygrid, nzgrid = ngrid[0], ngrid[1], ngrid[2]
    if np.isscalar(boxsize):
        xboxsize, yboxsize, zboxsize = boxsize, boxsize, boxsize
    else:
        xboxsize, yboxsize, zboxsize = boxsize[0], boxsize[1], boxsize[2]
    # knx = shift.cart.get_kn(xboxsize, nxgrid)
    # kny = shift.cart.get_kn(yboxsize, nygrid)
    # knz = shift.cart.get_kn(zboxsize, nzgrid)
    # Wk =  np.sinc(kx3d/(np.pi*2*knx))
    # Wk *= np.sinc(ky3d/(np.pi*2*kny))
    # Wk *= np.sinc(kz3d/(np.pi*2*knz))
    # Wk =  np.sinc(kx3d/(2*knx))
    # Wk *= np.sinc(ky3d/(2*kny))
    # Wk *= np.sinc(kz3d/(2*knz))
    dx = xboxsize/nxgrid
    dy = yboxsize/nygrid
    dz = zboxsize/nzgrid
    Wk = np.sinc(kx3d*dx/(2*np.pi))**get_deconvol_p(method)
    Wk *= np.sinc(ky3d*dy/(2*np.pi))**get_deconvol_p(method)
    Wk *= np.sinc(kz3d*dz/(2*np.pi))**get_deconvol_p(method)
    # Wk = Wk**get_deconvol_p(method)
    return Wk



def deconvolve_part2grid_2D(field: np.ndarray, boxsize: Union[float, List[float]], method: str = 'TSC') -> np.ndarray:
    """
    Deconvolve the grid assignment scheme in Fourier space.

    Parameters
    ----------
    field : ndarray
        Grid assigned field.
    boxsize : float or list
        Box size.
    method : str
        grid assignment scheme, either NGP, CIC, TSC or PCS.
    """
    fieldk = shift.cart.fft2D(field, boxsize)
    ngrid = np.shape(fieldk)
    kx2d, ky2d = shift.cart.kgrid2D(boxsize, ngrid)
    Wk = get_part2grid2D_kernel(kx2d, ky2d, ngrid, boxsize, method)
    fieldk /= Wk
    field = shift.cart.ifft2D(fieldk, boxsize)
    return field


def deconvolve_part2grid_3D(field: np.ndarray, boxsize: Union[float, List[float]], method: str = 'TSC') -> np.ndarray:
    """
    Deconvolve the grid assignment scheme in Fourier space.

    Parameters
    ----------
    field : ndarray
        Grid assigned field.
    boxsize : float or list
        Box size.
    method : str
        grid assignment scheme, either NGP, CIC, TSC or PCS.
    """
    fieldk = shift.cart.fft3D(field, boxsize)
    ngrid = np.shape(fieldk)
    kx3d, ky3d, kz3d = shift.cart.kgrid3D(boxsize, ngrid)
    Wk = get_part2grid3D_kernel(kx3d, ky3d, kz3d, ngrid, boxsize, method)
    fieldk /= Wk
    field = shift.cart.ifft3D(fieldk, boxsize)
    return field
