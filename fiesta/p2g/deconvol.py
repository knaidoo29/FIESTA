import numpy as np
import shift


def get_sinc(x):
    """Returns the sinc function for input x."""
    if np.isscalar(x) is True:
        if x == 0:
            sinc = 1.
        else:
            sinc = np.sin(x)/x
    else:
        sinc = np.ones(np.shape(x))
        condition = np.where(x != 0.)
        sinc[condition] = np.sin(x[condition])/x[condition]
    return sinc


def get_deconvol_p(method):
    """Returns the deconvolution power.

    Parameters
    ----------
    method : str
        grid assignment scheme, either NGP, CIC or TSC.

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


def deconvolve_part2grid_2D(field, boxsize, method='TSC'):
    """Deconvolve the grid assignment scheme in Fourier space.

    Parameters
    ----------
    field : ndarray
        Grid assigned field.
    boxsize : float
        Box size.
    method : str
        grid assignment scheme, either NGP, CIC, TSC or PCS.
    """
    fieldk = shift.cart.fft2D(field, boxsize)
    ngrid = len(fieldk)
    kx3d, ky3d = shift.cart.kgrid2D(boxsize, ngrid)
    kn = shift.cart.get_kn(boxsize, ngrid)
    if method == 'NGP':
        C = 1.
    elif method == 'CIC':
        C  = 1.-(2./3.)*(np.sin(np.pi*kx3d/(2*kn)))**2.
        C *= 1.-(2./3.)*(np.sin(np.pi*ky3d/(2*kn)))**2.
    elif method == 'TSC':
        C  = 1.-(np.sin(np.pi*kx3d/(2*kn)))**2. + (2./15.)*(np.sin(np.pi*kx3d/(2*kn)))**4.
        C *= 1.-(np.sin(np.pi*ky3d/(2*kn)))**2. + (2./15.)*(np.sin(np.pi*ky3d/(2*kn)))**4.
    elif method == 'PCS':
        C  = 1.-(4./3.)*(np.sin(np.pi*kx3d/(2*kn)))**2. + (2./5.)*(np.sin(np.pi*kx3d/(2*kn)))**4. - (4./315.)*(np.sin(np.pi*kx3d/(2*kn)))**6.
        C *= 1.-(4./3.)*(np.sin(np.pi*ky3d/(2*kn)))**2. + (2./5.)*(np.sin(np.pi*ky3d/(2*kn)))**4. - (4./315.)*(np.sin(np.pi*ky3d/(2*kn)))**6.
    fieldk /= np.sqrt(C)
    field = shift.cart.ifft2D(fieldk, boxsize)
    return field


def deconvolve_part2grid_3D(field, boxsize, method='TSC'):
    """Deconvolve the grid assignment scheme in Fourier space.

    Parameters
    ----------
    field : ndarray
        Grid assigned field.
    boxsize : float
        Box size.
    method : str
        grid assignment scheme, either NGP, CIC, TSC or PCS.
    """
    fieldk = shift.cart.fft3D(field, boxsize)
    ngrid = len(fieldk)
    kx3d, ky3d, kz3d = shift.cart.kgrid3D(boxsize, ngrid)
    kn = shift.cart.get_kn(boxsize, ngrid)
    if method == 'NGP':
        C = 1.
    elif method == 'CIC':
        C  = 1.-(2./3.)*(np.sin(np.pi*kx3d/(2*kn)))**2.
        C *= 1.-(2./3.)*(np.sin(np.pi*ky3d/(2*kn)))**2.
        C *= 1.-(2./3.)*(np.sin(np.pi*kz3d/(2*kn)))**2.
    elif method == 'TSC':
        C  = 1.-(np.sin(np.pi*kx3d/(2*kn)))**2. + (2./15.)*(np.sin(np.pi*kx3d/(2*kn)))**4.
        C *= 1.-(np.sin(np.pi*ky3d/(2*kn)))**2. + (2./15.)*(np.sin(np.pi*ky3d/(2*kn)))**4.
        C *= 1.-(np.sin(np.pi*kz3d/(2*kn)))**2. + (2./15.)*(np.sin(np.pi*kz3d/(2*kn)))**4.
    elif method == 'PCS':
        C  = 1.-(4./3.)*(np.sin(np.pi*kx3d/(2*kn)))**2. + (2./5.)*(np.sin(np.pi*kx3d/(2*kn)))**4. - (4./315.)*(np.sin(np.pi*kx3d/(2*kn)))**6.
        C *= 1.-(4./3.)*(np.sin(np.pi*ky3d/(2*kn)))**2. + (2./5.)*(np.sin(np.pi*ky3d/(2*kn)))**4. - (4./315.)*(np.sin(np.pi*ky3d/(2*kn)))**6.
        C *= 1.-(4./3.)*(np.sin(np.pi*kz3d/(2*kn)))**2. + (2./5.)*(np.sin(np.pi*kz3d/(2*kn)))**4. - (4./315.)*(np.sin(np.pi*kz3d/(2*kn)))**6.
    fieldk /= np.sqrt(C)
    field = shift.cart.ifft3D(fieldk, boxsize)
    return field
