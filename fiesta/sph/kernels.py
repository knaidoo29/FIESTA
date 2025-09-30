import numpy as np

from typing import Union 


def cubic_kernel(r: Union[float, np.ndarray], h: Union[float, np.ndarray], dim: int = 3) -> Union[float, np.ndarray]:
    """
    Cubic spline kernel.

    Parameters
    ----------
    r : float or array
        Distance to center point.
    h : float or array
        Kernel size.

    Returns
    -------
    w : float or array
        Kernel weights.
    """
    if dim == 1:
        sig3 = 2./(3.*h)
    elif dim == 2:
        sig3 = 10./(7*np.pi*h**2.)
    elif dim == 3:
        sig3 = 1./(np.pi*h**3.)
    q = r/h
    if np.isscalar(r) is True:
        if q >= 0. and q <= 1.:
            w = (1. - 1.5*(q**2.)*(1.-q/2.))
        elif q <= 2.:
            w = 0.25*(2.-q)**3.
        else:
            w = 0.
    else:
        w = np.zeros(len(q))
        cond = np.where((q >= 0.) & (q <= 1.))[0]
        w[cond] = (1. - 1.5*(q[cond]**2.)*(1.-q[cond]/2.))
        cond = np.where((q > 1.) & (q <= 2.))[0]
        w[cond] = 0.25*(2.-q[cond])**3.
    w *= sig3
    return w


def dcubic_kernel(r: Union[float, np.ndarray], h: Union[float, np.ndarray], dim: int = 3) -> Union[float, np.ndarray]:
    """
    Cubic derivative spline kernel.

    Parameters
    ----------
    r : float or array
        Distance to center point.
    h : float or array
        Kernel size.

    Returns
    -------
    w : float or array
        Kernel weights.
    """
    if dim == 1:
        sig3 = 2./(3.*h)
    elif dim == 2:
        sig3 = 10./(7*np.pi*h**2.)
    elif dim == 3:
        sig3 = 1./(np.pi*h**3.)
    q = r/h
    if np.isscalar(r) is True:
        if q >= 0. and q <= 1.:
            w = - 3.*q + 2.25*q**2.
        elif q <= 2.:
            w = -0.75*(2.-q)**2.
        else:
            w = 0.
    else:
        w = np.zeros(len(q))
        cond = np.where((q >= 0.) & (q <= 1.))[0]
        w[cond] = -3*q[cond] + 2.25*q[cond]**2.
        cond = np.where((q > 1.) & (q <= 2.))[0]
        w[cond] = -0.75*(2.-q[cond])**2.
    w *= sig3
    return w
