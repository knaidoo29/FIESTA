import numpy as np

from typing import Union


def cubic_kernel(
    r: Union[float, np.ndarray], h: Union[float, np.ndarray], dim: int = 3
) -> Union[float, np.ndarray]:
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
        sig3 = 2.0 / (3.0 * h)
    elif dim == 2:
        sig3 = 10.0 / (7 * np.pi * h**2.0)
    elif dim == 3:
        sig3 = 1.0 / (np.pi * h**3.0)
    q = r / h
    if np.isscalar(r) is True:
        if q >= 0.0 and q <= 1.0:
            w = 1.0 - 1.5 * (q**2.0) * (1.0 - q / 2.0)
        elif q <= 2.0:
            w = 0.25 * (2.0 - q) ** 3.0
        else:
            w = 0.0
    else:
        w = np.zeros(len(q))
        cond = np.where((q >= 0.0) & (q <= 1.0))[0]
        w[cond] = 1.0 - 1.5 * (q[cond] ** 2.0) * (1.0 - q[cond] / 2.0)
        cond = np.where((q > 1.0) & (q <= 2.0))[0]
        w[cond] = 0.25 * (2.0 - q[cond]) ** 3.0
    w *= sig3
    return w


def dcubic_kernel(
    r: Union[float, np.ndarray], h: Union[float, np.ndarray], dim: int = 3
) -> Union[float, np.ndarray]:
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
        sig3 = 2.0 / (3.0 * h)
    elif dim == 2:
        sig3 = 10.0 / (7 * np.pi * h**2.0)
    elif dim == 3:
        sig3 = 1.0 / (np.pi * h**3.0)
    q = r / h
    if np.isscalar(r) is True:
        if q >= 0.0 and q <= 1.0:
            w = -3.0 * q + 2.25 * q**2.0
        elif q <= 2.0:
            w = -0.75 * (2.0 - q) ** 2.0
        else:
            w = 0.0
    else:
        w = np.zeros(len(q))
        cond = np.where((q >= 0.0) & (q <= 1.0))[0]
        w[cond] = -3 * q[cond] + 2.25 * q[cond] ** 2.0
        cond = np.where((q > 1.0) & (q <= 2.0))[0]
        w[cond] = -0.75 * (2.0 - q[cond]) ** 2.0
    w *= sig3
    return w
