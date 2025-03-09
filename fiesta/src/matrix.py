import numpy as np
from numba import njit


@njit
def inv2by2(m):
    """
    Invert a 2x2 matrix.

    Parameters
    ----------
    m : array-like (4,)
        2x2 matrix stored in row-major order.

    Returns
    -------
    invm : array (4,)
        2x2 inverse matrix.
    """
    a, b, c, d = m[0], m[1], m[2], m[3]
    detm = a * d - b * c

    # if detm == 0:
    #     raise ValueError("Matrix is singular and cannot be inverted.")

    invm = np.zeros(4, dtype=np.float64)
    invm[0] = d / detm
    invm[1] = -b / detm
    invm[2] = -c / detm
    invm[3] = a / detm

    return invm


@njit
def inv3by3(m):
    """
    Invert a 3x3 matrix.

    Parameters
    ----------
    m : array-like (9,)
        3x3 matrix stored in row-major order.

    Returns
    -------
    invm : array (9,)
        3x3 inverse matrix.
    """
    a, b, c = m[0], m[1], m[2]
    d, e, f = m[3], m[4], m[5]
    g, h, i = m[6], m[7], m[8]

    aa = e * i - f * h
    bb = -(d * i - f * g)
    cc = d * h - e * g
    dd = -(b * i - c * h)
    ee = a * i - c * g
    ff = -(a * h - b * g)
    gg = b * f - c * e
    hh = -(a * f - c * d)
    ii = a * e - b * d

    detm = a * aa + b * bb + c * cc

    # if detm == 0:
    #     raise ValueError("Matrix is singular and cannot be inverted.")

    invm = np.zeros(9, dtype=np.float64)
    invm[0] = aa / detm
    invm[1] = dd / detm
    invm[2] = gg / detm
    invm[3] = bb / detm
    invm[4] = ee / detm
    invm[5] = hh / detm
    invm[6] = cc / detm
    invm[7] = ff / detm
    invm[8] = ii / detm

    return invm


@njit
def eig2by2(m):
    """
    Compute the eigenvalues of a 2x2 matrix.

    Parameters
    ----------
    m : array-like (4,)
        2x2 matrix stored in row-major order.

    Returns
    -------
    eig : array (2,)
        Eigenvalues sorted in ascending order.
    """
    m00, m01, m10, m11 = m[0], m[1], m[2], m[3]

    term = np.sqrt(m00**2 + m11**2 - 2 * m00 * m11 + 4 * m01 * m10)
    eig1 = 0.5 * (m00 + m11 + term)
    eig2 = 0.5 * (m00 + m11 - term)

    return np.sort(np.array([eig1, eig2]))


@njit
def symeig3by3(m):
    """
    Compute the eigenvalues of a 3x3 symmetric matrix.

    Parameters
    ----------
    m : array-like (9,)
        3x3 symmetric matrix stored in row-major order.

    Returns
    -------
    eig : array (3,)
        Eigenvalues sorted in ascending order.
    """
    pi = 4 * np.arctan(1.0)

    m00, m01, m02 = m[0], m[1], m[2]
    m11, m12, m22 = m[4], m[5], m[8]

    alpha = m00 + m11 + m22
    beta = m01**2 + m02**2 + m12**2 - m00 * m11 - m11 * m22 - m22 * m00
    gamma = m00 * m11 * m22 + 2 * m01 * m12 * m02 - m00 * m12**2 - m22 * m01**2 - m11 * m02**2

    p = - (3 * beta + alpha**2) / 3
    q = - (gamma + (2.0 / 27.0) * alpha**3 + alpha * beta / 3.0)
    phi = np.arccos(-q / (2 * (abs(p) / 3) ** 1.5))

    eig1 = alpha / 3 + 2 * np.sqrt(abs(p) / 3) * np.cos(phi / 3)
    eig2 = alpha / 3 - 2 * np.sqrt(abs(p) / 3) * np.cos((phi - pi) / 3)
    eig3 = alpha / 3 - 2 * np.sqrt(abs(p) / 3) * np.cos((phi + pi) / 3)

    return np.sort(np.array([eig1, eig2, eig3]))  # Ensure sorted order
