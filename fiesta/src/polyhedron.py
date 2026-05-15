import numpy as np
from numba import njit


@njit
def tetrahedron_volume(
    xa: float,
    ya: float,
    za: float,
    xb: float,
    yb: float,
    zb: float,
    xc: float,
    yc: float,
    zc: float,
    xd: float,
    yd: float,
    zd: float,
) -> float:  # pragma: no cover
    """
    Computes the volume of a tetrahedron from its vertices using the determinant method.

    Parameters
    ----------
    xa, ya, za : float
        Coordinates of point A.
    xb, yb, zb : float
        Coordinates of point B.
    xc, yc, zc : float
        Coordinates of point C.
    xd, yd, zd : float
        Coordinates of point D.

    Returns
    -------
    float
        Volume of the tetrahedron.
    """

    a, b, c = xa - xd, ya - yd, za - zd
    d, e, f = xb - xd, yb - yd, zb - zd
    g, h, i = xc - xd, yc - yd, zc - zd

    det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)

    return abs(det) / 6.0
