import numpy as np

from typing import Tuple


def det3by3(
    m00: Tuple[float, np.ndarray],
    m01: Tuple[float, np.ndarray],
    m02: Tuple[float, np.ndarray],
    m10: Tuple[float, np.ndarray],
    m11: Tuple[float, np.ndarray],
    m12: Tuple[float, np.ndarray],
    m20: Tuple[float, np.ndarray],
    m21: Tuple[float, np.ndarray],
    m22: Tuple[float, np.ndarray],
) -> Tuple[float, np.ndarray]:
    """
    Computes the determinant of a matrix or a series of matrices

    Parameters
    ----------
    m00 : float or array
        Matrix element for row 0 and column 0.
    m01 : float or array
        Matrix element for row 0 and column 1.
    m02 : float or array
        Matrix element for row 0 and column 2.
    m10 : float or array
        Matrix element for row 1 and column 0.
    m11 : float or array
        Matrix element for row 1 and column 1.
    m12 : float or array
        Matrix element for row 1 and column 2.
    m20 : float or array
        Matrix element for row 2 and column 0.
    m21 : float or array
        Matrix element for row 2 and column 1.
    m22 : float or array
        Matrix element for row 2 and column 2.

    Returns
    -------
    det : float or array
        Determinant for the matrix/matrices.
    """
    det = (
        m00 * m11 * m22
        + m01 * m12 * m20
        + m02 * m10 * m21
        - m02 * m11 * m20
        - m01 * m10 * m22
        - m00 * m12 * m21
    )
    return det


def det4by4(
    m00: Tuple[float, np.ndarray],
    m01: Tuple[float, np.ndarray],
    m02: Tuple[float, np.ndarray],
    m03: Tuple[float, np.ndarray],
    m10: Tuple[float, np.ndarray],
    m11: Tuple[float, np.ndarray],
    m12: Tuple[float, np.ndarray],
    m13: Tuple[float, np.ndarray],
    m20: Tuple[float, np.ndarray],
    m21: Tuple[float, np.ndarray],
    m22: Tuple[float, np.ndarray],
    m23: Tuple[float, np.ndarray],
    m30: Tuple[float, np.ndarray],
    m31: Tuple[float, np.ndarray],
    m32: Tuple[float, np.ndarray],
    m33: Tuple[float, np.ndarray],
) -> Tuple[float, np.ndarray]:
    """
    Computes the determinant of a matrix or a series of matrices

    Parameters
    ----------
    m00 : float or array
        Matrix element for row 0 and column 0.
    m01 : float or array
        Matrix element for row 0 and column 1.
    m02 : float or array
        Matrix element for row 0 and column 2.
    m03 : float or array
        Matrix element for row 0 and column 3.
    m10 : float or array
        Matrix element for row 1 and column 0.
    m11 : float or array
        Matrix element for row 1 and column 1.
    m12 : float or array
        Matrix element for row 1 and column 2.
    m13 : float or array
        Matrix element for row 1 and column 3.
    m20 : float or array
        Matrix element for row 2 and column 0.
    m21 : float or array
        Matrix element for row 2 and column 1.
    m22 : float or array
        Matrix element for row 2 and column 2.
    m23 : float or array
        Matrix element for row 2 and column 3.
    m30 : float or array
        Matrix element for row 3 and column 0.
    m31 : float or array
        Matrix element for row 3 and column 1.
    m32 : float or array
        Matrix element for row 3 and column 2.
    m33 : float or array
        Matrix element for row 3 and column 3.

    Returns
    -------
    det : float or array
        Determinant for the matrix/matrices.
    """
    det = m00 * det3by3(m11, m12, m13, m21, m22, m23, m31, m32, m33)
    det -= m01 * det3by3(m10, m12, m13, m20, m22, m23, m30, m32, m33)
    det += m02 * det3by3(m10, m11, m13, m20, m21, m23, m30, m31, m33)
    det -= m03 * det3by3(m10, m11, m12, m20, m21, m22, m30, m31, m32)
    return det
