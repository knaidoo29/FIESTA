import numpy as np

from typing import Union, Optional


def complex_mult(
    complex_array: np.ndarray, factors: Union[float, np.ndarray]
) -> np.ndarray:
    """Complex array multiplication.

    Parameters
    ----------
    complex_array : array
        Complex array.
    factors : float or array
        Factors to be multiplied by.
    """
    complex_array.real *= factors
    complex_array.imag *= factors
    return complex_array


def complex_div(
    complex_array: np.ndarray, factors: Union[float, np.ndarray]
) -> np.ndarray:
    """Complex array division.

    Parameters
    ----------
    complex_array : array
        Complex array.
    factors : float/array
        Factors to be dividied by.
    """
    complex_array.real /= factors
    complex_array.imag /= factors
    return complex_array
