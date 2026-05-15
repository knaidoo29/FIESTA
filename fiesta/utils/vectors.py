import numpy as np


def get_vector_magnitude_2D(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Returns the magnitude for vectors.

    Parameters
    ----------
    x : array
        X component.
    y : array
        Y component.

    Returns
    -------
    mag : array
        Magnitude of the vector.
    """
    mag = np.sqrt(x**2.0 + y**2.0)
    return mag


def get_vector_magnitude_3D(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    """Returns the magnitude for vectors.

    Parameters
    ----------
    x : array
        X component.
    y : array
        Y component.
    z : array
        Z component.

    Returns
    -------
    mag : array
        Magnitude of the vector.
    """
    mag = np.sqrt(x**2.0 + y**2.0 + z**2.0)
    return mag
