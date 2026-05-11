import numpy as np
import pytest
from fiesta.utils import (
    get_vector_magnitude_2D,
    get_vector_magnitude_3D,
    complex_mult,
    complex_div,
    flatten_list,
)


class TestUtils:
    def test_get_vector_magnitude_2D(self):
        x = np.array([1, 2, 3])
        y = np.array([4, 5, 6])
        mag = get_vector_magnitude_2D(x, y)
        expected = np.sqrt(x**2 + y**2)
        assert np.allclose(mag, expected)

    def test_get_vector_magnitude_3D(self):
        x = np.array([1, 2, 3])
        y = np.array([4, 5, 6])
        z = np.array([7, 8, 9])
        mag = get_vector_magnitude_3D(x, y, z)
        expected = np.sqrt(x**2 + y**2 + z**2)
        assert np.allclose(mag, expected)

    def test_complex_mult(self):
        arr = np.array([1+2j, 3+4j])
        factors = 2.0
        result = complex_mult(arr, factors)
        expected = arr * factors
        assert np.allclose(result, expected)

    def test_complex_div(self):
        arr = np.array([2+4j, 6+8j])
        factors = 2.0
        result = complex_div(arr, factors)
        expected = arr / factors
        assert np.allclose(result, expected)

    def test_flatten_list(self):
        nested = [1, [2, 3], [4, [5, 6]]]
        flat = flatten_list(nested)
        assert flat == [1, 2, 3, 4, 5, 6]