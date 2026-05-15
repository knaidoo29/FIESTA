import numpy as np
import pytest
from fiesta.randoms import random_uniform, random_box, random_cube


class TestRandoms:
    def test_random_uniform(self):
        result = random_uniform(10, 0.0, 1.0)
        assert len(result) == 10
        assert np.all((result >= 0.0) & (result <= 1.0))

    def test_random_box(self):
        x, y = random_box(10, 0.0, 1.0, 0.0, 1.0)
        assert len(x) == 10
        assert len(y) == 10
        assert np.all((x >= 0.0) & (x <= 1.0))
        assert np.all((y >= 0.0) & (y <= 1.0))

    def test_random_cube(self):
        x, y, z = random_cube(10, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
        assert len(x) == 10
        assert len(y) == 10
        assert len(z) == 10
        assert np.all((x >= 0.0) & (x <= 1.0))
        assert np.all((y >= 0.0) & (y <= 1.0))
        assert np.all((z >= 0.0) & (z <= 1.0))
