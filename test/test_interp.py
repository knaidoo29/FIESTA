import numpy as np
import pytest
from fiesta.interp import bilinear, trilinear


class TestInterp:
    def test_bilinear(self):
        fgrid = np.random.rand(10, 10)
        boxsize = 1.0
        x = np.array([0.5, 0.3])
        y = np.array([0.5, 0.7])
        result = bilinear(fgrid, boxsize, x, y)
        assert len(result) == 2
        assert not np.isnan(result).any()

    def test_bilinear_non_periodic(self):
        fgrid = np.random.rand(10, 10)
        boxsize = 1.0
        x = np.array([0.5, 0.3])
        y = np.array([0.5, 0.7])
        result = bilinear(fgrid, boxsize, x, y, periodic=False)
        assert len(result) == 2

    def test_trilinear(self):
        fgrid = np.random.rand(5, 5, 5)
        boxsize = 1.0
        x = np.array([0.5, 0.3])
        y = np.array([0.5, 0.7])
        z = np.array([0.2, 0.8])
        result = trilinear(fgrid, boxsize, x, y, z)
        assert len(result) == 2
        assert not np.isnan(result).any()

    def test_trilinear_non_periodic(self):
        fgrid = np.random.rand(5, 5, 5)
        boxsize = 1.0
        x = np.array([0.5, 0.3])
        y = np.array([0.5, 0.7])
        z = np.array([0.2, 0.8])
        result = trilinear(fgrid, boxsize, x, y, z, periodic=False)
        assert len(result) == 2