import numpy as np
import pytest
from fiesta.interp import bilinear, trilinear


class TestInterp:
    def test_bilinear(self):
        fgrid = np.random.rand(10, 10)
        boxsize = 1.0
        x = np.array([0.5, 0.3])
        y = np.array([0.5, 0.7])
        result = bilinear(fgrid, boxsize, x, y, periodic=True)
        result = bilinear(fgrid, boxsize, x, y, periodic=False)
        assert len(result) == 2
        assert not np.isnan(result).any()

    def test_bilinear_non_periodic(self):
        fgrid = np.random.rand(10, 10)
        boxsize = [1.0, 1.0]
        x = np.array([0.5, 0.3])
        y = np.array([0.5, 0.7])
        result = bilinear(
            fgrid, boxsize, x, y, periodic=[False, False], origin=[0.0, 0.0]
        )
        result = bilinear(
            fgrid, boxsize, x, y, periodic=[True, True], origin=[0.0, 0.0]
        )
        assert len(result) == 2

        x = np.array([0.5, 0.3, 10.0])
        y = np.array([0.5, 0.7, 10.0])
        result = bilinear(fgrid, boxsize, x, y, periodic=True, origin=[0.0, 0.0])
        result = bilinear(fgrid, boxsize, x, y, periodic=False, origin=[0.0, 0.0])
        result = bilinear(
            fgrid, boxsize, x, y, periodic=[False, False], origin=[0.0, 0.0]
        )
        result = bilinear(
            fgrid, boxsize, x, y, periodic=[True, True], origin=[0.0, 0.0]
        )
        assert len(result) == 3

    def test_trilinear(self):
        fgrid = np.random.rand(5, 5, 5)
        boxsize = 1.0
        x = np.array([0.5, 0.3])
        y = np.array([0.5, 0.7])
        z = np.array([0.2, 0.8])
        result = trilinear(fgrid, boxsize, x, y, z, periodic=True)
        result = trilinear(fgrid, boxsize, x, y, z, periodic=False)
        assert len(result) == 2
        assert not np.isnan(result).any()

    def test_trilinear_non_periodic(self):
        fgrid = np.random.rand(5, 5, 5)
        boxsize = [1.0, 1.0, 1.0]
        x = np.array([0.5, 0.3])
        y = np.array([0.5, 0.7])
        z = np.array([0.2, 0.8])
        result = trilinear(
            fgrid,
            boxsize,
            x,
            y,
            z,
            periodic=[False, True, False],
            origin=[0.0, 0.0, 0.0],
        )
        result = trilinear(
            fgrid,
            boxsize,
            x,
            y,
            z,
            periodic=[True, False, True],
            origin=[0.0, 0.0, 0.0],
        )
        assert len(result) == 2

        x = np.array([0.5, 0.3, 1.1])
        y = np.array([0.5, 0.7, 1.1])
        z = np.array([0.2, 0.8, 1.1])
        result = trilinear(
            fgrid, boxsize, x, y, z, periodic=True, origin=[0.0, 0.0, 0.0]
        )
        result = trilinear(
            fgrid, boxsize, x, y, z, periodic=False, origin=[0.0, 0.0, 0.0]
        )
        result = trilinear(
            fgrid,
            boxsize,
            x,
            y,
            z,
            periodic=[False, True, False],
            origin=[0.0, 0.0, 0.0],
        )
        result = trilinear(
            fgrid,
            boxsize,
            x,
            y,
            z,
            periodic=[True, False, True],
            origin=[0.0, 0.0, 0.0],
        )
        assert len(result) == 3
