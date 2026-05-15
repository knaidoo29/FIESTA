import numpy as np
import pytest
from fiesta.maths import dfdx, dfdy, dfdz, det3by3, det4by4


class TestMaths:
    def test_dfdx(self):
        xgrid = np.linspace(0, 1, 10)
        f = np.sin(xgrid)
        result = dfdx(xgrid, f)
        assert len(result) == len(f)

    def test_dfdx_periodic(self):
        xgrid = np.linspace(0, 1, 10)
        f = np.sin(xgrid)
        result = dfdx(xgrid, f, periodic=True)
        assert len(result) == len(f)

    def test_dfdy(self):
        ygrid = np.linspace(0, 1, 10)
        xgrid = np.linspace(0, 1, 5)
        X, Y = np.meshgrid(xgrid, ygrid)
        f = np.sin(Y)
        result = dfdy(ygrid, f)
        assert result.shape == f.shape

    def test_dfdz(self):
        zgrid = np.linspace(0, 1, 10)
        ygrid = np.linspace(0, 1, 8)
        xgrid = np.linspace(0, 1, 5)
        X, Y, Z = np.meshgrid(xgrid, ygrid, zgrid, indexing="ij")
        f = np.sin(Z)
        result = dfdz(zgrid, f)
        assert result.shape == f.shape

    def test_det3by3(self):
        det = det3by3(1, 0, 0, 0, 1, 0, 0, 0, 1)
        assert det == 1.0

    def test_det3by3_array(self):
        m00 = np.array([1, 2])
        det = det3by3(m00, 0, 0, 0, 1, 0, 0, 0, 1)
        assert len(det) == 2

    def test_det4by4(self):
        det = det4by4(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1)
        assert det == 1.0

    def test_det4by4_array(self):
        m00 = np.array([1, 2])
        det = det4by4(m00, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1)
        assert len(det) == 2
