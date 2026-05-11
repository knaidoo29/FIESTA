import numpy as np
import pytest
from fiesta.dtfe import Delaunay2D, Delaunay3D


class TestDTFE:
    def test_delaunay2d_setup(self):
        x = np.random.rand(10)
        y = np.random.rand(10)
        delaunay = Delaunay2D()
        delaunay.setup(x, y, boundary=[0, 1, 0, 1])
        assert delaunay.x is not None
        assert delaunay.y is not None

    def test_delaunay2d_construct(self):
        x = np.random.rand(10)
        y = np.random.rand(10)
        delaunay = Delaunay2D()
        delaunay.setup(x, y, boundary=[0, 1, 0, 1])
        delaunay.construct()
        # After construct, should have triangulation
        assert hasattr(delaunay, 'tri')

    def test_delaunay2d_get_area(self):
        x = np.random.rand(10)
        y = np.random.rand(10)
        delaunay = Delaunay2D()
        delaunay.setup(x, y, boundary=[0, 1, 0, 1])
        delaunay.construct()
        areas = delaunay.get_area()
        assert len(areas) == len(x)

    def test_delaunay2d_get_dens(self):
        x = np.random.rand(10)
        y = np.random.rand(10)
        delaunay = Delaunay2D()
        delaunay.setup(x, y, boundary=[0, 1, 0, 1])
        delaunay.construct()
        dens = delaunay.get_dens()
        assert len(dens) == len(x)

    def test_delaunay3d_setup(self):
        x = np.random.rand(10)
        y = np.random.rand(10)
        z = np.random.rand(10)
        delaunay = Delaunay3D()
        delaunay.setup(x, y, z, boundary=[0, 1, 0, 1, 0, 1])
        assert delaunay.x is not None
        assert delaunay.y is not None
        assert delaunay.z is not None

    def test_delaunay3d_construct(self):
        x = np.random.rand(10)
        y = np.random.rand(10)
        z = np.random.rand(10)
        delaunay = Delaunay3D()
        delaunay.setup(x, y, z, boundary=[0, 1, 0, 1, 0, 1])
        delaunay.construct()
        # After construct, should have triangulation
        assert hasattr(delaunay, 'tri')