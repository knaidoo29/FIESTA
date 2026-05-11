import numpy as np
import pytest
from fiesta.sph import cubic_kernel, dcubic_kernel, SPH1D, SPH2D, SPH3D


class TestSPH:
    def test_cubic_kernel_1d(self):
        r = np.array([0.0, 0.5, 1.0])
        h = 1.0
        result = cubic_kernel(r, h, dim=1)
        assert len(result) == 3

    def test_cubic_kernel_2d(self):
        r = np.array([0.0, 0.5, 1.0])
        h = 1.0
        result = cubic_kernel(r, h, dim=2)
        assert len(result) == 3

    def test_cubic_kernel_3d(self):
        r = np.array([0.0, 0.5, 1.0])
        h = 1.0
        result = cubic_kernel(r, h, dim=3)
        assert len(result) == 3

    def test_dcubic_kernel(self):
        r = np.array([0.0, 0.5, 1.0])
        h = 1.0
        result = dcubic_kernel(r, h, dim=3)
        assert len(result) == 3

    def test_sph1d(self):
        x = np.random.rand(100)
        sph = SPH1D()
        sph.setup(k=10)
        sph.build_tree(x)
        density = sph.get_density(x)
        assert len(density) == len(x)

    def test_sph2d(self):
        x = np.random.rand(100)
        y = np.random.rand(100)
        sph = SPH2D()
        sph.setup(k=10)
        sph.build_tree(x, y)
        density = sph.get_density(x, y)
        assert len(density) == len(x)

    def test_sph3d(self):
        x = np.random.rand(100)
        y = np.random.rand(100)
        z = np.random.rand(100)
        sph = SPH3D()
        sph.setup(k=10)
        sph.build_tree(x, y, z)
        density = sph.get_density(x, y, z)
        assert len(density) == len(x)