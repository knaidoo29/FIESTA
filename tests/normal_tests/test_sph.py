import numpy as np
import pytest
from fiesta.sph import cubic_kernel, dcubic_kernel, SPH1D, SPH2D, SPH3D


class TestSPH:
    def test_cubic_kernel_1d(self):
        r = 0.5
        h = 1.0
        result = cubic_kernel(r, h, dim=1)
        r = 1.5
        h = 1.0
        result = cubic_kernel(r, h, dim=1)
        r = 2.5
        h = 1.0
        result = cubic_kernel(r, h, dim=1)
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

    def test_dcubic_kernel_1d(self):
        r = np.array([0.0, 0.5, 1.0])
        h = 1.0
        result = dcubic_kernel(r, h, dim=1)
        assert len(result) == 3

    def test_dcubic_kernel_2d(self):
        r = np.array([0.0, 0.5, 1.0])
        h = 1.0
        result = dcubic_kernel(r, h, dim=2)
        assert len(result) == 3

    def test_dcubic_kernel_3d(self):
        r = np.array([0.0, 0.5, 1.0])
        h = 1.0
        result = dcubic_kernel(r, h, dim=3)
        assert len(result) == 3

    def test_dcubic_kernel(self):
        r = 0.5
        h = 1.0
        result = dcubic_kernel(r, h, dim=3)
        r = 1.5
        h = 1.0
        result = dcubic_kernel(r, h, dim=3)
        r = 2.5
        h = 1.0
        result = dcubic_kernel(r, h, dim=3)
        r = np.array([0.0, 0.5, 1.0])
        h = 1.0
        result = dcubic_kernel(r, h, dim=3)
        assert len(result) == 3

    def test_sph1d(self):
        x = np.random.rand(100)
        sph = SPH1D()
        sph.build_tree(x)
        sph.setup(k=10)
        density = sph.get_density(x)
        assert len(density) == len(x)

        x = np.random.rand(100)
        sph = SPH1D()
        sph.build_tree(x)
        sph.setup(k=10, mass=np.random.rand(100))
        density = sph.get_density(x)
        assert len(density) == len(x)
        field = sph.sph_estimate(x, f=None)
        field = sph.sph_estimate(x, f=np.random.rand(100))
        sph.set_field(np.random.rand(100))
        sph.estimate(x, dens=np.random.rand(100))

    def test_sph2d(self):
        x = np.random.rand(100)
        y = np.random.rand(100)
        sph = SPH2D()
        sph.build_tree(x, y)
        sph.setup(k=10)
        density = sph.get_density(x, y)
        assert len(density) == len(x)

        x = np.random.rand(100)
        y = np.random.rand(100)
        sph = SPH2D()
        sph.build_tree(x, y)
        sph.setup(k=10, mass=np.random.rand(100))
        density = sph.get_density(x, y)
        assert len(density) == len(x)
        field = sph.sph_estimate(x, y, f=None)
        field = sph.sph_estimate(x, y, f=np.random.rand(100))
        sph.set_field(np.random.rand(100))
        sph.estimate(x, y, dens=np.random.rand(100))

    def test_sph3d(self):
        x = np.random.rand(100)
        y = np.random.rand(100)
        z = np.random.rand(100)
        sph = SPH3D()
        sph.build_tree(x, y, z)
        sph.setup(k=10)
        density = sph.get_density(x, y, z)
        assert len(density) == len(x)

        x = np.random.rand(100)
        y = np.random.rand(100)
        z = np.random.rand(100)
        sph = SPH3D()
        sph.build_tree(x, y, z)
        sph.setup(k=10, mass=np.random.rand(100))
        density = sph.get_density(x, y, z)
        assert len(density) == len(x)
        field = sph.sph_estimate(x, y, z, f=None)
        field = sph.sph_estimate(x, y, z, f=np.random.rand(100))
        sph.set_field(np.random.rand(100))
        sph.estimate(x, y, z, dens=np.random.rand(100))
