import numpy as np
import pytest
from fiesta.p2g import part2grid2D, part2grid3D, get_deconvol_p


class TestP2G:
    def test_part2grid2D_ngp(self):
        x = np.random.rand(10)
        y = np.random.rand(10)
        f = np.random.rand(10)
        grid = part2grid2D(x, y, f, boxsize=1.0, ngrid=10, method='NGP')
        assert grid.shape == (10, 10)

    def test_part2grid2D_cic(self):
        x = np.random.rand(10)
        y = np.random.rand(10)
        f = np.random.rand(10)
        grid = part2grid2D(x, y, f, boxsize=1.0, ngrid=10, method='CIC')
        assert grid.shape == (10, 10)

    def test_part2grid2D_tsc(self):
        x = np.random.rand(10)
        y = np.random.rand(10)
        f = np.random.rand(10)
        grid = part2grid2D(x, y, f, boxsize=1.0, ngrid=10, method='TSC')
        assert grid.shape == (10, 10)

    def test_part2grid2D_pcs(self):
        x = np.random.rand(10)
        y = np.random.rand(10)
        f = np.random.rand(10)
        grid = part2grid2D(x, y, f, boxsize=1.0, ngrid=10, method='PCS')
        assert grid.shape == (10, 10)

    def test_part2grid3D_ngp(self):
        x = np.random.rand(10)
        y = np.random.rand(10)
        z = np.random.rand(10)
        f = np.random.rand(10)
        grid = part2grid3D(x, y, z, f, boxsize=1.0, ngrid=5, method='NGP')
        assert grid.shape == (5, 5, 5)

    def test_part2grid3D_cic(self):
        x = np.random.rand(10)
        y = np.random.rand(10)
        z = np.random.rand(10)
        f = np.random.rand(10)
        grid = part2grid3D(x, y, z, f, boxsize=1.0, ngrid=5, method='CIC')
        assert grid.shape == (5, 5, 5)

    def test_part2grid3D_tsc(self):
        x = np.random.rand(10)
        y = np.random.rand(10)
        z = np.random.rand(10)
        f = np.random.rand(10)
        grid = part2grid3D(x, y, z, f, boxsize=1.0, ngrid=5, method='TSC')
        assert grid.shape == (5, 5, 5)

    def test_part2grid3D_pcs(self):
        x = np.random.rand(10)
        y = np.random.rand(10)
        z = np.random.rand(10)
        f = np.random.rand(10)
        grid = part2grid3D(x, y, z, f, boxsize=1.0, ngrid=5, method='PCS')
        assert grid.shape == (5, 5, 5)

    def test_get_deconvol_p(self):
        assert get_deconvol_p('NGP') == 1
        assert get_deconvol_p('CIC') == 2
        assert get_deconvol_p('TSC') == 3
        assert get_deconvol_p('PCS') == 4