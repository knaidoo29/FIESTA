import numpy as np
import pytest
import shift
from fiesta.p2g import part2grid2D, part2grid3D, get_deconvol_p
from fiesta.p2g.deconvol import (
    get_part2grid2D_kernel,
    get_part2grid3D_kernel,
    deconvolve_part2grid_2D,
    deconvolve_part2grid_3D,
)


class TestP2G:
    def test_part2grid2D_ngp(self):
        x = np.random.rand(10)
        y = np.random.rand(10)
        f = np.random.rand(10)
        grid = part2grid2D(x, y, f, boxsize=1.0, ngrid=10, method="NGP")
        assert grid.shape == (10, 10)

    def test_part2grid2D_ngp_list(self):
        x = np.random.rand(10)
        y = np.random.rand(10)
        f = np.random.rand(10)
        grid = part2grid2D(
            x,
            y,
            f,
            boxsize=[1.0, 1.0],
            ngrid=[10, 10],
            origin=[0.0, 0.0],
            periodic=[True, False],
            method="NGP",
        )
        assert grid.shape == (10, 10)

    def test_part2grid2D_cic(self):
        x = np.random.rand(10)
        y = np.random.rand(10)
        f = np.random.rand(10)
        grid = part2grid2D(x, y, f, boxsize=1.0, ngrid=10, method="CIC")
        assert grid.shape == (10, 10)

    def test_part2grid2D_tsc(self):
        x = np.random.rand(10)
        y = np.random.rand(10)
        f = np.random.rand(10)
        grid = part2grid2D(x, y, f, boxsize=1.0, ngrid=10, method="TSC")
        assert grid.shape == (10, 10)

    def test_part2grid2D_pcs(self):
        x = np.random.rand(10)
        y = np.random.rand(10)
        f = np.random.rand(10)
        grid = part2grid2D(x, y, f, boxsize=1.0, ngrid=10, method="PCS")
        assert grid.shape == (10, 10)

    def test_part2grid3D_ngp(self):
        x = np.random.rand(10)
        y = np.random.rand(10)
        z = np.random.rand(10)
        f = np.random.rand(10)
        grid = part2grid3D(x, y, z, f, boxsize=1.0, ngrid=5, method="NGP")
        assert grid.shape == (5, 5, 5)

    def test_part2grid3D_ngp_list(self):
        x = np.random.rand(10)
        y = np.random.rand(10)
        z = np.random.rand(10)
        f = np.random.rand(10)
        grid = part2grid3D(
            x,
            y,
            z,
            f,
            boxsize=[1.0, 1.0, 1.0],
            ngrid=[5, 5, 5],
            origin=[0.0, 0.0, 0.0],
            periodic=[True, False, True],
            method="NGP",
        )
        assert grid.shape == (5, 5, 5)

    def test_part2grid3D_cic(self):
        x = np.random.rand(10)
        y = np.random.rand(10)
        z = np.random.rand(10)
        f = np.random.rand(10)
        grid = part2grid3D(x, y, z, f, boxsize=1.0, ngrid=5, method="CIC")
        assert grid.shape == (5, 5, 5)

    def test_part2grid3D_tsc(self):
        x = np.random.rand(10)
        y = np.random.rand(10)
        z = np.random.rand(10)
        f = np.random.rand(10)
        grid = part2grid3D(x, y, z, f, boxsize=1.0, ngrid=5, method="TSC")
        assert grid.shape == (5, 5, 5)

    def test_part2grid3D_pcs(self):
        x = np.random.rand(10)
        y = np.random.rand(10)
        z = np.random.rand(10)
        f = np.random.rand(10)
        grid = part2grid3D(x, y, z, f, boxsize=1.0, ngrid=5, method="PCS")
        assert grid.shape == (5, 5, 5)

    def test_get_deconvol_p(self):
        assert get_deconvol_p("NGP") == 1
        assert get_deconvol_p("CIC") == 2
        assert get_deconvol_p("TSC") == 3
        assert get_deconvol_p("PCS") == 4

    def test_get_part2grid2D_kernel_scalar(self):
        """Test get_part2grid2D_kernel with scalar ngrid and boxsize."""
        kx2d, ky2d = shift.cart.kgrid2D(1.0, 10)
        Wk = get_part2grid2D_kernel(kx2d, ky2d, ngrid=10, boxsize=1.0, method="TSC")
        assert Wk.shape == kx2d.shape
        assert np.all(np.isfinite(Wk))

    def test_get_part2grid2D_kernel_list(self):
        """Test get_part2grid2D_kernel with list ngrid and boxsize."""
        kx2d, ky2d = shift.cart.kgrid2D([1.0, 0.8], (10, 8))
        Wk = get_part2grid2D_kernel(
            kx2d, ky2d, ngrid=[10, 8], boxsize=[1.0, 0.8], method="TSC"
        )
        assert Wk.shape == kx2d.shape
        assert np.all(np.isfinite(Wk))

    def test_get_part2grid3D_kernel_scalar(self):
        """Test get_part2grid3D_kernel with scalar ngrid and boxsize."""
        kx3d, ky3d, kz3d = shift.cart.kgrid3D(1.0, 8)
        Wk = get_part2grid3D_kernel(
            kx3d, ky3d, kz3d, ngrid=8, boxsize=1.0, method="TSC"
        )
        assert Wk.shape == kx3d.shape
        assert np.all(np.isfinite(Wk))

    def test_get_part2grid3D_kernel_list(self):
        """Test get_part2grid3D_kernel with list ngrid and boxsize."""
        kx3d, ky3d, kz3d = shift.cart.kgrid3D([1.0, 0.8, 0.6], (8, 6, 4))
        Wk = get_part2grid3D_kernel(
            kx3d, ky3d, kz3d, ngrid=[8, 6, 4], boxsize=[1.0, 0.8, 0.6], method="TSC"
        )
        assert Wk.shape == kx3d.shape
        assert np.all(np.isfinite(Wk))

    def test_deconvolve_part2grid_2D(self):
        """Test deconvolve_part2grid_2D function."""
        # Create a simple 2D field
        field = np.random.rand(16, 16)
        deconvolved = deconvolve_part2grid_2D(field, boxsize=1.0, method="TSC")
        assert deconvolved.shape == field.shape
        assert np.all(np.isfinite(deconvolved))

    def test_deconvolve_part2grid_2D_list_params(self):
        """Test deconvolve_part2grid_2D with list boxsize."""
        field = np.random.rand(16, 12)
        deconvolved = deconvolve_part2grid_2D(field, boxsize=[1.0, 0.75], method="TSC")
        assert deconvolved.shape == field.shape
        assert np.all(np.isfinite(deconvolved))

    def test_deconvolve_part2grid_3D(self):
        """Test deconvolve_part2grid_3D function."""
        # Create a simple 3D field
        field = np.random.rand(8, 8, 8)
        deconvolved = deconvolve_part2grid_3D(field, boxsize=1.0, method="TSC")
        assert deconvolved.shape == field.shape
        assert np.all(np.isfinite(deconvolved))

    def test_deconvolve_part2grid_3D_list_params(self):
        """Test deconvolve_part2grid_3D with list boxsize."""
        field = np.random.rand(8, 6, 4)
        deconvolved = deconvolve_part2grid_3D(
            field, boxsize=[1.0, 0.75, 0.5], method="TSC"
        )
        assert deconvolved.shape == field.shape
        assert np.all(np.isfinite(deconvolved))
