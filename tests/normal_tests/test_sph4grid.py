import numpy as np
import pytest
from fiesta.sph import sph4grid2D, sph4grid3D

RNG = np.random.default_rng(42)
N = 100
BOXSIZE = 1.0

# Shared random data
X2, Y2 = RNG.uniform(0, BOXSIZE, N), RNG.uniform(0, BOXSIZE, N)
X3, Y3, Z3 = RNG.uniform(0, BOXSIZE, N), RNG.uniform(0, BOXSIZE, N), RNG.uniform(0, BOXSIZE, N)
F2 = RNG.standard_normal(N)
F3 = RNG.standard_normal(N)


# ===========================================================================
# sph4grid2D
# ===========================================================================

class TestSph4grid2DShape:

    def test_scalar_ngrid(self):
        result = sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=8)
        assert result.shape == (8, 8)

    def test_list_ngrid(self):
        result = sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=[4, 6])
        assert result.shape == (4, 6)

    def test_field_shape(self):
        result = sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=8, f=F2)
        assert result.shape == (8, 8)

    def test_outputgrid_density_returns_tuple(self):
        dgrid, x2d, y2d = sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=8, outputgrid=True)
        assert dgrid.shape == x2d.shape == y2d.shape == (8, 8)

    def test_outputgrid_field_returns_tuple(self):
        fgrid, x2d, y2d = sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=8, f=F2, outputgrid=True)
        assert fgrid.shape == x2d.shape == y2d.shape == (8, 8)

    def test_list_boxsize(self):
        result = sph4grid2D(X2, Y2, boxsize=[1.0, 2.0], ngrid=8)
        assert result.shape == (8, 8)

    def test_list_ngrid_list_boxsize(self):
        result = sph4grid2D(X2, Y2, boxsize=[1.0, 2.0], ngrid=[4, 6])
        assert result.shape == (4, 6)


class TestSph4grid2DValues:

    def test_density_is_finite(self):
        result = sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=8)
        assert np.all(np.isfinite(result))

    def test_density_is_nonnegative(self):
        result = sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=8)
        assert np.all(result >= 0)

    def test_field_is_finite(self):
        result = sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=8, f=F2)
        assert np.all(np.isfinite(result))

    def test_subsampling_does_not_change_shape(self):
        r1 = sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=8, subsampling=1)
        r2 = sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=8, subsampling=2)
        assert r1.shape == r2.shape

    def test_explicit_mass(self):
        mass = RNG.uniform(0.5, 1.5, N)
        result = sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=8, mass=mass)
        assert result.shape == (8, 8)
        assert np.all(np.isfinite(result))

    def test_default_mass_matches_unit_mass(self):
        r_default = sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=8)
        r_ones = sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=8, mass=np.ones(N))
        np.testing.assert_array_equal(r_default, r_ones)


class TestSph4grid2DOptions:

    def test_periodic_false(self):
        result = sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=8, periodic=False)
        assert result.shape == (8, 8)

    def test_scalar_origin(self):
        result = sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=8, origin=0.0)
        assert result.shape == (8, 8)

    def test_list_origin(self):
        result = sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=8, origin=[0.0, 0.0])
        assert result.shape == (8, 8)

    def test_custom_k(self):
        result = sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=8, k=10)
        assert result.shape == (8, 8)


class TestSph4grid2DListParams:

    # --- boxsize ---
    def test_list_boxsize_symmetric(self):
        result = sph4grid2D(X2, Y2, boxsize=[1.0, 1.0], ngrid=8)
        assert result.shape == (8, 8)

    def test_list_boxsize_asymmetric(self):
        result = sph4grid2D(X2, Y2, boxsize=[1.0, 2.0], ngrid=8)
        assert result.shape == (8, 8)

    def test_list_boxsize_matches_scalar_boxsize(self):
        r_scalar = sph4grid2D(X2, Y2, boxsize=1.0, ngrid=8)
        r_list   = sph4grid2D(X2, Y2, boxsize=[1.0, 1.0], ngrid=8)
        np.testing.assert_array_equal(r_scalar, r_list)

    # --- ngrid ---
    def test_list_ngrid_symmetric(self):
        result = sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=[8, 8])
        assert result.shape == (8, 8)

    def test_list_ngrid_asymmetric(self):
        result = sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=[4, 6])
        assert result.shape == (4, 6)

    def test_list_ngrid_matches_scalar_ngrid(self):
        r_scalar = sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=8)
        r_list   = sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=[8, 8])
        np.testing.assert_array_equal(r_scalar, r_list)

    # --- origin ---
    def test_list_origin_zeros(self):
        result = sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=8, origin=[0.0, 0.0])
        assert result.shape == (8, 8)

    def test_list_origin_asymmetric(self):
        result = sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=8, origin=[0.1, 0.2])
        assert result.shape == (8, 8)

    def test_list_origin_matches_scalar_origin(self):
        r_scalar = sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=8, origin=0.0)
        r_list   = sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=8, origin=[0.0, 0.0])
        np.testing.assert_array_equal(r_scalar, r_list)

    # --- all list ---
    def test_all_list_params(self):
        result = sph4grid2D(X2, Y2, boxsize=[1.0, 2.0], ngrid=[4, 6], origin=[0.1, 0.2], subsampling=2)
        assert result.shape == (4, 6)

    def test_all_list_params_finite(self):
        result = sph4grid2D(X2, Y2, boxsize=[1.0, 2.0], ngrid=[4, 6], origin=[0.1, 0.2], subsampling=2)
        assert np.all(np.isfinite(result))


# ===========================================================================
# sph4grid3D
# ===========================================================================

class TestSph4grid3DShape:

    def test_scalar_ngrid(self):
        result = sph4grid3D(X3, Y3, Z3, boxsize=BOXSIZE, ngrid=4)
        assert result.shape == (4, 4, 4)

    def test_list_ngrid(self):
        result = sph4grid3D(X3, Y3, Z3, boxsize=BOXSIZE, ngrid=[3, 4, 5])
        assert result.shape == (3, 4, 5)

    def test_field_shape(self):
        result = sph4grid3D(X3, Y3, Z3, boxsize=BOXSIZE, ngrid=4, f=F3)
        assert result.shape == (4, 4, 4)

    def test_outputgrid_density_returns_tuple(self):
        dgrid, x3d, y3d, z3d = sph4grid3D(X3, Y3, Z3, boxsize=BOXSIZE, ngrid=4, outputgrid=True)
        assert dgrid.shape == x3d.shape == y3d.shape == z3d.shape == (4, 4, 4)

    def test_outputgrid_field_returns_tuple(self):
        fgrid, x3d, y3d, z3d = sph4grid3D(X3, Y3, Z3, boxsize=BOXSIZE, ngrid=4, f=F3, outputgrid=True)
        assert fgrid.shape == x3d.shape == y3d.shape == z3d.shape == (4, 4, 4)

    def test_list_boxsize(self):
        result = sph4grid3D(X3, Y3, Z3, boxsize=[1.0, 2.0, 3.0], ngrid=4)
        assert result.shape == (4, 4, 4)

    def test_list_ngrid_list_boxsize(self):
        result = sph4grid3D(X3, Y3, Z3, boxsize=[1.0, 2.0, 3.0], ngrid=[3, 4, 5])
        assert result.shape == (3, 4, 5)


class TestSph4grid3DValues:

    def test_density_is_finite(self):
        result = sph4grid3D(X3, Y3, Z3, boxsize=BOXSIZE, ngrid=4)
        assert np.all(np.isfinite(result))

    def test_density_is_nonnegative(self):
        result = sph4grid3D(X3, Y3, Z3, boxsize=BOXSIZE, ngrid=4)
        assert np.all(result >= 0)

    def test_field_is_finite(self):
        result = sph4grid3D(X3, Y3, Z3, boxsize=BOXSIZE, ngrid=4, f=F3)
        assert np.all(np.isfinite(result))

    def test_subsampling_does_not_change_shape(self):
        r1 = sph4grid3D(X3, Y3, Z3, boxsize=BOXSIZE, ngrid=4, subsampling=1)
        r2 = sph4grid3D(X3, Y3, Z3, boxsize=BOXSIZE, ngrid=4, subsampling=2)
        assert r1.shape == r2.shape

    def test_explicit_mass(self):
        mass = RNG.uniform(0.5, 1.5, N)
        result = sph4grid3D(X3, Y3, Z3, boxsize=BOXSIZE, ngrid=4, mass=mass)
        assert result.shape == (4, 4, 4)
        assert np.all(np.isfinite(result))

    def test_default_mass_matches_unit_mass(self):
        r_default = sph4grid3D(X3, Y3, Z3, boxsize=BOXSIZE, ngrid=4)
        r_ones = sph4grid3D(X3, Y3, Z3, boxsize=BOXSIZE, ngrid=4, mass=np.ones(N))
        np.testing.assert_array_equal(r_default, r_ones)


class TestSph4grid3DOptions:

    def test_periodic_false(self):
        result = sph4grid3D(X3, Y3, Z3, boxsize=BOXSIZE, ngrid=4, periodic=False)
        assert result.shape == (4, 4, 4)

    def test_scalar_origin(self):
        result = sph4grid3D(X3, Y3, Z3, boxsize=BOXSIZE, ngrid=4, origin=0.0)
        assert result.shape == (4, 4, 4)

    def test_list_origin(self):
        result = sph4grid3D(X3, Y3, Z3, boxsize=BOXSIZE, ngrid=4, origin=[0.0, 0.0, 0.0])
        assert result.shape == (4, 4, 4)

    def test_custom_k(self):
        result = sph4grid3D(X3, Y3, Z3, boxsize=BOXSIZE, ngrid=4, k=16)
        assert result.shape == (4, 4, 4)


class TestSph4grid3DListParams:

    # --- boxsize ---
    def test_list_boxsize_symmetric(self):
        result = sph4grid3D(X3, Y3, Z3, boxsize=[1.0, 1.0, 1.0], ngrid=4)
        assert result.shape == (4, 4, 4)

    def test_list_boxsize_asymmetric(self):
        result = sph4grid3D(X3, Y3, Z3, boxsize=[1.0, 2.0, 3.0], ngrid=4)
        assert result.shape == (4, 4, 4)

    def test_list_boxsize_matches_scalar_boxsize(self):
        r_scalar = sph4grid3D(X3, Y3, Z3, boxsize=1.0, ngrid=4)
        r_list   = sph4grid3D(X3, Y3, Z3, boxsize=[1.0, 1.0, 1.0], ngrid=4)
        np.testing.assert_array_equal(r_scalar, r_list)

    # --- ngrid ---
    def test_list_ngrid_symmetric(self):
        result = sph4grid3D(X3, Y3, Z3, boxsize=BOXSIZE, ngrid=[4, 4, 4])
        assert result.shape == (4, 4, 4)

    def test_list_ngrid_asymmetric(self):
        result = sph4grid3D(X3, Y3, Z3, boxsize=BOXSIZE, ngrid=[3, 4, 5])
        assert result.shape == (3, 4, 5)

    def test_list_ngrid_matches_scalar_ngrid(self):
        r_scalar = sph4grid3D(X3, Y3, Z3, boxsize=BOXSIZE, ngrid=4)
        r_list   = sph4grid3D(X3, Y3, Z3, boxsize=BOXSIZE, ngrid=[4, 4, 4])
        np.testing.assert_array_equal(r_scalar, r_list)

    # --- origin ---
    def test_list_origin_zeros(self):
        result = sph4grid3D(X3, Y3, Z3, boxsize=BOXSIZE, ngrid=4, origin=[0.0, 0.0, 0.0])
        assert result.shape == (4, 4, 4)

    def test_list_origin_asymmetric(self):
        result = sph4grid3D(X3, Y3, Z3, boxsize=BOXSIZE, ngrid=4, origin=[0.1, 0.2, 0.3])
        assert result.shape == (4, 4, 4)

    def test_list_origin_matches_scalar_origin(self):
        r_scalar = sph4grid3D(X3, Y3, Z3, boxsize=BOXSIZE, ngrid=4, origin=0.0)
        r_list   = sph4grid3D(X3, Y3, Z3, boxsize=BOXSIZE, ngrid=4, origin=[0.0, 0.0, 0.0])
        np.testing.assert_array_equal(r_scalar, r_list)

    # --- all list ---
    def test_all_list_params(self):
        result = sph4grid3D(X3, Y3, Z3, boxsize=[1.0, 2.0, 3.0], ngrid=[3, 4, 5], origin=[0.1, 0.2, 0.3], subsampling=2)
        assert result.shape == (3, 4, 5)

    def test_all_list_params_finite(self):
        result = sph4grid3D(X3, Y3, Z3, boxsize=[1.0, 2.0, 3.0], ngrid=[3, 4, 5], origin=[0.1, 0.2, 0.3], subsampling=2)
        assert np.all(np.isfinite(result))