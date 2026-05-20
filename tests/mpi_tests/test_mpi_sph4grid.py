import numpy as np
import pytest

pytest.importorskip("mpi4py")
from shift.mpiutils import MPI

mpi = MPI()

from fiesta.sph import mpi_sph4grid2D, mpi_sph4grid3D  # adjust to your package path

RNG = np.random.default_rng(42 + mpi.rank)  # different seed per rank
N = 100
BOXSIZE = 1.0

# Shared random data — each rank holds its own particles
X2, Y2 = RNG.uniform(0, BOXSIZE, N), RNG.uniform(0, BOXSIZE, N)
X3, Y3, Z3 = (
    RNG.uniform(0, BOXSIZE, N),
    RNG.uniform(0, BOXSIZE, N),
    RNG.uniform(0, BOXSIZE, N),
)
F2 = RNG.standard_normal(N)
F3 = RNG.standard_normal(N)


# ===========================================================================
# mpi_sph4grid2D
# ===========================================================================


class TestMpiSph4grid2DShape:

    def test_scalar_ngrid(self):
        result = mpi_sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=9, MPI=mpi)
        assert result.ndim == 2

    def test_list_ngrid(self):
        result = mpi_sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=[9, 6], MPI=mpi)
        assert result.ndim == 2

    def test_field_shape(self):
        result = mpi_sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=9, MPI=mpi, f=F2)
        assert result.ndim == 2

    def test_outputgrid_density_returns_tuple(self):
        dgrid, x2d, y2d = mpi_sph4grid2D(
            X2, Y2, boxsize=BOXSIZE, ngrid=9, MPI=mpi, outputgrid=True
        )
        assert dgrid.shape == x2d.shape == y2d.shape

    def test_outputgrid_field_returns_tuple(self):
        fgrid, x2d, y2d = mpi_sph4grid2D(
            X2, Y2, boxsize=BOXSIZE, ngrid=9, MPI=mpi, f=F2, outputgrid=True
        )
        assert fgrid.shape == x2d.shape == y2d.shape

    def test_list_boxsize(self):
        result = mpi_sph4grid2D(X2, Y2, boxsize=[1.0, 2.0], ngrid=9, MPI=mpi)
        assert result.ndim == 2

    def test_list_ngrid_list_boxsize(self):
        result = mpi_sph4grid2D(X2, Y2, boxsize=[1.0, 2.0], ngrid=[9, 6], MPI=mpi)
        assert result.ndim == 2


class TestMpiSph4grid2DValues:

    def test_density_is_finite(self):
        result = mpi_sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=9, MPI=mpi)
        assert np.all(np.isfinite(result))

    def test_density_is_nonnegative(self):
        result = mpi_sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=9, MPI=mpi)
        assert np.all(result >= 0)

    def test_field_is_finite(self):
        result = mpi_sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=9, MPI=mpi, f=F2)
        assert np.all(np.isfinite(result))

    def test_explicit_mass(self):
        mass = RNG.uniform(0.5, 1.5, N)
        result = mpi_sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=9, MPI=mpi, mass=mass)
        assert np.all(np.isfinite(result))

    def test_explicit_mass_and_field(self):
        mass = RNG.uniform(0.5, 1.5, N)
        result = mpi_sph4grid2D(
            X2, Y2, boxsize=BOXSIZE, ngrid=9, MPI=mpi, mass=mass, f=F2
        )
        assert np.all(np.isfinite(result))

    def test_default_mass_matches_unit_mass(self):
        r_default = mpi_sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=9, MPI=mpi)
        r_ones = mpi_sph4grid2D(
            X2, Y2, boxsize=BOXSIZE, ngrid=9, MPI=mpi, mass=np.ones(N)
        )
        np.testing.assert_array_equal(r_default, r_ones)

    def test_subsampling_does_not_change_shape(self):
        r1 = mpi_sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=9, MPI=mpi, subsampling=1)
        r2 = mpi_sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=9, MPI=mpi, subsampling=2)
        assert r1.shape == r2.shape


class TestMpiSph4grid2DOptions:

    def test_periodic_false(self):
        result = mpi_sph4grid2D(
            X2, Y2, boxsize=BOXSIZE, ngrid=9, MPI=mpi, periodic=False
        )
        assert result.ndim == 2

    def test_custom_k(self):
        result = mpi_sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=9, MPI=mpi, k=10)
        assert result.ndim == 2

    def test_buffer_length(self):
        result = mpi_sph4grid2D(
            X2, Y2, boxsize=BOXSIZE, ngrid=9, MPI=mpi, buffer_length=0.1
        )
        assert result.ndim == 2


class TestMpiSph4grid2DListParams:

    # --- boxsize ---
    def test_list_boxsize_symmetric(self):
        result = mpi_sph4grid2D(X2, Y2, boxsize=[1.0, 1.0], ngrid=9, MPI=mpi)
        assert result.ndim == 2

    def test_list_boxsize_asymmetric(self):
        result = mpi_sph4grid2D(X2, Y2, boxsize=[1.0, 2.0], ngrid=9, MPI=mpi)
        assert result.ndim == 2

    def test_list_boxsize_matches_scalar(self):
        r_scalar = mpi_sph4grid2D(X2, Y2, boxsize=1.0, ngrid=9, MPI=mpi)
        r_list = mpi_sph4grid2D(X2, Y2, boxsize=[1.0, 1.0], ngrid=9, MPI=mpi)
        np.testing.assert_array_equal(r_scalar, r_list)

    # --- ngrid ---
    def test_list_ngrid_symmetric(self):
        result = mpi_sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=[9, 9], MPI=mpi)
        assert result.ndim == 2

    def test_list_ngrid_asymmetric(self):
        result = mpi_sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=[9, 6], MPI=mpi)
        assert result.ndim == 2

    def test_list_ngrid_matches_scalar(self):
        r_scalar = mpi_sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=9, MPI=mpi)
        r_list = mpi_sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=[9, 9], MPI=mpi)
        np.testing.assert_array_equal(r_scalar, r_list)

    # --- origin ---
    def test_list_origin_zeros(self):
        result = mpi_sph4grid2D(
            X2, Y2, boxsize=BOXSIZE, ngrid=9, MPI=mpi, origin=[0.0, 0.0]
        )
        assert result.ndim == 2

    def test_list_origin_asymmetric(self):
        result = mpi_sph4grid2D(
            X2, Y2, boxsize=BOXSIZE, ngrid=9, MPI=mpi, origin=[0.1, 0.2]
        )
        assert result.ndim == 2

    def test_list_origin_matches_scalar(self):
        r_scalar = mpi_sph4grid2D(X2, Y2, boxsize=BOXSIZE, ngrid=9, MPI=mpi, origin=0.0)
        r_list = mpi_sph4grid2D(
            X2, Y2, boxsize=BOXSIZE, ngrid=9, MPI=mpi, origin=[0.0, 0.0]
        )
        np.testing.assert_array_equal(r_scalar, r_list)

    # --- all list ---
    def test_all_list_params(self):
        result = mpi_sph4grid2D(
            X2, Y2, boxsize=[1.0, 2.0], ngrid=[9, 6], MPI=mpi, origin=[0.1, 0.2]
        )
        assert result.ndim == 2
        assert np.all(np.isfinite(result))


# ===========================================================================
# mpi_sph4grid3D
# ===========================================================================


class TestMpiSph4grid3DShape:

    def test_scalar_ngrid(self):
        result = mpi_sph4grid3D(X3, Y3, Z3, boxsize=BOXSIZE, ngrid=9, MPI=mpi)
        assert result.ndim == 3

    def test_list_ngrid(self):
        result = mpi_sph4grid3D(X3, Y3, Z3, boxsize=BOXSIZE, ngrid=[9, 6, 6], MPI=mpi)
        assert result.ndim == 3

    def test_field_shape(self):
        result = mpi_sph4grid3D(X3, Y3, Z3, boxsize=BOXSIZE, ngrid=9, MPI=mpi, f=F3)
        assert result.ndim == 3

    def test_outputgrid_density_returns_tuple(self):
        dgrid, x3d, y3d, z3d = mpi_sph4grid3D(
            X3, Y3, Z3, boxsize=BOXSIZE, ngrid=9, MPI=mpi, outputgrid=True
        )
        assert dgrid.shape == x3d.shape == y3d.shape == z3d.shape

    def test_outputgrid_field_returns_tuple(self):
        fgrid, x3d, y3d, z3d = mpi_sph4grid3D(
            X3, Y3, Z3, boxsize=BOXSIZE, ngrid=9, MPI=mpi, f=F3, outputgrid=True
        )
        assert fgrid.shape == x3d.shape == y3d.shape == z3d.shape

    def test_list_boxsize(self):
        result = mpi_sph4grid3D(X3, Y3, Z3, boxsize=[1.0, 2.0, 3.0], ngrid=9, MPI=mpi)
        assert result.ndim == 3

    def test_list_ngrid_list_boxsize(self):
        result = mpi_sph4grid3D(
            X3, Y3, Z3, boxsize=[1.0, 2.0, 3.0], ngrid=[9, 6, 6], MPI=mpi
        )
        assert result.ndim == 3


class TestMpiSph4grid3DValues:

    def test_density_is_finite(self):
        result = mpi_sph4grid3D(X3, Y3, Z3, boxsize=BOXSIZE, ngrid=9, MPI=mpi)
        assert np.all(np.isfinite(result))

    def test_density_is_nonnegative(self):
        result = mpi_sph4grid3D(X3, Y3, Z3, boxsize=BOXSIZE, ngrid=9, MPI=mpi)
        assert np.all(result >= 0)

    def test_field_is_finite(self):
        result = mpi_sph4grid3D(X3, Y3, Z3, boxsize=BOXSIZE, ngrid=9, MPI=mpi, f=F3)
        assert np.all(np.isfinite(result))

    def test_explicit_mass(self):
        mass = RNG.uniform(0.5, 1.5, N)
        result = mpi_sph4grid3D(
            X3, Y3, Z3, boxsize=BOXSIZE, ngrid=9, MPI=mpi, mass=mass
        )
        assert np.all(np.isfinite(result))

    def test_explicit_mass_and_field(self):
        mass = RNG.uniform(0.5, 1.5, N)
        result = mpi_sph4grid3D(
            X3, Y3, Z3, boxsize=BOXSIZE, ngrid=9, MPI=mpi, mass=mass, f=F2
        )
        assert np.all(np.isfinite(result))

    def test_default_mass_matches_unit_mass(self):
        r_default = mpi_sph4grid3D(X3, Y3, Z3, boxsize=BOXSIZE, ngrid=9, MPI=mpi)
        r_ones = mpi_sph4grid3D(
            X3, Y3, Z3, boxsize=BOXSIZE, ngrid=9, MPI=mpi, mass=np.ones(N)
        )
        np.testing.assert_array_equal(r_default, r_ones)

    def test_subsampling_does_not_change_shape(self):
        r1 = mpi_sph4grid3D(
            X3, Y3, Z3, boxsize=BOXSIZE, ngrid=9, MPI=mpi, subsampling=1
        )
        r2 = mpi_sph4grid3D(
            X3, Y3, Z3, boxsize=BOXSIZE, ngrid=9, MPI=mpi, subsampling=2
        )
        assert r1.shape == r2.shape


class TestMpiSph4grid3DOptions:

    def test_periodic_false(self):
        result = mpi_sph4grid3D(
            X3, Y3, Z3, boxsize=BOXSIZE, ngrid=9, MPI=mpi, periodic=False
        )
        assert result.ndim == 3

    def test_custom_k(self):
        result = mpi_sph4grid3D(X3, Y3, Z3, boxsize=BOXSIZE, ngrid=9, MPI=mpi, k=10)
        assert result.ndim == 3

    def test_buffer_length(self):
        result = mpi_sph4grid3D(
            X3, Y3, Z3, boxsize=BOXSIZE, ngrid=9, MPI=mpi, buffer_length=0.1
        )
        assert result.ndim == 3


class TestMpiSph4grid3DListParams:

    # --- boxsize ---
    def test_list_boxsize_symmetric(self):
        result = mpi_sph4grid3D(X3, Y3, Z3, boxsize=[1.0, 1.0, 1.0], ngrid=9, MPI=mpi)
        assert result.ndim == 3

    def test_list_boxsize_asymmetric(self):
        result = mpi_sph4grid3D(X3, Y3, Z3, boxsize=[1.0, 2.0, 3.0], ngrid=9, MPI=mpi)
        assert result.ndim == 3

    def test_list_boxsize_matches_scalar(self):
        r_scalar = mpi_sph4grid3D(X3, Y3, Z3, boxsize=1.0, ngrid=9, MPI=mpi)
        r_list = mpi_sph4grid3D(X3, Y3, Z3, boxsize=[1.0, 1.0, 1.0], ngrid=9, MPI=mpi)
        np.testing.assert_array_equal(r_scalar, r_list)

    # --- ngrid ---
    def test_list_ngrid_symmetric(self):
        result = mpi_sph4grid3D(X3, Y3, Z3, boxsize=BOXSIZE, ngrid=[9, 9, 9], MPI=mpi)
        assert result.ndim == 3

    def test_list_ngrid_asymmetric(self):
        result = mpi_sph4grid3D(X3, Y3, Z3, boxsize=BOXSIZE, ngrid=[9, 6, 6], MPI=mpi)
        assert result.ndim == 3

    def test_list_ngrid_matches_scalar(self):
        r_scalar = mpi_sph4grid3D(X3, Y3, Z3, boxsize=BOXSIZE, ngrid=9, MPI=mpi)
        r_list = mpi_sph4grid3D(X3, Y3, Z3, boxsize=BOXSIZE, ngrid=[9, 9, 9], MPI=mpi)
        np.testing.assert_array_equal(r_scalar, r_list)

    # --- origin ---
    def test_list_origin_zeros(self):
        result = mpi_sph4grid3D(
            X3, Y3, Z3, boxsize=BOXSIZE, ngrid=9, MPI=mpi, origin=[0.0, 0.0, 0.0]
        )
        assert result.ndim == 3

    def test_list_origin_asymmetric(self):
        result = mpi_sph4grid3D(
            X3, Y3, Z3, boxsize=BOXSIZE, ngrid=9, MPI=mpi, origin=[0.1, 0.2, 0.3]
        )
        assert result.ndim == 3

    def test_list_origin_matches_scalar(self):
        r_scalar = mpi_sph4grid3D(
            X3, Y3, Z3, boxsize=BOXSIZE, ngrid=9, MPI=mpi, origin=0.0
        )
        r_list = mpi_sph4grid3D(
            X3, Y3, Z3, boxsize=BOXSIZE, ngrid=9, MPI=mpi, origin=[0.0, 0.0, 0.0]
        )
        np.testing.assert_array_equal(r_scalar, r_list)

    # --- all list ---
    def test_all_list_params(self):
        result = mpi_sph4grid3D(
            X3,
            Y3,
            Z3,
            boxsize=[1.0, 2.0, 3.0],
            ngrid=[9, 6, 6],
            MPI=mpi,
            origin=[0.1, 0.2, 0.3],
        )
        assert result.ndim == 3
        assert np.all(np.isfinite(result))
