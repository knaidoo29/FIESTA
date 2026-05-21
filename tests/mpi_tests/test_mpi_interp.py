import numpy as np
import pytest

import shift

pytest.importorskip("mpi4py")
from shift.mpiutils import MPI
mpi = MPI()

from fiesta import coords
from fiesta.interp import mpi_bilinear, mpi_trilinear

RNG = np.random.default_rng(42 + mpi.rank)  # different seed per rank
N = 100
BOXSIZE = 1.0
NGRID = 9  # divisible by 3 ranks

def make_fgrid2D(ngrid=NGRID):
    """Build a correctly-shaped distributed 2D fgrid for the local rank."""
    return shift.cart.mpi_grid2D(BOXSIZE, ngrid, mpi, origin=0.0)[0]


def make_fgrid3D(ngrid=NGRID):
    """Build a correctly-shaped distributed 3D fgrid for the local rank."""
    return shift.cart.mpi_grid3D(BOXSIZE, ngrid, mpi, origin=0.0)[0]

if mpi.rank == 0:
    # Random query points on each rank
    X2, Y2 = RNG.uniform(0, BOXSIZE, N), RNG.uniform(0, BOXSIZE, N)
    X3, Y3, Z3 = RNG.uniform(0, BOXSIZE, N), RNG.uniform(0, BOXSIZE, N), RNG.uniform(0, BOXSIZE, N)
else:
    X2 = None
    Y2 = None
    X3 = None
    Y3 = None
    Z3 = None

if mpi.rank == 0:
    data = coords.xy2points(X2, Y2)
else:
    data = None
data = coords.distribute_points_by_x(data, BOXSIZE, NGRID, origin=0.0, MPI=mpi)
X2, Y2 = coords.points2xy(data)

if mpi.rank == 0:
    data = coords.xyz2points(X3, Y3, Z3)
else:
    data = None
data = coords.distribute_points_by_x(data, BOXSIZE, NGRID, origin=0.0, MPI=mpi)
X3, Y3, Z3 = coords.points2xyz(data)

# ===========================================================================
# mpi_bilinear
# ===========================================================================

class TestMpiBilinearReturnStructure:

    def test_returns_three_tuple(self):
        fgrid = make_fgrid2D()
        result = mpi_bilinear(fgrid, NGRID, BOXSIZE, X2, Y2, mpi)
        assert isinstance(result, tuple) and len(result) == 3

    def test_x_y_f_same_length(self):
        fgrid = make_fgrid2D()
        x, y, f = mpi_bilinear(fgrid, NGRID, BOXSIZE, X2, Y2, mpi)
        assert len(x) == len(y) == len(f)

    def test_f_is_array(self):
        fgrid = make_fgrid2D()
        x, y, f = mpi_bilinear(fgrid, NGRID, BOXSIZE, X2, Y2, mpi)
        assert isinstance(f, np.ndarray)


class TestMpiBilinearValues:

    def test_f_is_finite_inside_box(self):
        fgrid = make_fgrid2D()
        x, y, f = mpi_bilinear(fgrid, ngrid=NGRID, boxsize=BOXSIZE, x=X2, y=Y2, MPI=mpi)
        assert np.all(np.isfinite(f))

    def test_fill_value_outside_box(self):
        fgrid = make_fgrid2D()
        fill = -999.0
        x_out = RNG.uniform(2.0, 3.0, N)   # outside [0, BOXSIZE]
        y_out = RNG.uniform(2.0, 3.0, N)
        x, y, f = mpi_bilinear(
            fgrid, ngrid=NGRID, boxsize=BOXSIZE, x=x_out, y=y_out,
            MPI=mpi, periodic=False, fill_value=fill,
        )
        assert np.all(f == fill)

    def test_periodic_false_runs(self):
        fgrid = make_fgrid2D()
        x, y, f = mpi_bilinear(
            fgrid, ngrid=NGRID, boxsize=BOXSIZE, x=X2, y=Y2, MPI=mpi, periodic=False
        )
        assert len(f) == len(x)


class TestMpiBilinearWrongShape:

    def test_wrong_fgrid_shape_returns_nan(self):
        bad_fgrid = RNG.uniform(0, 1, (2, 2))  # wrong shape
        x, y, f = mpi_bilinear(bad_fgrid, ngrid=NGRID, boxsize=BOXSIZE, x=X2, y=Y2, MPI=mpi)
        assert np.all(np.isnan([x, y, f]))


class TestMpiBilinearListParams:

    # --- boxsize ---
    def test_list_boxsize_symmetric(self):
        fgrid = make_fgrid2D()
        x, y, f = mpi_bilinear(fgrid, ngrid=NGRID, boxsize=[1.0, 1.0], x=X2, y=Y2, MPI=mpi)
        assert len(f) == len(x)

    def test_list_boxsize_asymmetric(self):
        fgrid = make_fgrid2D()
        x, y, f = mpi_bilinear(fgrid, ngrid=NGRID, boxsize=[1.0, 2.0], x=X2, y=Y2, MPI=mpi)
        assert len(f) == len(x)

    def test_list_boxsize_matches_scalar(self):
        fgrid = make_fgrid2D()
        x1, y1, f1 = mpi_bilinear(fgrid, ngrid=NGRID, boxsize=1.0,        x=X2, y=Y2, MPI=mpi)
        x2, y2, f2 = mpi_bilinear(fgrid, ngrid=NGRID, boxsize=[1.0, 1.0], x=X2, y=Y2, MPI=mpi)
        np.testing.assert_array_equal(f1, f2)

    # --- ngrid ---
    def test_list_ngrid_symmetric(self):
        fgrid = make_fgrid2D()
        x, y, f = mpi_bilinear(fgrid, ngrid=[NGRID, NGRID], boxsize=BOXSIZE, x=X2, y=Y2, MPI=mpi)
        assert len(f) == len(x)

    def test_list_ngrid_asymmetric(self):
        nygrid = 6
        fgrid = make_fgrid2D(ngrid=[NGRID, nygrid])
        x, y, f = mpi_bilinear(fgrid, ngrid=[NGRID, nygrid], boxsize=BOXSIZE, x=X2, y=Y2, MPI=mpi)
        assert len(f) == len(x)

    def test_list_ngrid_matches_scalar(self):
        fgrid = make_fgrid2D()
        x1, y1, f1 = mpi_bilinear(fgrid, ngrid=NGRID,            boxsize=BOXSIZE, x=X2, y=Y2, MPI=mpi)
        x2, y2, f2 = mpi_bilinear(fgrid, ngrid=[NGRID, NGRID],   boxsize=BOXSIZE, x=X2, y=Y2, MPI=mpi)
        np.testing.assert_array_equal(f1, f2)

    # --- origin ---
    def test_list_origin_zeros(self):
        fgrid = make_fgrid2D()
        x, y, f = mpi_bilinear(fgrid, ngrid=NGRID, boxsize=BOXSIZE, x=X2, y=Y2, MPI=mpi, origin=[0.0, 0.0])
        assert len(f) == len(x)

    def test_list_origin_matches_scalar(self):
        fgrid = make_fgrid2D()
        x1, y1, f1 = mpi_bilinear(fgrid, ngrid=NGRID, boxsize=BOXSIZE, x=X2, y=Y2, MPI=mpi, origin=0.0)
        x2, y2, f2 = mpi_bilinear(fgrid, ngrid=NGRID, boxsize=BOXSIZE, x=X2, y=Y2, MPI=mpi, origin=[0.0, 0.0])
        np.testing.assert_array_equal(f1, f2)

    # --- periodic ---
    def test_list_periodic_both_true(self):
        fgrid = make_fgrid2D()
        x, y, f = mpi_bilinear(fgrid, ngrid=NGRID, boxsize=BOXSIZE, x=X2, y=Y2, MPI=mpi, periodic=[True, True])
        assert len(f) == len(x)

    def test_list_periodic_mixed(self):
        fgrid = make_fgrid2D()
        x, y, f = mpi_bilinear(fgrid, ngrid=NGRID, boxsize=BOXSIZE, x=X2, y=Y2, MPI=mpi, periodic=[True, False])
        assert len(f) == len(x)

    def test_list_periodic_matches_scalar(self):
        fgrid = make_fgrid2D()
        x1, y1, f1 = mpi_bilinear(fgrid, ngrid=NGRID, boxsize=BOXSIZE, x=X2, y=Y2, MPI=mpi, periodic=True)
        x2, y2, f2 = mpi_bilinear(fgrid, ngrid=NGRID, boxsize=BOXSIZE, x=X2, y=Y2, MPI=mpi, periodic=[True, True])
        np.testing.assert_array_equal(f1, f2)

    # --- all list ---
    def test_all_list_params(self):
        nygrid = 6
        fgrid = make_fgrid2D(ngrid=[NGRID, nygrid])
        x, y, f = mpi_bilinear(
            fgrid, ngrid=[NGRID, nygrid], boxsize=[1.0, 2.0], x=X2, y=Y2,
            MPI=mpi, origin=[0.0, 0.0], periodic=[True, False],
        )
        assert len(f) == len(x)
        assert np.all(np.isfinite(f))

    def test_None_at_mpi_rank_no_0(self):
        if mpi.rank == 0:
            # Random query points on each rank
            _X2, _Y2 = RNG.uniform(0, BOXSIZE, N), RNG.uniform(0, BOXSIZE, N)
        else:
            _X2 = None
            _Y2 = None
        x, y, f = mpi_bilinear(make_fgrid2D(), ngrid=NGRID, boxsize=BOXSIZE, x=_X2, y=_Y2, MPI=mpi)
        assert len(f) == len(x)
        assert np.all(np.isfinite(f))

# ===========================================================================
# mpi_trilinear
# ===========================================================================

class TestMpiTrilinearReturnStructure:

    def test_returns_four_tuple(self):
        fgrid = make_fgrid3D()
        result = mpi_trilinear(fgrid, ngrid=NGRID, boxsize=BOXSIZE, x=X3, y=Y3, z=Z3, MPI=mpi)
        assert isinstance(result, tuple) and len(result) == 4

    def test_x_y_z_f_same_length(self):
        fgrid = make_fgrid3D()
        x, y, z, f = mpi_trilinear(fgrid, ngrid=NGRID, boxsize=BOXSIZE, x=X3, y=Y3, z=Z3, MPI=mpi)
        assert len(x) == len(y) == len(z) == len(f)

    def test_f_is_array(self):
        fgrid = make_fgrid3D()
        x, y, z, f = mpi_trilinear(fgrid, ngrid=NGRID, boxsize=BOXSIZE, x=X3, y=Y3, z=Z3, MPI=mpi)
        assert isinstance(f, np.ndarray)


class TestMpiTrilinearValues:

    def test_f_is_finite_inside_box(self):
        fgrid = make_fgrid3D()
        x, y, z, f = mpi_trilinear(fgrid, ngrid=NGRID, boxsize=BOXSIZE, x=X3, y=Y3, z=Z3, MPI=mpi)
        assert np.all(np.isfinite(f))

    def test_fill_value_outside_box(self):
        fgrid = make_fgrid3D()
        fill = -999.0
        x_out = RNG.uniform(2.0, 3.0, N)
        y_out = RNG.uniform(2.0, 3.0, N)
        z_out = RNG.uniform(2.0, 3.0, N)
        x, y, z, f = mpi_trilinear(
            fgrid, ngrid=NGRID, boxsize=BOXSIZE, x=x_out, y=y_out, z=z_out,
            MPI=mpi, periodic=False, fill_value=fill,
        )
        assert np.all(f == fill)

    def test_periodic_false_runs(self):
        fgrid = make_fgrid3D()
        x, y, z, f = mpi_trilinear(
            fgrid, ngrid=NGRID, boxsize=BOXSIZE, x=X3, y=Y3, z=Z3, MPI=mpi, periodic=False
        )
        assert len(f) == len(x)


class TestMpiTrilinearWrongShape:

    def test_wrong_fgrid_shape_returns_nan(self):
        bad_fgrid = RNG.uniform(0, 1, (2, 2, 2))  # wrong shape
        x, y, z, f = mpi_trilinear(bad_fgrid, ngrid=NGRID, boxsize=BOXSIZE, x=X3, y=Y3, z=Z3, MPI=mpi)
        assert np.all(np.isnan([x, y, z, f]))


class TestMpiTrilinearListParams:

    # --- boxsize ---
    def test_list_boxsize_symmetric(self):
        fgrid = make_fgrid3D()
        x, y, z, f = mpi_trilinear(fgrid, ngrid=NGRID, boxsize=[1.0, 1.0, 1.0], x=X3, y=Y3, z=Z3, MPI=mpi)
        assert len(f) == len(x)

    def test_list_boxsize_asymmetric(self):
        fgrid = make_fgrid3D()
        x, y, z, f = mpi_trilinear(fgrid, ngrid=NGRID, boxsize=[1.0, 2.0, 3.0], x=X3, y=Y3, z=Z3, MPI=mpi)
        assert len(f) == len(x)

    def test_list_boxsize_matches_scalar(self):
        fgrid = make_fgrid3D()
        x1, y1, z1, f1 = mpi_trilinear(fgrid, ngrid=NGRID, boxsize=1.0,             x=X3, y=Y3, z=Z3, MPI=mpi)
        x2, y2, z2, f2 = mpi_trilinear(fgrid, ngrid=NGRID, boxsize=[1.0, 1.0, 1.0], x=X3, y=Y3, z=Z3, MPI=mpi)
        np.testing.assert_array_equal(f1, f2)

    # --- ngrid ---
    def test_list_ngrid_symmetric(self):
        fgrid = make_fgrid3D()
        x, y, z, f = mpi_trilinear(fgrid, ngrid=[NGRID, NGRID, NGRID], boxsize=BOXSIZE, x=X3, y=Y3, z=Z3, MPI=mpi)
        assert len(f) == len(x)

    def test_list_ngrid_asymmetric(self):
        nygrid, nzgrid = 6, 6
        fgrid = make_fgrid3D(ngrid=[NGRID, nygrid, nzgrid])
        x, y, z, f = mpi_trilinear(fgrid, ngrid=[NGRID, nygrid, nzgrid], boxsize=BOXSIZE, x=X3, y=Y3, z=Z3, MPI=mpi)
        assert len(f) == len(x)

    def test_list_ngrid_matches_scalar(self):
        fgrid = make_fgrid3D()
        x1, y1, z1, f1 = mpi_trilinear(fgrid, ngrid=NGRID,                  boxsize=BOXSIZE, x=X3, y=Y3, z=Z3, MPI=mpi)
        x2, y2, z2, f2 = mpi_trilinear(fgrid, ngrid=[NGRID, NGRID, NGRID],  boxsize=BOXSIZE, x=X3, y=Y3, z=Z3, MPI=mpi)
        np.testing.assert_array_equal(f1, f2)

    # --- origin ---
    def test_list_origin_zeros(self):
        fgrid = make_fgrid3D()
        x, y, z, f = mpi_trilinear(fgrid, ngrid=NGRID, boxsize=BOXSIZE, x=X3, y=Y3, z=Z3, MPI=mpi, origin=[0.0, 0.0, 0.0])
        assert len(f) == len(x)

    def test_list_origin_matches_scalar(self):
        fgrid = make_fgrid3D()
        x1, y1, z1, f1 = mpi_trilinear(fgrid, ngrid=NGRID, boxsize=BOXSIZE, x=X3, y=Y3, z=Z3, MPI=mpi, origin=0.0)
        x2, y2, z2, f2 = mpi_trilinear(fgrid, ngrid=NGRID, boxsize=BOXSIZE, x=X3, y=Y3, z=Z3, MPI=mpi, origin=[0.0, 0.0, 0.0])
        np.testing.assert_array_equal(f1, f2)

    # --- periodic ---
    def test_list_periodic_both_true(self):
        fgrid = make_fgrid3D()
        x, y, z, f = mpi_trilinear(fgrid, ngrid=NGRID, boxsize=BOXSIZE, x=X3, y=Y3, z=Z3, MPI=mpi, periodic=[True, True, True])
        assert len(f) == len(x)

    def test_list_periodic_mixed(self):
        fgrid = make_fgrid3D()
        x, y, z, f = mpi_trilinear(fgrid, ngrid=NGRID, boxsize=BOXSIZE, x=X3, y=Y3, z=Z3, MPI=mpi, periodic=[True, False, True])
        assert len(f) == len(x)

    def test_list_periodic_matches_scalar(self):
        fgrid = make_fgrid3D()
        x1, y1, z1, f1 = mpi_trilinear(fgrid, ngrid=NGRID, boxsize=BOXSIZE, x=X3, y=Y3, z=Z3, MPI=mpi, periodic=True)
        x2, y2, z2, f2 = mpi_trilinear(fgrid, ngrid=NGRID, boxsize=BOXSIZE, x=X3, y=Y3, z=Z3, MPI=mpi, periodic=[True, True, True])
        np.testing.assert_array_equal(f1, f2)

    # --- all list ---
    def test_all_list_params(self):
        nygrid, nzgrid = 6, 6
        fgrid = make_fgrid3D(ngrid=[NGRID, nygrid, nzgrid])
        x, y, z, f = mpi_trilinear(
            fgrid, ngrid=[NGRID, nygrid, nzgrid], boxsize=[1.0, 2.0, 3.0], x=X3, y=Y3, z=Z3,
            MPI=mpi, origin=[0.0, 0.0, 0.0], periodic=[True, False, True],
        )
        assert len(f) == len(x)
        assert np.all(np.isfinite(f))
    
    def test_None_at_mpi_rank_no_0(self):
        if mpi.rank == 0:
            # Random query points on each rank
            _X3, _Y3, _Z3 = RNG.uniform(0, BOXSIZE, N), RNG.uniform(0, BOXSIZE, N), RNG.uniform(0, BOXSIZE, N)
        else:
            _X3 = None
            _Y3 = None
            _Z3 = None
        x, y, z, f = mpi_trilinear(make_fgrid3D(), ngrid=NGRID, boxsize=BOXSIZE, x=_X3, y=_Y3, z=_Z3, MPI=mpi)
        assert len(f) == len(x)
        assert np.all(np.isfinite(f))