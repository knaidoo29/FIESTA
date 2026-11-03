import numpy as np
import pytest

pytest.importorskip("mpi4py")
import shift
from shift.mpiutils import MPI

from fiesta.dtfe import (
    mpi_delaunay_density4grid2D,
    mpi_delaunay_field4grid2D,
    mpi_delaunay_density4grid3D,
    mpi_delaunay_field4grid3D,
)


def assert_mpi_local_shape(grid, ngrid, mpi):
    assert grid.ndim == len(ngrid)
    assert grid.shape[1:] == tuple(ngrid[1:])
    expected_base = ngrid[0] // mpi.size
    assert grid.shape[0] in (expected_base, expected_base + 1)


def local_partition_mask(x, xmin, xmax, rank, size):
    if rank == size - 1:
        return (x >= xmin) & (x <= xmax)
    return (x >= xmin) & (x < xmax)


def test_mpi_delaunay_density4grid2D_outputgrid_returns_local_grid_and_coords():
    mpi = MPI()
    rng = np.random.default_rng(1234)

    x = rng.random(24)
    y = rng.random(24)
    xedges, _ = shift.cart.mpi_grid1D(1.0, 8, mpi, origin=0.0)
    local_mask = local_partition_mask(x, xedges[0], xedges[-1], mpi.rank, mpi.size)
    x_local = x[local_mask]
    y_local = y[local_mask]

    dens, x2d, y2d = mpi_delaunay_density4grid2D(
        x_local,
        y_local,
        1.0,
        8,
        mpi,
        origin=0.0,
        partition=1,
        periodic=[False, False],
        fbuffer=0.25,
        subsampling=1,
        outputgrid=True,
    )

    assert_mpi_local_shape(dens, [8, 8], mpi)
    assert np.all(np.isfinite(dens))
    assert x2d.ndim == 1 and y2d.ndim == 1
    assert x2d.shape == y2d.shape
    assert x2d.size == dens.size
    assert np.all(np.isfinite(x2d))
    assert np.all(np.isfinite(y2d))

    dens, x2d, y2d = mpi_delaunay_density4grid2D(
        x_local,
        y_local,
        [1.0, 1.0],
        [8, 8],
        mpi,
        origin=[0.0, 0.0],
        partition=1,
        periodic=[False, False],
        fbuffer=0.25,
        subsampling=[1, 1],
        outputgrid=True,
        mass=np.ones_like(x_local),
    )

    assert_mpi_local_shape(dens, [8, 8], mpi)
    assert np.all(np.isfinite(dens))
    assert x2d.ndim == 1 and y2d.ndim == 1
    assert x2d.shape == y2d.shape
    assert x2d.size == dens.size
    assert np.all(np.isfinite(x2d))
    assert np.all(np.isfinite(y2d))


def test_mpi_delaunay_field4grid2D_outputgrid_returns_local_grid_and_coords():
    mpi = MPI()
    rng = np.random.default_rng(2345)

    x = rng.random(24)
    y = rng.random(24)
    f = rng.normal(size=24)
    xedges, _ = shift.cart.mpi_grid1D(1.0, 8, mpi, origin=0.0)
    local_mask = local_partition_mask(x, xedges[0], xedges[-1], mpi.rank, mpi.size)
    x_local = x[local_mask]
    y_local = y[local_mask]
    f_local = f[local_mask]

    field, x2d, y2d = mpi_delaunay_field4grid2D(
        x_local,
        y_local,
        f_local,
        1.0,
        8,
        mpi,
        partition=1,
        periodic=[False, False],
        fbuffer=0.25,
        subsampling=1,
        outputgrid=True,
    )

    field, x2d, y2d = mpi_delaunay_field4grid2D(
        x_local,
        y_local,
        f_local,
        [1.0, 1.0],
        [8, 8],
        mpi,
        mass=np.ones_like(x_local),
        origin=[0.0, 0.0],
        partition=1,
        periodic=[False, False],
        fbuffer=0.25,
        subsampling=[1, 1],
        outputgrid=True,
    )

    assert_mpi_local_shape(field, [8, 8], mpi)
    assert np.all(np.isfinite(field))
    assert x2d.ndim == 1 and y2d.ndim == 1
    assert x2d.shape == y2d.shape
    assert x2d.size == field.size
    assert np.all(np.isfinite(x2d))
    assert np.all(np.isfinite(y2d))


def test_mpi_delaunay_density4grid3D_outputgrid_returns_local_grid_and_coords():
    mpi = MPI()
    rng = np.random.default_rng(3456)

    x = rng.random(32)
    y = rng.random(32)
    z = rng.random(32)
    xedges, _ = shift.cart.mpi_grid1D(1.0, 6, mpi, origin=0.0)
    local_mask = local_partition_mask(x, xedges[0], xedges[-1], mpi.rank, mpi.size)
    x_local = x[local_mask]
    y_local = y[local_mask]
    z_local = z[local_mask]

    dens, x3d, y3d, z3d = mpi_delaunay_density4grid3D(
        x_local,
        y_local,
        z_local,
        1.0,
        6,
        mpi,
        partition=1,
        periodic=[False, False, False],
        fbuffer=0.25,
        subsampling=1,
        outputgrid=True,
    )

    dens, x3d, y3d, z3d = mpi_delaunay_density4grid3D(
        x_local,
        y_local,
        z_local,
        [1.0, 1.0, 1.0],
        [6, 6, 6],
        mpi,
        origin=[0.0, 0.0, 0.0],
        mass=np.ones_like(x_local),
        partition=1,
        periodic=[False, False, False],
        fbuffer=0.25,
        subsampling=[1, 1, 1],
        outputgrid=True,
    )

    assert_mpi_local_shape(dens, [6, 6, 6], mpi)
    assert np.all(np.isfinite(dens))
    assert x3d.ndim == 1 and y3d.ndim == 1 and z3d.ndim == 1
    assert x3d.shape == y3d.shape == z3d.shape
    assert x3d.size == dens.size
    assert np.all(np.isfinite(x3d))
    assert np.all(np.isfinite(y3d))
    assert np.all(np.isfinite(z3d))


def test_mpi_delaunay_field4grid3D_outputgrid_returns_local_grid_and_coords():
    mpi = MPI()
    rng = np.random.default_rng(4567)

    x = rng.random(32)
    y = rng.random(32)
    z = rng.random(32)
    f = rng.normal(size=32)
    xedges, _ = shift.cart.mpi_grid1D(1.0, 6, mpi, origin=0.0)
    local_mask = local_partition_mask(x, xedges[0], xedges[-1], mpi.rank, mpi.size)
    x_local = x[local_mask]
    y_local = y[local_mask]
    z_local = z[local_mask]
    f_local = f[local_mask]

    field, x3d, y3d, z3d = mpi_delaunay_field4grid3D(
        x_local,
        y_local,
        z_local,
        f_local,
        1.0,
        6,
        mpi,
        partition=1,
        periodic=[False, False, False],
        fbuffer=0.25,
        subsampling=1,
        outputgrid=True,
    )

    field, x3d, y3d, z3d = mpi_delaunay_field4grid3D(
        x_local,
        y_local,
        z_local,
        f_local,
        [1.0, 1.0, 1.0],
        [6, 6, 6],
        mpi,
        mass=np.ones_like(x_local),
        origin=[0.0, 0.0, 0.0],
        partition=1,
        periodic=[False, False, False],
        fbuffer=0.25,
        subsampling=[1, 1, 1],
        outputgrid=True,
    )

    assert_mpi_local_shape(field, [6, 6, 6], mpi)
    assert np.all(np.isfinite(field))
    assert x3d.ndim == 1 and y3d.ndim == 1 and z3d.ndim == 1
    assert x3d.shape == y3d.shape == z3d.shape
    assert x3d.size == field.size
    assert np.all(np.isfinite(x3d))
    assert np.all(np.isfinite(y3d))
    assert np.all(np.isfinite(z3d))


@pytest.mark.parametrize(
    "ndim, field_mode",
    [
        (2, False),
        (2, True),
        (3, False),
        (3, True),
    ],
)
def test_mpi_delaunay_root_only_input(ndim, field_mode):
    """Test particle redistribution when only rank 0 supplies data."""
    mpi = MPI()

    if mpi.size == 1:
        pytest.skip("Requires at least two MPI ranks.")

    rng = np.random.default_rng(1234)
    ngrid = max(8, 2 * mpi.size)

    if mpi.rank == 0:
        npart = 256 if ndim == 2 else 800
        x = rng.random(npart)
        y = rng.random(npart)
        z = rng.random(npart) if ndim == 3 else None
        f = np.ones(npart)
    else:
        x = None
        y = None
        z = None
        f = None

    kwargs = {
        "partition": 1,
        "periodic": [False] * ndim,
        "fbuffer": 0.25,
        "subsampling": 1,
        "outputgrid": True,
        "check_distributed": True,
    }

    if ndim == 2:
        if field_mode:
            result = mpi_delaunay_field4grid2D(
                x, y, f, 1.0, ngrid, mpi, **kwargs
            )
        else:
            result = mpi_delaunay_density4grid2D(
                x, y, 1.0, ngrid, mpi, **kwargs
            )
    else:
        if field_mode:
            result = mpi_delaunay_field4grid3D(
                x, y, z, f, 1.0, ngrid, mpi, **kwargs
            )
        else:
            result = mpi_delaunay_density4grid3D(
                x, y, z, 1.0, ngrid, mpi, **kwargs
            )

    grid, *coordinates = result

    assert_mpi_local_shape(grid, [ngrid] * ndim, mpi)
    assert np.all(np.isfinite(grid))
    assert len(coordinates) == ndim

    for coord in coordinates:
        assert coord.ndim == 1
        assert coord.size == grid.size
        assert np.all(np.isfinite(coord))

    # Ensure that the redistributed particles produce a nonzero field.
    local_nonzero = np.count_nonzero(grid)
    total_nonzero = mpi.Sum(local_nonzero)
    total_nonzero = mpi.bcast(total_nonzero)

    assert total_nonzero > 0
