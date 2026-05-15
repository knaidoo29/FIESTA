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
from fiesta.p2g.mpi_part2grid import mpi_part2grid2D, mpi_part2grid3D


def assert_mpi_local_shape(grid, ngrid, mpi):
    assert grid.ndim == len(ngrid)
    assert grid.shape[1:] == tuple(ngrid[1:])
    expected_base = ngrid[0] // mpi.size
    assert grid.shape[0] in (expected_base, expected_base + 1)


def get_local_partition(x, xmin, xmax, rank, size):
    if rank == size - 1:
        return (x >= xmin) & (x <= xmax)
    return (x >= xmin) & (x < xmax)


def test_mpi_part2grid2D_ngp():
    mpi = MPI()
    boxsize = [1.0, 1.0]
    ngrid = [8, 8]

    x = np.linspace(0.05, 0.95, 8)
    y = np.linspace(0.05, 0.95, 8)
    x, y = np.meshgrid(x, y, indexing="ij")
    x = x.ravel()
    y = y.ravel()
    f = np.ones_like(x)

    xedges, _ = shift.cart.mpi_grid1D(boxsize[0], ngrid[0], mpi, origin=0.0)
    local_mask = get_local_partition(x, xedges[0], xedges[-1], mpi.rank, mpi.size)
    x_local = x[local_mask]
    y_local = y[local_mask]
    f_local = f[local_mask]

    grid = mpi_part2grid2D(
        x_local, y_local, f_local, boxsize, ngrid, mpi, method="NGP", periodic=True
    )

    assert_mpi_local_shape(grid, ngrid, mpi)
    assert np.all(np.isfinite(grid))
    dx = (xedges[-1] - xedges[0]) / len(grid)
    dy = boxsize[1] / ngrid[1]
    assert np.isclose(grid.sum() * dx * dy, float(len(x_local)))


def test_mpi_part2grid3D_ngp():
    mpi = MPI()
    boxsize = [1.0, 1.0, 1.0]
    ngrid = [6, 6, 6]

    coords = np.linspace(0.05, 0.95, 6)
    x, y, z = np.meshgrid(coords, coords, coords, indexing="ij")
    x = x.ravel()
    y = y.ravel()
    z = z.ravel()
    f = np.ones_like(x)

    xedges, _ = shift.cart.mpi_grid1D(boxsize[0], ngrid[0], mpi, origin=0.0)
    local_mask = get_local_partition(x, xedges[0], xedges[-1], mpi.rank, mpi.size)
    x_local = x[local_mask]
    y_local = y[local_mask]
    z_local = z[local_mask]
    f_local = f[local_mask]

    grid = mpi_part2grid3D(
        x_local,
        y_local,
        z_local,
        f_local,
        boxsize,
        ngrid,
        mpi,
        method="NGP",
        periodic=True,
    )

    assert_mpi_local_shape(grid, ngrid, mpi)
    assert np.all(np.isfinite(grid))
    dx = (xedges[-1] - xedges[0]) / len(grid)
    dy = boxsize[1] / ngrid[1]
    dz = boxsize[2] / ngrid[2]
    assert np.isclose(grid.sum() * dx * dy * dz, float(len(x_local)))


def test_mpi_delaunay_density4grid2D():
    mpi = MPI()
    rng = np.random.default_rng(1234)

    x = rng.random(24)
    y = rng.random(24)
    xedges, _ = shift.cart.mpi_grid1D(1.0, 8, mpi, origin=0.0)
    local_mask = get_local_partition(x, xedges[0], xedges[-1], mpi.rank, mpi.size)
    x_local = x[local_mask]
    y_local = y[local_mask]

    dens = mpi_delaunay_density4grid2D(
        x_local,
        y_local,
        [1.0, 1.0],
        [8, 8],
        mpi,
        partition=1,
        periodic=True,
        fbuffer=0.25,
        subsampling=1,
    )

    assert_mpi_local_shape(dens, [8, 8], mpi)
    assert np.all(np.isfinite(dens))


def test_mpi_delaunay_field4grid2D():
    mpi = MPI()
    rng = np.random.default_rng(2345)

    x = rng.random(24)
    y = rng.random(24)
    f = rng.normal(size=24)
    xedges, _ = shift.cart.mpi_grid1D(1.0, 8, mpi, origin=0.0)
    local_mask = get_local_partition(x, xedges[0], xedges[-1], mpi.rank, mpi.size)
    x_local = x[local_mask]
    y_local = y[local_mask]
    f_local = f[local_mask]

    field = mpi_delaunay_field4grid2D(
        x_local,
        y_local,
        f_local,
        [1.0, 1.0],
        [8, 8],
        mpi,
        partition=1,
        periodic=True,
        fbuffer=0.25,
        subsampling=1,
    )

    assert_mpi_local_shape(field, [8, 8], mpi)
    assert np.all(np.isfinite(field))


def test_mpi_delaunay_density4grid3D():
    mpi = MPI()
    rng = np.random.default_rng(3456)

    x = rng.random(32)
    y = rng.random(32)
    z = rng.random(32)
    xedges, _ = shift.cart.mpi_grid1D(1.0, 6, mpi, origin=0.0)
    local_mask = get_local_partition(x, xedges[0], xedges[-1], mpi.rank, mpi.size)
    x_local = x[local_mask]
    y_local = y[local_mask]
    z_local = z[local_mask]

    dens = mpi_delaunay_density4grid3D(
        x_local,
        y_local,
        z_local,
        [1.0, 1.0, 1.0],
        [6, 6, 6],
        mpi,
        partition=1,
        periodic=True,
        fbuffer=0.25,
        subsampling=1,
    )

    assert_mpi_local_shape(dens, [6, 6, 6], mpi)
    assert np.all(np.isfinite(dens))


def test_mpi_delaunay_field4grid3D():
    mpi = MPI()
    rng = np.random.default_rng(4567)

    x = rng.random(32)
    y = rng.random(32)
    z = rng.random(32)
    f = rng.normal(size=32)
    xedges, _ = shift.cart.mpi_grid1D(1.0, 6, mpi, origin=0.0)
    local_mask = get_local_partition(x, xedges[0], xedges[-1], mpi.rank, mpi.size)
    x_local = x[local_mask]
    y_local = y[local_mask]
    z_local = z[local_mask]
    f_local = f[local_mask]

    field = mpi_delaunay_field4grid3D(
        x_local,
        y_local,
        z_local,
        f_local,
        [1.0, 1.0, 1.0],
        [6, 6, 6],
        mpi,
        partition=1,
        periodic=True,
        fbuffer=0.25,
        subsampling=1,
    )

    assert_mpi_local_shape(field, [6, 6, 6], mpi)
    assert np.all(np.isfinite(field))
