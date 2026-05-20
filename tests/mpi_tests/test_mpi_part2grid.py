import numpy as np
import pytest

pytest.importorskip("mpi4py")
import shift
from shift.mpiutils import MPI

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


@pytest.mark.parametrize("method", ["NGP", "CIC", "TSC", "PCS"])
@pytest.mark.parametrize("periodic", [True, False, [True, False], [False, True]])
def test_mpi_part2grid2d_returns_valid_local_grid(method, periodic):
    mpi = MPI()
    boxsize = 1.0
    ngrid = 8
    origin = 0.0

    x = np.linspace(0.05, 0.95, 24)
    y = np.linspace(0.05, 0.95, 12)
    x2d, y2d = np.meshgrid(x, y, indexing="ij")
    x_flat = x2d.ravel()
    y_flat = y2d.ravel()
    f = np.ones_like(x_flat)

    xedges, _ = shift.cart.mpi_grid1D(boxsize, ngrid, mpi, origin=origin)
    local_mask = get_local_partition(x_flat, xedges[0], xedges[-1], mpi.rank, mpi.size)
    x_local = x_flat[local_mask]
    y_local = y_flat[local_mask]
    f_local = f[local_mask]

    grid = mpi_part2grid2D(
        x_local,
        y_local,
        f_local,
        boxsize,
        ngrid,
        mpi,
        method=method,
        periodic=periodic,
        origin=origin,
    )

    if method == "NGP":
        dx = xedges[1] - xedges[0]
        dy = boxsize / ngrid
        assert np.isclose(
            grid.sum() * dx * dy, float(len(f_local)), atol=1e-8, rtol=1e-8
        )
    else:
        assert grid.sum() > 0


@pytest.mark.parametrize("method", ["NGP", "CIC", "TSC", "PCS"])
@pytest.mark.parametrize("periodic", [True, False, [True, False], [False, True]])
def test_mpi_part2grid2d_returns_valid_local_grid_list(method, periodic):
    mpi = MPI()
    boxsize = [1.0, 1.0]
    ngrid = [8, 8]
    origin = [0.0, 0.0]

    x = np.linspace(0.05, 0.95, 24)
    y = np.linspace(0.05, 0.95, 12)
    x2d, y2d = np.meshgrid(x, y, indexing="ij")
    x_flat = x2d.ravel()
    y_flat = y2d.ravel()
    f = np.ones_like(x_flat)

    xedges, _ = shift.cart.mpi_grid1D(boxsize[0], ngrid[0], mpi, origin=origin[0])
    local_mask = get_local_partition(x_flat, xedges[0], xedges[-1], mpi.rank, mpi.size)
    # x_local = x_flat[local_mask]
    # y_local = y_flat[local_mask]
    f_local2 = f[local_mask]

    if mpi.rank == 0:
        x_local = x_flat
        y_local = y_flat
        f_local = f
    else:
        x_local = None
        y_local = None
        f_local = None

    grid = mpi_part2grid2D(
        x_local,
        y_local,
        f_local,
        boxsize,
        ngrid,
        mpi,
        method=method,
        periodic=periodic,
        origin=origin,
    )

    assert_mpi_local_shape(grid, ngrid, mpi)
    assert np.all(np.isfinite(grid))
    if method == "NGP":
        dx = xedges[1] - xedges[0]
        dy = boxsize[1] / ngrid[1]
        assert np.isclose(
            grid.sum() * dx * dy, float(len(f_local2)), atol=1e-8, rtol=1e-8
        )
    else:
        assert grid.sum() > 0


@pytest.mark.parametrize("method", ["NGP", "CIC", "TSC", "PCS"])
@pytest.mark.parametrize(
    "periodic", [True, False, [True, False, True], [False, True, False]]
)
def test_mpi_part2grid3d_returns_valid_local_grid(method, periodic):
    mpi = MPI()
    boxsize = 1.0
    ngrid = 6
    origin = 0.0

    coords = np.linspace(0.05, 0.95, 18)
    x3d, y3d, z3d = np.meshgrid(coords, coords[:12], coords[:9], indexing="ij")
    x_flat = x3d.ravel()
    y_flat = y3d.ravel()
    z_flat = z3d.ravel()
    f = np.ones_like(x_flat)

    xedges, _ = shift.cart.mpi_grid1D(boxsize, ngrid, mpi, origin=origin)
    local_mask = get_local_partition(x_flat, xedges[0], xedges[-1], mpi.rank, mpi.size)
    # x_local = x_flat[local_mask]
    # y_local = y_flat[local_mask]
    # z_local = z_flat[local_mask]
    f_local2 = f[local_mask]

    if mpi.rank == 0:
        x_local = x_flat
        y_local = y_flat
        z_local = z_flat
        f_local = f
    else:
        x_local = None
        y_local = None
        z_local = None
        f_local = None
    
    grid = mpi_part2grid3D(
        x_local,
        y_local,
        z_local,
        f_local,
        boxsize,
        ngrid,
        mpi,
        method=method,
        periodic=periodic,
        origin=origin,
    )

    if method == "NGP":
        dx = xedges[1] - xedges[0]
        dy = boxsize / ngrid
        dz = boxsize / ngrid
        assert np.isclose(
            grid.sum() * dx * dy * dz, float(len(f_local2)), atol=1e-8, rtol=1e-8
        )
    else:
        assert grid.sum() > 0


@pytest.mark.parametrize("method", ["NGP", "CIC", "TSC", "PCS"])
@pytest.mark.parametrize(
    "periodic", [True, False, [True, False, True], [False, True, False]]
)
def test_mpi_part2grid3d_returns_valid_local_grid_list(method, periodic):
    mpi = MPI()
    boxsize = [1.0, 1.0, 1.0]
    ngrid = [6, 6, 6]
    origin = [0.0, 0.0, 0.0]

    coords = np.linspace(0.05, 0.95, 18)
    x3d, y3d, z3d = np.meshgrid(coords, coords[:12], coords[:9], indexing="ij")
    x_flat = x3d.ravel()
    y_flat = y3d.ravel()
    z_flat = z3d.ravel()
    f = np.ones_like(x_flat)

    xedges, _ = shift.cart.mpi_grid1D(boxsize[0], ngrid[0], mpi, origin=origin[0])
    local_mask = get_local_partition(x_flat, xedges[0], xedges[-1], mpi.rank, mpi.size)
    x_local = x_flat[local_mask]
    y_local = y_flat[local_mask]
    z_local = z_flat[local_mask]
    f_local = f[local_mask]

    grid = mpi_part2grid3D(
        x_local,
        y_local,
        z_local,
        f_local,
        boxsize,
        ngrid,
        mpi,
        method=method,
        periodic=periodic,
        origin=origin,
    )

    assert_mpi_local_shape(grid, ngrid, mpi)
    assert np.all(np.isfinite(grid))
    if method == "NGP":
        dx = xedges[1] - xedges[0]
        dy = boxsize[1] / ngrid[1]
        dz = boxsize[2] / ngrid[2]
        assert np.isclose(
            grid.sum() * dx * dy * dz, float(len(f_local)), atol=1e-8, rtol=1e-8
        )
    else:
        assert grid.sum() > 0
