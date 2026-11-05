import numpy as np

from fiesta import p2g
from fiesta.gridsph.mpi_gridsph import mpi_gridSPH2D, mpi_gridSPH3D
from shift.mpiutils import MPI


def uneven_grid_size(mpi, ngrid=20):
    """Choose a grid size that gives uneven, nonempty MPI partitions."""
    ngrid = max(ngrid, mpi.size + 1)

    if mpi.size > 1 and ngrid % mpi.size == 0:
        ngrid += 1

    return ngrid


def assert_mpi_grid(grid, mpi, ngrid, ndim):
    """Check the shape and finiteness of a distributed grid."""
    base, remainder = divmod(ngrid, mpi.size)
    expected_nx = base + (mpi.rank < remainder)

    expected_shape = (expected_nx,) + (ngrid,) * (ndim - 1)

    assert isinstance(grid, np.ndarray)
    assert grid.shape == expected_shape
    assert np.all(np.isfinite(grid))


def test_mpi_gridSPH2D():
    mpi = MPI()
    ngrid = uneven_grid_size(mpi)

    x = np.linspace(0.1, 0.9, 50)
    y = np.linspace(0.1, 0.9, 50)

    x2d, y2d = np.meshgrid(x, y, indexing="ij")

    dgrid = mpi_gridSPH2D(
        x2d.ravel(),
        y2d.ravel(),
        1.0,
        ngrid,
        mpi,
        minpart=1,
    )

    assert_mpi_grid(dgrid, mpi, ngrid, ndim=2)

    dgrid = mpi_gridSPH2D(
        x2d.ravel(),
        y2d.ravel(),
        [1.0, 1.0],
        [ngrid, ngrid],
        mpi,
        minpart=1,
        periodic=[False, True],
    )

    assert_mpi_grid(dgrid, mpi, ngrid, ndim=2)


def test_mpi_gridSPH2D_pregrid():
    mpi = MPI()
    ngrid = uneven_grid_size(mpi)

    x = np.linspace(0.1, 0.9, 50)
    y = np.linspace(0.1, 0.9, 50)

    x2d, y2d = np.meshgrid(x, y, indexing="ij")

    dgrid_per = p2g.mpi_part2grid2D(
        x2d.ravel(),
        y2d.ravel(),
        1.0,
        ngrid,
        mpi,
        f=np.ones(len(x2d.ravel())),
        method="NGP",
        periodic=[False, True],
        origin=0.0,
    )

    dgrid = mpi_gridSPH2D(
        x2d.ravel(),
        y2d.ravel(),
        1.0,
        ngrid,
        mpi,
        minpart=1,
        dgrid=dgrid_per,
    )

    assert_mpi_grid(dgrid, mpi, ngrid, ndim=2)


def test_mpi_gridSPH2D_field():
    mpi = MPI()
    ngrid = uneven_grid_size(mpi)

    x = np.linspace(0.1, 0.9, 50)
    y = np.linspace(0.1, 0.9, 50)

    x2d, y2d = np.meshgrid(x, y, indexing="ij")

    f = np.cos(2 * np.pi * x2d.ravel())

    fgrid = mpi_gridSPH2D(
        x2d.ravel(),
        y2d.ravel(),
        1.0,
        ngrid,
        mpi,
        minpart=1,
        f=f,
    )

    assert_mpi_grid(fgrid, mpi, ngrid, ndim=2)

    fgrid = mpi_gridSPH2D(
        x2d.ravel(),
        y2d.ravel(),
        1.0,
        ngrid,
        mpi,
        minpart=1,
        f=f,
        periodic=[False, True],
    )

    assert_mpi_grid(fgrid, mpi, ngrid, ndim=2)


def test_mpi_gridSPH2D_pregrid_field():
    mpi = MPI()
    ngrid = uneven_grid_size(mpi)

    x = np.linspace(0.1, 0.9, 50)
    y = np.linspace(0.1, 0.9, 50)

    x2d, y2d = np.meshgrid(x, y, indexing="ij")

    f = np.cos(2 * np.pi * x2d.ravel())

    dgrid_pre = p2g.mpi_part2grid2D(
        x2d.ravel(),
        y2d.ravel(),
        1.0,
        ngrid,
        mpi,
        f=np.ones(len(x2d.ravel())),
        method="NGP",
        periodic=[False, True],
        origin=0.0,
    )

    fgrid_pre = p2g.mpi_part2grid2D(
        x2d.ravel(),
        y2d.ravel(),
        1.0,
        ngrid,
        mpi,
        f=f,
        method="NGP",
        periodic=[False, True],
        origin=0.0,
    )

    fgrid = mpi_gridSPH2D(
        x2d.ravel(),
        y2d.ravel(),
        1.0,
        ngrid,
        mpi,
        minpart=1,
        f=f,
        dgrid=dgrid_pre,
        fgrid=fgrid_pre,
    )

    assert_mpi_grid(fgrid, mpi, ngrid, ndim=2)


def test_mpi_gridSPH3D():
    mpi = MPI()
    ngrid = uneven_grid_size(mpi)

    x = np.linspace(0.1, 0.9, 50)
    y = np.linspace(0.1, 0.9, 50)
    z = np.linspace(0.1, 0.9, 50)

    x3d, y3d, z3d = np.meshgrid(x, y, z, indexing="ij")

    dgrid = mpi_gridSPH3D(
        x3d.ravel(),
        y3d.ravel(),
        z3d.ravel(),
        1.0,
        ngrid,
        mpi,
        minpart=1,
    )

    assert_mpi_grid(dgrid, mpi, ngrid, ndim=3)

    dgrid = mpi_gridSPH3D(
        x3d.ravel(),
        y3d.ravel(),
        z3d.ravel(),
        [1.0, 1.0, 1.0],
        [ngrid, ngrid, ngrid],
        mpi,
        minpart=1,
        periodic=[False, True, True],
    )

    assert_mpi_grid(dgrid, mpi, ngrid, ndim=3)


def test_mpi_pregrid_gridSPH3D():
    mpi = MPI()
    ngrid = uneven_grid_size(mpi)

    x = np.linspace(0.1, 0.9, 50)
    y = np.linspace(0.1, 0.9, 50)
    z = np.linspace(0.1, 0.9, 50)

    x3d, y3d, z3d = np.meshgrid(x, y, z, indexing="ij")

    dgrid_per = p2g.mpi_part2grid3D(
        x3d.ravel(),
        y3d.ravel(),
        z3d.ravel(),
        1.0,
        ngrid,
        mpi,
        method="NGP",
        periodic=[False, True, True],
        origin=0.0,
    )

    dgrid = mpi_gridSPH3D(
        x3d.ravel(),
        y3d.ravel(),
        z3d.ravel(),
        1.0,
        ngrid,
        mpi,
        minpart=1,
        dgrid=dgrid_per,
    )

    assert_mpi_grid(dgrid, mpi, ngrid, ndim=3)


def test_mpi_gridSPH3D_field():
    mpi = MPI()
    ngrid = uneven_grid_size(mpi)

    x = np.linspace(0.1, 0.9, 50)
    y = np.linspace(0.1, 0.9, 50)
    z = np.linspace(0.1, 0.9, 50)

    x3d, y3d, z3d = np.meshgrid(x, y, z, indexing="ij")

    f = np.sin(2 * np.pi * x3d.ravel())

    fgrid = mpi_gridSPH3D(
        x3d.ravel(),
        y3d.ravel(),
        z3d.ravel(),
        1.0,
        ngrid,
        mpi,
        minpart=1,
        f=f,
    )

    assert_mpi_grid(fgrid, mpi, ngrid, ndim=3)

    fgrid = mpi_gridSPH3D(
        x3d.ravel(),
        y3d.ravel(),
        z3d.ravel(),
        1.0,
        ngrid,
        mpi,
        minpart=1,
        f=f,
        periodic=[False, True, True],
    )

    assert_mpi_grid(fgrid, mpi, ngrid, ndim=3)


def test_mpi_gridSPH3D_pregrid_field():
    mpi = MPI()
    ngrid = uneven_grid_size(mpi)

    x = np.linspace(0.1, 0.9, 50)
    y = np.linspace(0.1, 0.9, 50)
    z = np.linspace(0.1, 0.9, 50)

    x3d, y3d, z3d = np.meshgrid(x, y, z, indexing="ij")

    f = np.sin(2 * np.pi * x3d.ravel())

    dgrid_pre = p2g.mpi_part2grid3D(
        x3d.ravel(),
        y3d.ravel(),
        z3d.ravel(),
        1.0,
        ngrid,
        mpi,
        method="NGP",
        origin=0.0,
    )

    grid_pre = p2g.mpi_part2grid3D(
        x3d.ravel(),
        y3d.ravel(),
        z3d.ravel(),
        1.0,
        ngrid,
        mpi,
        f=f,
        method="NGP",
        origin=0.0,
    )

    fgrid = mpi_gridSPH3D(
        x3d.ravel(),
        y3d.ravel(),
        z3d.ravel(),
        1.0,
        ngrid,
        mpi,
        minpart=1,
        f=f,
        dgrid=dgrid_pre,
        fgrid=grid_pre,
    )

    assert_mpi_grid(fgrid, mpi, ngrid, ndim=3)
