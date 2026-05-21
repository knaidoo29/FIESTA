import numpy as np

from fiesta import p2g
from fiesta.gridsph.mpi_gridsph import mpi_gridSPH2D, mpi_gridSPH3D
from shift.mpiutils import MPI


def test_mpi_gridSPH2D():
    mpi = MPI()

    x = np.linspace(0.1, 0.9, 50)
    y = np.linspace(0.1, 0.9, 50)

    x2d, y2d = np.meshgrid(x, y, indexing="ij")

    dgrid = mpi_gridSPH2D(
        x2d.ravel(),
        y2d.ravel(),
        1.0,
        20,
        mpi,
        minpart=1,
    )

    assert isinstance(dgrid, np.ndarray)
    if mpi.rank <= 1:
        assert dgrid.shape == (7, 20)
    else:
        assert dgrid.shape == (6, 20)
    assert np.all(np.isfinite(dgrid))

    dgrid = mpi_gridSPH2D(
        x2d.ravel(),
        y2d.ravel(),
        [1.0, 1.0],
        [20, 20],
        mpi,
        minpart=1,
        periodic=[False, True],
    )

    assert isinstance(dgrid, np.ndarray)
    if mpi.rank <= 1:
        assert dgrid.shape == (7, 20)
    else:
        assert dgrid.shape == (6, 20)
    assert np.all(np.isfinite(dgrid))


def test_mpi_gridSPH2D_pregrid():
    mpi = MPI()

    x = np.linspace(0.1, 0.9, 50)
    y = np.linspace(0.1, 0.9, 50)

    x2d, y2d = np.meshgrid(x, y, indexing="ij")

    dgrid_per = p2g.mpi_part2grid2D(
        x2d.ravel(),
        y2d.ravel(),
        np.ones(len(x2d.ravel())),
        1.0,
        20,
        mpi,
        method="NGP",
        periodic=[False, True],
        origin=0.0,
    )

    dgrid = mpi_gridSPH2D(
        x2d.ravel(),
        y2d.ravel(),
        1.0,
        20,
        mpi,
        minpart=1,
        dgrid=dgrid_per,
    )

    assert isinstance(dgrid, np.ndarray)
    if mpi.rank <= 1:
        assert dgrid.shape == (7, 20)
    else:
        assert dgrid.shape == (6, 20)
    assert np.all(np.isfinite(dgrid))


def test_mpi_gridSPH2D_field():
    mpi = MPI()

    x = np.linspace(0.1, 0.9, 50)
    y = np.linspace(0.1, 0.9, 50)

    x2d, y2d = np.meshgrid(x, y, indexing="ij")

    f = np.cos(2 * np.pi * x2d.ravel())

    fgrid = mpi_gridSPH2D(
        x2d.ravel(),
        y2d.ravel(),
        1.0,
        20,
        mpi,
        minpart=1,
        f=f,
    )

    assert isinstance(fgrid, np.ndarray)
    if mpi.rank <= 1:
        assert fgrid.shape == (7, 20)
    else:
        assert fgrid.shape == (6, 20)
    assert np.all(np.isfinite(fgrid))

    fgrid = mpi_gridSPH2D(
        x2d.ravel(),
        y2d.ravel(),
        1.0,
        20,
        mpi,
        minpart=1,
        f=f,
        periodic=[False, True]
    )

    assert isinstance(fgrid, np.ndarray)
    if mpi.rank <= 1:
        assert fgrid.shape == (7, 20)
    else:
        assert fgrid.shape == (6, 20)
    assert np.all(np.isfinite(fgrid))


def test_mpi_gridSPH2D_pregrid_field():
    mpi = MPI()

    x = np.linspace(0.1, 0.9, 50)
    y = np.linspace(0.1, 0.9, 50)

    x2d, y2d = np.meshgrid(x, y, indexing="ij")

    f = np.cos(2 * np.pi * x2d.ravel())

    dgrid_pre = p2g.mpi_part2grid2D(
        x2d.ravel(),
        y2d.ravel(),
        np.ones(len(x2d.ravel())),
        1.0,
        20,
        mpi,
        method="NGP",
        periodic=[False, True],
        origin=0.0,
    )

    fgrid_pre = p2g.mpi_part2grid2D(
        x2d.ravel(),
        y2d.ravel(),
        f,
        1.0,
        20,
        mpi,
        method="NGP",
        periodic=[False, True],
        origin=0.0,
    )

    fgrid = mpi_gridSPH2D(
        x2d.ravel(),
        y2d.ravel(),
        1.0,
        20,
        mpi,
        minpart=1,
        f=f,
        dgrid=dgrid_pre,
        fgrid=fgrid_pre,
    )

    assert isinstance(fgrid, np.ndarray)
    if mpi.rank <= 1:
        assert fgrid.shape == (7, 20)
    else:
        assert fgrid.shape == (6, 20)
    assert np.all(np.isfinite(fgrid))


def test_mpi_gridSPH3D():
    mpi = MPI()

    x = np.linspace(0.1, 0.9, 50)
    y = np.linspace(0.1, 0.9, 50)
    z = np.linspace(0.1, 0.9, 50)

    x3d, y3d, z3d = np.meshgrid(x, y, z, indexing="ij")

    dgrid = mpi_gridSPH3D(
        x3d.ravel(),
        y3d.ravel(),
        z3d.ravel(),
        1.0,
        20,
        mpi,
        minpart=1,
    )

    assert isinstance(dgrid, np.ndarray)
    if mpi.rank <= 1:
        assert dgrid.shape == (7, 20, 20)
    else:
        assert dgrid.shape == (6, 20, 20)
    assert np.all(np.isfinite(dgrid))

    dgrid = mpi_gridSPH3D(
        x3d.ravel(),
        y3d.ravel(),
        z3d.ravel(),
        [1.0, 1.0, 1.0],
        [20, 20, 20],
        mpi,
        minpart=1,
        periodic=[False, True, True],
    )

    assert isinstance(dgrid, np.ndarray)
    if mpi.rank <= 1:
        assert dgrid.shape == (7, 20, 20)
    else:
        assert dgrid.shape == (6, 20, 20)
    assert np.all(np.isfinite(dgrid))


def test_mpi_pregrid_gridSPH3D():
    mpi = MPI()

    x = np.linspace(0.1, 0.9, 50)
    y = np.linspace(0.1, 0.9, 50)
    z = np.linspace(0.1, 0.9, 50)

    x3d, y3d, z3d = np.meshgrid(x, y, z, indexing="ij")

    dgrid_per = p2g.mpi_part2grid3D(
        x3d.ravel(),
        y3d.ravel(),
        z3d.ravel(),
        np.ones(len(x3d.ravel())),
        1.0,
        20,
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
        20,
        mpi,
        minpart=1,
        dgrid=dgrid_per,
    )

    assert isinstance(dgrid, np.ndarray)
    if mpi.rank <= 1:
        assert dgrid.shape == (7, 20, 20)
    else:
        assert dgrid.shape == (6, 20, 20)
    assert np.all(np.isfinite(dgrid))


def test_mpi_gridSPH3D_field():
    mpi = MPI()

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
        20,
        mpi,
        minpart=1,
        f=f,
    )

    assert isinstance(fgrid, np.ndarray)
    if mpi.rank <= 1:
        assert fgrid.shape == (7, 20, 20)
    else:
        assert fgrid.shape == (6, 20, 20)
    assert np.all(np.isfinite(fgrid))

    fgrid = mpi_gridSPH3D(
        x3d.ravel(),
        y3d.ravel(),
        z3d.ravel(),
        1.0,
        20,
        mpi,
        minpart=1,
        f=f,
        periodic=[False, True, True]
    )

    assert isinstance(fgrid, np.ndarray)
    if mpi.rank <= 1:
        assert fgrid.shape == (7, 20, 20)
    else:
        assert fgrid.shape == (6, 20, 20)
    assert np.all(np.isfinite(fgrid))



def test_mpi_gridSPH3D_pregrid_field():
    mpi = MPI()

    x = np.linspace(0.1, 0.9, 50)
    y = np.linspace(0.1, 0.9, 50)
    z = np.linspace(0.1, 0.9, 50)

    x3d, y3d, z3d = np.meshgrid(x, y, z, indexing="ij")

    f = np.sin(2 * np.pi * x3d.ravel())

    dgrid_pre = p2g.mpi_part2grid3D(
        x3d.ravel(),
        y3d.ravel(),
        z3d.ravel(),
        np.ones(len(x3d.ravel())),
        1.0,
        20,
        mpi,
        method="NGP",
        origin=0.0,
    )

    grid_pre = p2g.mpi_part2grid3D(
        x3d.ravel(),
        y3d.ravel(),
        z3d.ravel(),
        f,
        1.0,
        20,
        mpi,
        method="NGP",
        origin=0.0,
    )

    fgrid = mpi_gridSPH3D(
        x3d.ravel(),
        y3d.ravel(),
        z3d.ravel(),
        1.0,
        20,
        mpi,
        minpart=1,
        f=f,
        dgrid=dgrid_pre,
        fgrid=grid_pre,
    )

    assert isinstance(fgrid, np.ndarray)
    if mpi.rank <= 1:
        assert fgrid.shape == (7, 20, 20)
    else:
        assert fgrid.shape == (6, 20, 20)
    assert np.all(np.isfinite(fgrid))
