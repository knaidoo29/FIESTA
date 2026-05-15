import numpy as np

from fiesta.gridsph.mpi_gridsph import mpi_gridSPH2D, mpi_gridSPH3D
from shift.mpiutils import MPI


def test_mpi_gridSPH2D():
    mpi = MPI()

    x = np.linspace(0.1, 0.9, 50)
    y = np.linspace(0.1, 0.9, 50)

    x2d, y2d = np.meshgrid(x, y, indexing="ij")

    dgrid = mpi_gridSPH2D(
        1.0,
        20,
        x2d.ravel(),
        y2d.ravel(),
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
        [1.0, 1.0],
        [20,20],
        x2d.ravel(),
        y2d.ravel(),
        mpi,
        minpart=1,
        periodic=[True, True],
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
        1.0,
        20,
        x2d.ravel(),
        y2d.ravel(),
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


def test_mpi_gridSPH3D():
    mpi = MPI()

    x = np.linspace(0.1, 0.9, 50)
    y = np.linspace(0.1, 0.9, 50)
    z = np.linspace(0.1, 0.9, 50)

    x3d, y3d, z3d = np.meshgrid(x, y, z, indexing="ij")

    dgrid = mpi_gridSPH3D(
        1.0,
        20,
        x3d.ravel(),
        y3d.ravel(),
        z3d.ravel(),
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
        [1.0, 1.0, 1.0],
        [20, 20, 20],
        x3d.ravel(),
        y3d.ravel(),
        z3d.ravel(),
        mpi,
        minpart=1,
        periodic=[True, True, True],
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
        1.0,
        20,
        x3d.ravel(),
        y3d.ravel(),
        z3d.ravel(),
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