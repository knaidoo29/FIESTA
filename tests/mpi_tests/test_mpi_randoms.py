import numpy as np
import pytest

pytest.importorskip("mpi4py")
from shift.mpiutils import MPI

from fiesta.boundary.mpi_randoms import (
    mpi_buffer_random_2D,
    mpi_buffer_random_3D,
    mpi_buffer_random_utils,
)


def test_mpi_buffer_random_2D():
    mpi = MPI()
    np.random.seed(123)

    npart = 20
    boxsize = 2.0
    limits = [0.0, 1.0, 0.0, 1.0]
    buffer_length = 0.2

    x, y = mpi_buffer_random_2D(npart, boxsize, limits, buffer_length, mpi)

    assert isinstance(x, np.ndarray)
    assert isinstance(y, np.ndarray)
    assert x.shape == y.shape
    assert x.ndim == 1
    assert x.size >= 0
    assert np.all(np.isfinite(x))
    assert np.all(np.isfinite(y))


def test_mpi_buffer_random_2D_with_list_boxsize():
    mpi = MPI()
    np.random.seed(456)

    npart = 15
    boxsize = [2.0, 2.0]
    limits = [0.0, 1.0, 0.0, 1.0]
    buffer_length = 0.25

    x, y = mpi_buffer_random_2D(npart, boxsize, limits, buffer_length, mpi)

    assert isinstance(x, np.ndarray)
    assert x.shape == y.shape
    assert np.all(np.isfinite(x))
    assert np.all(np.isfinite(y))


def test_mpi_buffer_random_3D():
    mpi = MPI()
    np.random.seed(789)

    npart = 18
    boxsize = 2.0
    limits = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
    buffer_length = 0.2

    x, y, z = mpi_buffer_random_3D(npart, boxsize, limits, buffer_length, mpi)

    assert isinstance(x, np.ndarray)
    assert x.shape == y.shape == z.shape
    assert x.ndim == 1
    assert x.size >= 0
    assert np.all(np.isfinite(x))
    assert np.all(np.isfinite(y))
    assert np.all(np.isfinite(z))


def test_mpi_buffer_random_3D_with_list_boxsize():
    mpi = MPI()
    np.random.seed(987)

    npart = 12
    boxsize = [2.0, 2.0, 2.0]
    limits = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
    buffer_length = 0.15

    x, y, z = mpi_buffer_random_3D(npart, boxsize, limits, buffer_length, mpi)

    assert isinstance(x, np.ndarray)
    assert x.shape == y.shape == z.shape
    assert np.all(np.isfinite(x))
    assert np.all(np.isfinite(y))
    assert np.all(np.isfinite(z))


def test_mpi_buffer_random_utils():
    mpi = MPI()
    np.random.seed(321)

    data = np.column_stack(
        [
            np.linspace(0.1, 0.9, 10),
            np.linspace(0.2, 0.8, 10),
        ]
    )
    limits = [0.0, 1.0, 0.0, 1.0]
    buffer_length = 0.1

    result = mpi_buffer_random_utils(data, limits, buffer_length, mpi)

    assert isinstance(result, np.ndarray)
    assert result.ndim == 2
    assert result.shape[1] == 2
    assert result.shape[0] >= data.shape[0]
    assert np.all(np.isfinite(result))
