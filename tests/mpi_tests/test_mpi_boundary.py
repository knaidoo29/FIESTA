import numpy as np
import pytest

pytest.importorskip("mpi4py")
import shift
from shift.mpiutils import MPI

from fiesta.boundary.mpi_periodic import (
    mpi_buffer_periodic_2D,
    mpi_buffer_periodic_3D,
    mpi_buffer_internal_3D,
)


def test_mpi_buffer_periodic_2D():
    """Test 2D periodic buffer generation."""
    mpi = MPI()
    rng = np.random.default_rng(42)

    boxsize = 1.0
    ngrid = 8
    buffer_length = 0.1
    origin = 0.0

    x = rng.random(20)
    y = rng.random(20)
    data = np.column_stack([x, y])

    result = mpi_buffer_periodic_2D(data, boxsize, buffer_length, mpi, origin=origin)

    # Check that result is a numpy array
    assert isinstance(result, np.ndarray)
    # Original data should be included
    assert result.shape[0] >= data.shape[0]
    # Check columns
    assert result.shape[1] == 2
    # Check values are finite
    assert np.all(np.isfinite(result))


def test_mpi_buffer_periodic_2D_with_list_boxsize():
    """Test 2D periodic buffer with list boxsize."""
    mpi = MPI()
    rng = np.random.default_rng(123)

    boxsize = [1.0, 1.0]
    buffer_length = 0.1
    origin = [0.0, 0.0]

    x = rng.random(15)
    y = rng.random(15)
    data = np.column_stack([x, y])

    result = mpi_buffer_periodic_2D(data, boxsize, buffer_length, mpi, origin=origin)

    assert isinstance(result, np.ndarray)
    assert result.shape[1] == 2
    assert np.all(np.isfinite(result))


def test_mpi_buffer_periodic_3D():
    """Test 3D periodic buffer generation."""
    mpi = MPI()
    rng = np.random.default_rng(234)

    boxsize = 1.0
    buffer_length = 0.1
    origin = 0.0

    x = rng.random(20)
    y = rng.random(20)
    z = rng.random(20)
    data = np.column_stack([x, y, z])

    result = mpi_buffer_periodic_3D(data, boxsize, buffer_length, mpi, origin=origin)

    # Check that result is a numpy array
    assert isinstance(result, np.ndarray)
    # Original data should be included
    assert result.shape[0] >= data.shape[0]
    # Check columns (x, y, z)
    assert result.shape[1] == 3
    # Check values are finite
    assert np.all(np.isfinite(result))


def test_mpi_buffer_periodic_3D_with_list_boxsize():
    """Test 3D periodic buffer with list boxsize."""
    mpi = MPI()
    rng = np.random.default_rng(345)

    boxsize = [1.0, 1.0, 1.0]
    buffer_length = 0.15
    origin = [0.0, 0.0, 0.0]

    x = rng.random(20)
    y = rng.random(20)
    z = rng.random(20)
    data = np.column_stack([x, y, z])

    result = mpi_buffer_periodic_3D(data, boxsize, buffer_length, mpi, origin=origin)

    assert isinstance(result, np.ndarray)
    assert result.shape[1] == 3
    assert np.all(np.isfinite(result))


def test_mpi_buffer_internal_3D():
    """Test 3D internal buffer particle exchange."""
    mpi = MPI()
    rng = np.random.default_rng(456)

    boxsize = 1.0
    buffer_length = 0.1
    origin = 0.0

    x = rng.random(20)
    y = rng.random(20)
    z = rng.random(20)
    data = np.column_stack([x, y, z])

    result = mpi_buffer_internal_3D(data, boxsize, buffer_length, mpi, origin=origin)

    # Check that result is a numpy array
    assert isinstance(result, np.ndarray)
    # Original data should be included
    assert result.shape[0] >= data.shape[0]
    # Check columns
    assert result.shape[1] == 3
    # Check values are finite
    assert np.all(np.isfinite(result))


def test_mpi_buffer_internal_3D_with_list_boxsize():
    """Test 3D internal buffer with list boxsize."""
    mpi = MPI()
    rng = np.random.default_rng(567)

    boxsize = [1.0, 1.0, 1.0]
    buffer_length = 0.12
    origin = [0.0, 0.0, 0.0]

    x = rng.random(20)
    y = rng.random(20)
    z = rng.random(20)
    data = np.column_stack([x, y, z])

    result = mpi_buffer_internal_3D(data, boxsize, buffer_length, mpi, origin=origin)

    assert isinstance(result, np.ndarray)
    assert result.shape[1] == 3
    assert np.all(np.isfinite(result))


def test_mpi_buffer_periodic_2D_with_nonzero_origin():
    """Test 2D buffer with non-zero origin."""
    mpi = MPI()
    rng = np.random.default_rng(678)

    boxsize = 10.0
    buffer_length = 0.5
    origin = 5.0

    x = rng.uniform(origin, origin + boxsize, 15)
    y = rng.uniform(origin, origin + boxsize, 15)
    data = np.column_stack([x, y])

    result = mpi_buffer_periodic_2D(data, boxsize, buffer_length, mpi, origin=origin)

    assert isinstance(result, np.ndarray)
    assert result.shape[1] == 2
    assert np.all(np.isfinite(result))


def test_mpi_buffer_periodic_3D_boundary_particles():
    """Test that 3D buffer correctly adds particles near boundaries."""
    mpi = MPI()

    boxsize = 1.0
    buffer_length = 0.2
    origin = 0.0

    # Create particles near y and z boundaries
    x = np.array([0.5, 0.5, 0.5, 0.5])
    y = np.array([0.02, 0.98, 0.5, 0.5])
    z = np.array([0.02, 0.02, 0.98, 0.5])
    data = np.column_stack([x, y, z])

    result = mpi_buffer_periodic_3D(data, boxsize, buffer_length, mpi, origin=origin)

    # Should have original particles plus periodic copies
    assert result.shape[0] > data.shape[0]
    assert result.shape[1] == 3
    assert np.all(np.isfinite(result))
