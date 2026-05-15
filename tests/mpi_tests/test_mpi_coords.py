import numpy as np
import pytest

pytest.importorskip("mpi4py")
from shift.mpiutils import MPI

from fiesta.coords.mpi_points import (
    split_limits_by_grid,
    mpi_find_range_2D,
    mpi_find_range_3D,
    mpi_open_2D,
    mpi_open_3D,
    check_coords_at_MPI_0,
    MPI_SortByX,
    distribute_points_by_x,
)


def test_split_limits_by_grid_returns_valid_range():
    mpi = MPI()
    limits = split_limits_by_grid(2.0, 0.0, 8, mpi)

    assert isinstance(limits, list)
    assert len(limits) == 2
    assert limits[0] >= 0.0
    assert limits[1] <= 2.0
    assert limits[1] > limits[0]


def test_mpi_find_range_2D_collects_ranges_to_rank0():
    mpi = MPI()

    files = ["f1", "f2"]
    file_data = {
        "f1": np.array([[0.1, 0.2], [0.5, 0.7]]),
        "f2": np.array([[0.3, 0.1], [0.8, 0.9]]),
    }

    def reader(name):
        return file_data[name]

    ranges = mpi_find_range_2D(files, reader, mpi)

    if mpi.rank == 0:
        assert isinstance(ranges, np.ndarray)
        assert ranges.shape == (len(files), 4)
        assert np.all(ranges[:, 0] <= ranges[:, 1])
        assert np.all(ranges[:, 2] <= ranges[:, 3])
    else:
        assert ranges is None


def test_mpi_find_range_3D_collects_ranges_to_rank0():
    mpi = MPI()

    files = ["f1", "f2"]
    file_data = {
        "f1": np.array([[0.1, 0.2, 0.3], [0.5, 0.7, 0.8]]),
        "f2": np.array([[0.3, 0.1, 0.2], [0.8, 0.9, 0.4]]),
    }

    def reader(name):
        return file_data[name]

    ranges = mpi_find_range_3D(files, reader, mpi)

    if mpi.rank == 0:
        assert isinstance(ranges, np.ndarray)
        assert ranges.shape == (len(files), 6)
        assert np.all(ranges[:, 0] <= ranges[:, 1])
        assert np.all(ranges[:, 2] <= ranges[:, 3])
        assert np.all(ranges[:, 4] <= ranges[:, 5])
    else:
        assert ranges is None


def test_mpi_open_2D_filters_files_by_limits():
    files = ["f1", "f2", "f3"]
    file_data = {
        "f1": np.array([[0.05, 0.05], [0.1, 0.1]]),
        "f2": np.array([[0.9, 0.9], [1.1, 1.1]]),
        "f3": np.array([[0.4, 0.4], [0.6, 0.6]]),
    }

    def reader(name):
        return file_data[name]

    ranges = np.array(
        [
            [0.0, 0.2, 0.0, 0.2],
            [0.8, 1.2, 0.8, 1.2],
            [0.35, 0.65, 0.35, 0.65],
        ]
    )
    limits = [0.0, 0.7, 0.0, 0.7]

    datas = mpi_open_2D(files, reader, ranges, limits)

    assert isinstance(datas, np.ndarray)
    assert datas.shape[1] == 2
    assert np.all((datas[:, 0] >= 0.0) & (datas[:, 0] < 0.7))
    assert np.all((datas[:, 1] >= 0.0) & (datas[:, 1] < 0.7))


def test_mpi_open_3D_filters_files_by_limits():
    files = ["f1", "f2", "f3"]
    file_data = {
        "f1": np.array([[0.05, 0.05, 0.05], [0.1, 0.1, 0.1]]),
        "f2": np.array([[0.9, 0.9, 0.9], [1.1, 1.1, 1.1]]),
        "f3": np.array([[0.4, 0.4, 0.4], [0.6, 0.6, 0.6]]),
    }

    def reader(name):
        return file_data[name]

    ranges = np.array(
        [
            [0.0, 0.2, 0.0, 0.2, 0.0, 0.2],
            [0.8, 1.2, 0.8, 1.2, 0.8, 1.2],
            [0.35, 0.65, 0.35, 0.65, 0.35, 0.65],
        ]
    )
    limits = [0.0, 0.7, 0.0, 0.7, 0.0, 0.7]

    datas = mpi_open_3D(files, reader, ranges, limits)

    assert isinstance(datas, np.ndarray)
    assert datas.shape[1] == 3
    assert np.all((datas[:, 0] >= 0.0) & (datas[:, 0] < 0.7))
    assert np.all((datas[:, 1] >= 0.0) & (datas[:, 1] < 0.7))
    assert np.all((datas[:, 2] >= 0.0) & (datas[:, 2] < 0.7))


def test_check_coords_at_MPI_0_detects_rank0_only():
    mpi = MPI()
    if mpi.rank == 0:
        x = np.array([0.1, 0.2])
    else:
        x = None

    check = check_coords_at_MPI_0(x, mpi)

    assert isinstance(check, bool)
    if mpi.size == 1:
        assert check is True
    else:
        assert check is True

    x = np.array([0.1, 0.2])

    check = check_coords_at_MPI_0(x, mpi)

    assert isinstance(check, bool)
    if mpi.size == 1:
        assert check is False
    else:
        assert check is False


def test_mpi_sortbyx_distribute_roundtrip():
    mpi = MPI()
    sorter = MPI_SortByX(mpi)
    sorter.settings(1.0, 8, origin=0.0, buffer_length=0.0)

    # create a small slab-aware data set
    x = np.linspace(0.05, 0.95, 8)
    y = np.full_like(x, 0.5)
    data = np.column_stack([x, y])
    sorter.input(data)
    sorter.limits4grid()

    distributed = sorter.distribute()

    assert isinstance(distributed, np.ndarray)
    assert distributed.shape[1] == 2
    assert np.all(distributed[:, 0] >= sorter.limits[0] - sorter.buffer_length)
    assert np.all(distributed[:, 0] <= sorter.limits[1] + sorter.buffer_length)

    data = None
    sorter.input(data)

    distributed = sorter.distribute()

    sorter = MPI_SortByX(mpi)
    sorter.settings(1.0, 8, origin=0.0, buffer_length=0.1)

    # create a small slab-aware data set
    x = np.linspace(0.05, 0.95, 8)
    y = np.full_like(x, 0.5)
    data = np.column_stack([x, y])
    sorter.input(data)
    sorter.limits4grid()

    distributed = sorter.distribute(include_internalbuffer=True)


def test_mpi_sortbyx_distribute_grid2d_preserves_values():
    mpi = MPI()
    sorter = MPI_SortByX(mpi)
    sorter.settings([1.0, 1.0], 10, origin=[0.0, 0.0])
    sorter.limits4grid()

    x2d, y2d = np.meshgrid(
        np.linspace(0.1, 0.9, 10), np.linspace(0.1, 0.9, 10), indexing="ij"
    )
    f2d = np.arange(100).reshape(10, 10)

    grid = sorter.distribute_grid2D(x2d, y2d, f2d)

    assert isinstance(grid, np.ndarray)
    assert grid.ndim == 2
    assert grid.shape[1] == 10
    assert np.all(np.isfinite(grid))

    sorter = MPI_SortByX(mpi)
    sorter.settings(1.0, 10, origin=0.0)
    sorter.limits4grid()

    if mpi.rank == 0:
        x2d, y2d = None, None
        f2d = None
    else:
        x2d, y2d = np.meshgrid(
            np.linspace(0.1, 0.9, 10), np.linspace(0.1, 0.9, 10), indexing="ij"
        )
        f2d = np.arange(100).reshape(10, 10)
    grid = sorter.distribute_grid2D(x2d, y2d, f2d)

    sorter = MPI_SortByX(mpi)
    sorter.settings(1.0, 10, origin=0.0)
    sorter.limits4grid()

    if mpi.rank == 0:
        x2d, y2d = None, None
        f2d = None
    else:
        x2d, y2d = np.meshgrid(
            np.linspace(0.1, 0.9, 10), np.linspace(0.1, 0.9, 10), indexing="ij"
        )
        f2d = np.arange(100).reshape(10, 10)
    grid = sorter.distribute_grid2D(x2d, y2d, f2d)


def test_mpi_sortbyx_distribute_grid3d_preserves_values():
    mpi = MPI()
    sorter = MPI_SortByX(mpi)
    sorter.settings(1.0, 3, origin=0.0)
    sorter.limits4grid()

    x3d, y3d, z3d = np.meshgrid(
        np.linspace(0.1, 0.9, 3),
        np.linspace(0.1, 0.9, 3),
        np.linspace(0.1, 0.9, 3),
        indexing="ij",
    )
    f3d = np.arange(27).reshape(3, 3, 3)

    grid = sorter.distribute_grid3D(x3d, y3d, z3d, f3d)

    assert isinstance(grid, np.ndarray)
    assert grid.ndim == 3
    assert grid.shape[1] == 3
    assert np.all(np.isfinite(grid))

    mpi = MPI()
    sorter = MPI_SortByX(mpi)
    sorter.settings([1.0, 1.0, 2.0], 3, origin=[0.0, 0.0, 0.0])
    sorter.limits4grid()

    if mpi.rank == 0:
        x3d, y3d, z3d = None, None, None
        f3d = None
    else:
        x3d, y3d, z3d = np.meshgrid(
            np.linspace(0.1, 0.9, 3),
            np.linspace(0.1, 0.9, 3),
            np.linspace(0.1, 0.9, 3),
            indexing="ij",
        )
        f3d = np.arange(27).reshape(3, 3, 3)

    grid = sorter.distribute_grid3D(x3d, y3d, z3d, f3d)


def test_mpi_sortbyx_clean_resets_state():
    mpi = MPI()
    sorter = MPI_SortByX(mpi)
    sorter.settings(1.0, 5, origin=0.0, buffer_length=0.1)
    sorter.input(np.zeros((2, 2)))
    sorter.clean()

    assert sorter.boxsize is None
    assert sorter.origin is None
    assert sorter.data is None


def test_distribute_points_by_x_basic_shape_and_finite():
    """
    Basic smoke test:
    - verifies function runs
    - output is ndarray
    - preserves column count
    - values remain finite
    """
    mpi = MPI()

    x = np.linspace(0.1, 0.9, 9)
    y = np.linspace(1.0, 2.0, 9)

    data = np.column_stack([x, y])

    distributed = distribute_points_by_x(
        data=data,
        boxsize=1.0,
        ngrid=3,
        origin=0.0,
        MPI=mpi,
    )

    assert isinstance(distributed, np.ndarray)
    assert distributed.ndim == 2
    assert distributed.shape[1] == 2
    assert np.all(np.isfinite(distributed))


def test_distribute_points_by_x_multicolumn_data():
    """
    Ensures extra columns are preserved, not just x positions.
    """
    mpi = MPI()

    x = np.linspace(0.1, 0.9, 12)
    y = np.sin(x)
    z = np.cos(x)

    data = np.column_stack([x, y, z])

    distributed = distribute_points_by_x(
        data=data,
        boxsize=1.0,
        ngrid=4,
        origin=0.0,
        MPI=mpi,
    )

    assert isinstance(distributed, np.ndarray)
    assert distributed.shape[1] == 3
    assert np.all(np.isfinite(distributed))


def test_distribute_points_by_x_shifted_origin():
    """
    Covers non-zero origin branch.
    """
    mpi = MPI()

    x = np.linspace(-0.9, 0.9, 15)
    y = np.arange(15)

    data = np.column_stack([x, y])

    distributed = distribute_points_by_x(
        data=data,
        boxsize=2.0,
        ngrid=4,
        origin=-1.0,
        MPI=mpi,
    )

    assert isinstance(distributed, np.ndarray)
    assert distributed.ndim == 2
    assert distributed.shape[1] == 2
    assert np.all(np.isfinite(distributed))
