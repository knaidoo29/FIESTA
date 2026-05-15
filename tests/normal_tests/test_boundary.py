import numpy as np
import pytest
from fiesta.boundary import (
    buffer_periodic,
    subbox_buffer_periodic,
    buffer_periodic_2D,
    subbox_buffer_periodic_2D,
    buffer_periodic_3D,
    subbox_buffer_periodic_3D,
    buffer_random_2D,
    buffer_random_3D,
)


class TestBoundary:
    def test_buffer_periodic(self):
        data = np.random.rand(10, 2)
        result = buffer_periodic(data, axis=0, boxsize=1.0, buffer_length=0.1)
        assert result.shape[1] == data.shape[1]
        assert result.shape[0] >= data.shape[0]

    def test_subbox_buffer_periodic(self):
        data = np.random.rand(10, 2)
        result = subbox_buffer_periodic(
            data, axis=0, boxsize=1.0, buffer_length=0.1, subboxsize=0.5
        )
        assert result.shape[1] == data.shape[1]
        # May not add particles if subbox is small

    def test_buffer_periodic_2D(self):
        data = np.random.rand(10, 2)
        result = buffer_periodic_2D(data, boxsize=1.0, buffer_length=0.1)
        assert result.shape[1] == 2
        assert result.shape[0] >= data.shape[0]

    def test_buffer_periodic_2D_list_params(self):
        data = np.random.rand(10, 2)
        result = buffer_periodic_2D(
            data, boxsize=[1.0, 2.0], buffer_length=0.1, origin=[0.0, 1.0]
        )
        assert result.shape[1] == 2
        assert result.shape[0] >= data.shape[0]

    def test_subbox_buffer_periodic_2D(self):
        data = np.random.rand(10, 2)
        result = subbox_buffer_periodic_2D(
            data, boxsize=1.0, buffer_length=0.1, subboxsize=0.5
        )
        assert result.shape[1] == 2
        # May not add particles

    def test_subbox_buffer_periodic_2D_list_params(self):
        data = np.random.rand(10, 2)
        result = subbox_buffer_periodic_2D(
            data,
            boxsize=[1.0, 2.0],
            buffer_length=0.1,
            subboxsize=[0.5, 1.0],
            origin=[0.0, 1.0],
            subbox_origin=[0.0, 1.0],
        )
        assert result.shape[1] == 2
        # May not add particles

    def test_buffer_periodic_3D(self):
        data = np.random.rand(10, 3)
        result = buffer_periodic_3D(data, boxsize=1.0, buffer_length=0.1)
        assert result.shape[1] == 3
        assert result.shape[0] >= data.shape[0]

    def test_buffer_periodic_3D_list_params(self):
        data = np.random.rand(10, 3)
        result = buffer_periodic_3D(
            data, boxsize=[1.0, 2.0, 3.0], buffer_length=0.1, origin=[0.0, 1.0, 2.0]
        )
        assert result.shape[1] == 3
        assert result.shape[0] >= data.shape[0]

    def test_subbox_buffer_periodic_3D(self):
        data = np.random.rand(10, 3)
        result = subbox_buffer_periodic_3D(
            data, boxsize=1.0, buffer_length=0.1, subboxsize=0.5
        )
        assert result.shape[1] == 3
        # May not add particles

    def test_subbox_buffer_periodic_3D_list_params(self):
        data = np.random.rand(10, 3)
        result = subbox_buffer_periodic_3D(
            data,
            boxsize=[1.0, 2.0, 3.0],
            buffer_length=0.1,
            subboxsize=[0.5, 1.0, 1.5],
            origin=[0.0, 1.0, 2.0],
            subbox_origin=[0.0, 1.0, 2.0],
        )
        assert result.shape[1] == 3
        # May not add particles

    def test_buffer_random_2D(self):
        x, y = buffer_random_2D(npart=5, boxsize=1.0, buffer_length=0.1)
        # May return empty if buffer_length is small
        assert isinstance(x, np.ndarray)
        assert isinstance(y, np.ndarray)

    def test_buffer_random_2D_list_params(self):
        x, y = buffer_random_2D(
            npart=5, boxsize=[1.0, 2.0], buffer_length=0.1, origin=[0.0, 1.0]
        )
        # May return empty if buffer_length is small
        assert isinstance(x, np.ndarray)
        assert isinstance(y, np.ndarray)

    def test_buffer_random_3D(self):
        x, y, z = buffer_random_3D(npart=5, boxsize=1.0, buffer_length=0.1)
        assert isinstance(x, np.ndarray)
        assert isinstance(y, np.ndarray)
        assert isinstance(z, np.ndarray)

    def test_buffer_random_3D_list_params(self):
        x, y, z = buffer_random_3D(
            npart=5, boxsize=[1.0, 2.0, 3.0], buffer_length=0.1, origin=[0.0, 1.0, 2.0]
        )
        assert isinstance(x, np.ndarray)
        assert isinstance(y, np.ndarray)
        assert isinstance(z, np.ndarray)
