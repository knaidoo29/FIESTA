import numpy as np
import pytest
from fiesta.coords import (
    x2points,
    points2x,
    xy2points,
    points2xy,
    xyz2points,
    points2xyz,
    coord2points,
)


class TestCoords:
    def test_x2points(self):
        x = np.array([1, 2, 3])
        points = x2points(x)
        assert points.shape == (3, 1)
        assert np.array_equal(points.flatten(), x)

    def test_x2points_scalar(self):
        x = 1.0
        points = x2points(x)
        assert points.shape == (1, 1)

    def test_points2x(self):
        points = np.array([[1], [2], [3]])
        x = points2x(points)
        assert len(x) == 3
        assert np.array_equal(x, [1, 2, 3])

    def test_xy2points(self):
        x = np.array([1, 2, 3])
        y = np.array([4, 5, 6])
        points = xy2points(x, y)
        assert points.shape == (3, 2)
        assert np.array_equal(points[:, 0], x)
        assert np.array_equal(points[:, 1], y)

    def test_xy2points_scalar(self):
        x = 1.0
        y = 2.0
        points = xy2points(x, y)
        assert points.shape == (1, 2)
        assert np.array_equal(points[:, 0], [x])
        assert np.array_equal(points[:, 1], [y])

    def test_points2xy(self):
        points = np.array([[1, 4], [2, 5], [3, 6]])
        x, y = points2xy(points)
        assert np.array_equal(x, [1, 2, 3])
        assert np.array_equal(y, [4, 5, 6])

    def test_xyz2points(self):
        x = np.array([1, 2, 3])
        y = np.array([4, 5, 6])
        z = np.array([7, 8, 9])
        points = xyz2points(x, y, z)
        assert points.shape == (3, 3)
        assert np.array_equal(points[:, 0], x)
        assert np.array_equal(points[:, 1], y)
        assert np.array_equal(points[:, 2], z)

    def test_xyz2points_scalar(self):
        x = 1.0
        y = 2.0
        z = 3.0
        points = xyz2points(x, y, z)
        assert points.shape == (1, 3)
        assert np.array_equal(points[:, 0], [x])
        assert np.array_equal(points[:, 1], [y])
        assert np.array_equal(points[:, 2], [z])

    def test_points2xyz(self):
        points = np.array([[1, 4, 7], [2, 5, 8], [3, 6, 9]])
        x, y, z = points2xyz(points)
        assert np.array_equal(x, [1, 2, 3])
        assert np.array_equal(y, [4, 5, 6])
        assert np.array_equal(z, [7, 8, 9])

    def test_coord2points(self):
        x = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]])
        points = coord2points(x)
        assert points.shape == (3, 4)
        assert np.array_equal(points.flatten(), x.T.flatten())

    def test_coord2points_scalar(self):
        x = np.array([1.0, 2.0, 3.0, 4.0])
        points = coord2points(x)
        assert points.shape == (1, 4)
        assert np.array_equal(points.flatten(), x.T.flatten())
