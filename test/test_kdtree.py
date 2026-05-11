import numpy as np
import pytest
from fiesta.kdtree import KDTree1D, KDTree2D, KDTree3D


class TestKDTree:
    def test_kdtree1d(self):
        x = np.random.rand(100)
        tree = KDTree1D()
        tree.build_tree(x)
        nearest = tree.nearest(np.array([0.5]), k=1)
        assert len(nearest) == 1

    def test_kdtree1d_nearest_k(self):
        x = np.random.rand(100)
        tree = KDTree1D()
        tree.build_tree(x)
        nearest = tree.nearest(np.array([0.5]), k=5)
        assert len(nearest) == 5

    def test_kdtree1d_find_in_r(self):
        x = np.random.rand(100)
        tree = KDTree1D()
        tree.build_tree(x)
        points = tree.find_points_in_r(np.array([0.5]), r=0.1)
        assert isinstance(points, list)

    def test_kdtree2d(self):
        x = np.random.rand(100)
        y = np.random.rand(100)
        tree = KDTree2D()
        tree.build_tree(x, y)
        nearest = tree.nearest(np.array([0.5]), np.array([0.5]), k=1)
        assert len(nearest) == 1

    def test_kdtree2d_nearest_k(self):
        x = np.random.rand(100)
        y = np.random.rand(100)
        tree = KDTree2D()
        tree.build_tree(x, y)
        nearest = tree.nearest(np.array([0.5]), np.array([0.5]), k=5)
        assert len(nearest) == 5

    def test_kdtree2d_find_in_r(self):
        x = np.random.rand(100)
        y = np.random.rand(100)
        tree = KDTree2D()
        tree.build_tree(x, y)
        points = tree.find_points_in_r(np.array([0.5]), np.array([0.5]), r=0.1)
        assert isinstance(points, list)

    def test_kdtree3d(self):
        x = np.random.rand(100)
        y = np.random.rand(100)
        z = np.random.rand(100)
        tree = KDTree3D()
        tree.build_tree(x, y, z)
        nearest = tree.nearest(np.array([0.5]), np.array([0.5]), np.array([0.5]), k=1)
        assert len(nearest) == 1

    def test_kdtree3d_nearest_k(self):
        x = np.random.rand(100)
        y = np.random.rand(100)
        z = np.random.rand(100)
        tree = KDTree3D()
        tree.build_tree(x, y, z)
        nearest = tree.nearest(np.array([0.5]), np.array([0.5]), np.array([0.5]), k=5)
        assert len(nearest) == 5

    def test_kdtree3d_find_in_r(self):
        x = np.random.rand(100)
        y = np.random.rand(100)
        z = np.random.rand(100)
        tree = KDTree3D()
        tree.build_tree(x, y, z)
        points = tree.find_points_in_r(np.array([0.5]), np.array([0.5]), np.array([0.5]), r=0.1)
        assert isinstance(points, list)