import numpy as np
import pytest
from fiesta.gridsph import (
    get_integral_image2D,
    get_integral_image3D,
    gridSPH2D,
    gridSPH3D,
)

from fiesta.src import (
    sum_from_integral_image_2D,
    get_volume_enclosing_box_2D,
    sum_from_integral_image_3D,
    get_volume_enclosing_box_3D,
)


class TestVEFE:
    def test_get_integral_image2D(self):
        fgrid = np.random.rand(5, 5)
        iimg = get_integral_image2D(fgrid)
        assert iimg.shape == (6, 6)  # integral image is padded

    def test_get_integral_image3D(self):
        fgrid = np.random.rand(3, 3, 3)
        iimg = get_integral_image3D(fgrid)
        assert iimg.shape == (4, 4, 4)

    def test_sum_from_integral_image_2D(self):
        iimg = np.random.rand(6, 6)
        sum_val = sum_from_integral_image_2D(iimg, 1, 3, 1, 3)
        assert isinstance(sum_val, (float, np.ndarray))

    def test_sum_from_integral_image_2D_periodic_wrap(self):
        fgrid = np.arange(9, dtype=np.float64).reshape(3, 3)
        iimg = get_integral_image2D(fgrid)
        wrapped = sum_from_integral_image_2D(
            iimg, -1, 1, 1, 2, xperiodic=True, yperiodic=True
        )
        assert np.isfinite(wrapped)

    def test_sum_from_integral_image_2D_nonperiodic_clamp(self):
        fgrid = np.arange(16, dtype=np.float64).reshape(4, 4)
        iimg = get_integral_image2D(fgrid)
        clipped = sum_from_integral_image_2D(
            iimg, -2, 2, -1, 1, xperiodic=False, yperiodic=False
        )
        expected = sum_from_integral_image_2D(
            iimg, 0, 2, 0, 1, xperiodic=False, yperiodic=False
        )
        assert np.isclose(clipped, expected)

    def test_sum_from_integral_image_2D_py_func(self):
        fgrid = np.arange(9, dtype=np.float64).reshape(3, 3)
        iimg = get_integral_image2D(fgrid)
        wrapped = sum_from_integral_image_2D(
            iimg, -1, 1, 1, 2, xperiodic=True, yperiodic=True
        )
        py_result = sum_from_integral_image_2D.py_func(
            iimg, -1, 1, 1, 2, xperiodic=True, yperiodic=True
        )
        assert np.isclose(wrapped, py_result)

    def test_sum_from_integral_image_2D_py_func_direct(self):
        fgrid = np.arange(16, dtype=np.float64).reshape(4, 4)
        iimg = get_integral_image2D(fgrid)
        result = sum_from_integral_image_2D.py_func(
            iimg, 1, 3, 1, 3, xperiodic=True, yperiodic=True
        )
        assert result == 30.0

    def test_get_volume_enclosing_box_2D(self):
        dgrid = np.random.rand(5, 5)
        idgrid = get_integral_image2D(dgrid)
        volume = get_volume_enclosing_box_2D(1.0, 1.0, dgrid, idgrid, minpart=1)
        assert volume.shape == (5, 5)

    def test_get_volume_enclosing_box_2D_field(self):
        dgrid = np.full((4, 4), 0.5, dtype=np.float64)
        idgrid = get_integral_image2D(dgrid)
        fgrid = np.arange(16, dtype=np.float64).reshape(4, 4)
        ifgrid = get_integral_image2D(fgrid)
        fvolume = get_volume_enclosing_box_2D(
            4.0,
            4.0,
            dgrid,
            idgrid,
            minpart=1,
            xperiodic=False,
            yperiodic=False,
            ifgrid=ifgrid,
        )
        assert fvolume.shape == (4, 4)
        assert np.all(np.isfinite(fvolume))

    def test_get_volume_enclosing_box_2D_py_func_field(self):
        dgrid = np.array(
            [[10.0, 0.2, 0.2], [0.2, 10.0, 0.2], [0.2, 0.2, 10.0]], dtype=np.float64
        )
        idgrid = get_integral_image2D(dgrid)
        fgrid = np.arange(9, dtype=np.float64).reshape(3, 3)
        ifgrid = get_integral_image2D(fgrid)
        fvolume_py = get_volume_enclosing_box_2D.py_func(
            3.0,
            3.0,
            dgrid,
            idgrid,
            minpart=1,
            xperiodic=False,
            yperiodic=False,
            ifgrid=ifgrid,
        )
        assert fvolume_py.shape == (3, 3)
        assert np.all(np.isfinite(fvolume_py))

    def test_get_volume_enclosing_box_2D_small_density(self):
        dgrid = np.full((3, 3), 0.2, dtype=np.float64)
        idgrid = get_integral_image2D(dgrid)
        volume = get_volume_enclosing_box_2D(
            3.0, 3.0, dgrid, idgrid, minpart=1, xperiodic=True, yperiodic=True
        )
        assert volume.shape == (3, 3)
        assert np.all(np.isfinite(volume))

    def test_sum_from_integral_image_3D(self):
        iimg = np.random.rand(4, 4, 4)
        sum_val = sum_from_integral_image_3D(iimg, 1, 2, 1, 2, 1, 2)
        assert isinstance(sum_val, (float, np.ndarray))

    def test_sum_from_integral_image_3D_periodic_wrap(self):
        fgrid = np.arange(27, dtype=np.float64).reshape(3, 3, 3)
        iimg = get_integral_image3D(fgrid)
        wrapped = sum_from_integral_image_3D(
            iimg, -2, 2, 1, 2, 1, 2, xperiodic=True, yperiodic=True, zperiodic=True
        )
        clipped = sum_from_integral_image_3D(
            iimg, -2, 2, 1, 2, 1, 2, xperiodic=False, yperiodic=False, zperiodic=False
        )
        assert np.isfinite(wrapped)
        assert wrapped != clipped

    def test_sum_from_integral_image_3D_nonperiodic_clamp(self):
        fgrid = np.arange(64, dtype=np.float64).reshape(4, 4, 4)
        iimg = get_integral_image3D(fgrid)
        clipped = sum_from_integral_image_3D(
            iimg, -2, 2, -1, 1, -1, 1, xperiodic=False, yperiodic=False, zperiodic=False
        )
        expected = sum_from_integral_image_3D(
            iimg, 0, 2, 0, 1, 0, 1, xperiodic=False, yperiodic=False, zperiodic=False
        )
        assert np.isclose(clipped, expected)

    def test_sum_from_integral_image_3D_py_func(self):
        fgrid = np.arange(27, dtype=np.float64).reshape(3, 3, 3)
        iimg = get_integral_image3D(fgrid)
        wrapped = sum_from_integral_image_3D(
            iimg, -1, 1, 1, 2, 1, 2, xperiodic=True, yperiodic=True, zperiodic=True
        )
        py_result = sum_from_integral_image_3D.py_func(
            iimg, -1, 1, 1, 2, 1, 2, xperiodic=True, yperiodic=True, zperiodic=True
        )
        assert np.isclose(wrapped, py_result)

    def test_get_volume_enclosing_box_3D(self):
        dgrid = np.random.rand(3, 3, 3)
        idgrid = get_integral_image3D(dgrid)
        volume = get_volume_enclosing_box_3D(1.0, 1.0, 1.0, dgrid, idgrid, minpart=1)
        assert volume.shape == (3, 3, 3)

    def test_get_volume_enclosing_box_3D_field(self):
        dgrid = np.full((3, 3, 3), 0.5, dtype=np.float64)
        idgrid = get_integral_image3D(dgrid)
        fgrid = np.arange(27, dtype=np.float64).reshape(3, 3, 3)
        ifgrid = get_integral_image3D(fgrid)
        fvolume = get_volume_enclosing_box_3D(
            3.0,
            3.0,
            3.0,
            dgrid,
            idgrid,
            minpart=1,
            xperiodic=False,
            yperiodic=False,
            zperiodic=False,
            ifgrid=ifgrid,
        )
        assert fvolume.shape == (3, 3, 3)
        assert np.all(np.isfinite(fvolume))

    def test_get_volume_enclosing_box_3D_py_func_field(self):
        dgrid = np.array(
            [[[10.0, 0.2], [0.2, 10.0]], [[0.2, 0.2], [10.0, 0.2]]], dtype=np.float64
        )
        idgrid = get_integral_image3D(dgrid)
        fgrid = np.arange(8, dtype=np.float64).reshape(2, 2, 2)
        ifgrid = get_integral_image3D(fgrid)
        fvolume_py = get_volume_enclosing_box_3D.py_func(
            2.0,
            2.0,
            2.0,
            dgrid,
            idgrid,
            minpart=1,
            xperiodic=False,
            yperiodic=False,
            zperiodic=False,
            ifgrid=ifgrid,
        )
        assert fvolume_py.shape == (2, 2, 2)
        assert np.all(np.isfinite(fvolume_py))

    def test_get_volume_enclosing_box_3D_small_density(self):
        dgrid = np.full((2, 2, 2), 0.2, dtype=np.float64)
        idgrid = get_integral_image3D(dgrid)
        volume = get_volume_enclosing_box_3D(
            2.0,
            2.0,
            2.0,
            dgrid,
            idgrid,
            minpart=1,
            xperiodic=True,
            yperiodic=True,
            zperiodic=True,
        )
        assert volume.shape == (2, 2, 2)
        assert np.all(np.isfinite(volume))

    def test_gridSPH2D(self):
        x = np.linspace(0.1, 0.9, 9)
        y = np.linspace(0.1, 0.9, 9)
        x2d, y2d = np.meshgrid(x, y)
        dgrid = gridSPH2D(1.0, 3, x2d.ravel(), y2d.ravel(), minpart=1)
        assert dgrid.shape == (3, 3)
        dgrid = gridSPH2D(
            [1.0, 1.0], 3, x2d.ravel(), y2d.ravel(), minpart=1, periodic=[True, True]
        )
        assert dgrid.shape == (3, 3)

    def test_gridSPH2D_field(self):
        x = np.linspace(0.1, 0.9, 9)
        y = np.linspace(0.1, 0.9, 9)
        x2d, y2d = np.meshgrid(x, y)
        f = np.cos(2 * np.pi * x2d.ravel())
        fgrid = gridSPH2D(1.0, 3, x2d.ravel(), y2d.ravel(), minpart=1, f=f)
        assert fgrid.shape == (3, 3)
        assert np.all(np.isfinite(fgrid))

    def test_gridSPH3D(self):
        x = np.linspace(0.2, 0.8, 8)
        y = np.linspace(0.2, 0.8, 8)
        z = np.linspace(0.2, 0.8, 8)
        x3d, y3d, z3d = np.meshgrid(x, y, z, indexing="ij")
        dgrid = gridSPH3D(1.0, 2, x3d.ravel(), y3d.ravel(), z3d.ravel(), minpart=1)
        assert dgrid.shape == (2, 2, 2)
        dgrid = gridSPH3D(
            [1.0, 1.0, 1.0],
            2,
            x3d.ravel(),
            y3d.ravel(),
            z3d.ravel(),
            minpart=1,
            periodic=[True, True, True],
        )
        assert dgrid.shape == (2, 2, 2)

    def test_gridSPH3D_field(self):
        x = np.linspace(0.2, 0.8, 8)
        y = np.linspace(0.2, 0.8, 8)
        z = np.linspace(0.2, 0.8, 8)
        x3d, y3d, z3d = np.meshgrid(x, y, z, indexing="ij")
        f = np.sin(2 * np.pi * x3d.ravel())
        fgrid = gridSPH3D(1.0, 2, x3d.ravel(), y3d.ravel(), z3d.ravel(), minpart=1, f=f)
        assert fgrid.shape == (2, 2, 2)
        assert np.all(np.isfinite(fgrid))
