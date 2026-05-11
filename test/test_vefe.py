import numpy as np
import pytest
from fiesta.vefe import (
    get_integral_image2D,
    get_integral_image3D,
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

    def test_get_volume_enclosing_box_2D(self):
        dgrid = np.random.rand(5, 5)
        idgrid = get_integral_image2D(dgrid)
        volume = get_volume_enclosing_box_2D(1.0, 5, dgrid, idgrid, minpart=1)
        assert volume.shape == (5, 5)

    def test_sum_from_integral_image_3D(self):
        iimg = np.random.rand(4, 4, 4)
        sum_val = sum_from_integral_image_3D(iimg, 1, 2, 1, 2, 1, 2)
        assert isinstance(sum_val, (float, np.ndarray))

    def test_get_volume_enclosing_box_3D(self):
        dgrid = np.random.rand(3, 3, 3)
        idgrid = get_integral_image3D(dgrid)
        volume = get_volume_enclosing_box_3D(1.0, 3, dgrid, idgrid, minpart=1)
        assert volume.shape == (3, 3, 3)