import numpy as np
import pytest
from fiesta.dtfe import (
    Delaunay2D,
    DelaunayMerger2D,
    Delaunay3D,
    DelaunayMerger3D,
    delaunay_density4grid2D,
    delaunay_field4grid2D,
    delaunay_density4grid3D,
    delaunay_field4grid3D,
)


class TestDTFE:
    def test_delaunay2d_setup(self):
        x = np.random.rand(10)
        y = np.random.rand(10)
        delaunay = Delaunay2D()
        delaunay.setup(x, y, boundary=[0, 1, 0, 1])
        assert delaunay.x is not None
        assert delaunay.y is not None

    def test_delaunay2d_construct(self):
        x = np.random.rand(10)
        y = np.random.rand(10)
        delaunay = Delaunay2D()
        delaunay.setup(x, y, boundary=[0, 1, 0, 1])
        delaunay.construct()
        # After construct, should have triangulation
        assert delaunay.simplices is not None

    def test_delaunay2d_get_area(self):
        x = np.random.rand(10)
        y = np.random.rand(10)
        delaunay = Delaunay2D()
        delaunay.setup(x, y, boundary=[0, 1, 0, 1])
        delaunay.construct()
        delaunay.get_area()
        assert delaunay.delaunay_area is not None

    def test_delaunay2d_get_dens(self):
        x = np.random.rand(10)
        y = np.random.rand(10)
        delaunay = Delaunay2D()
        delaunay.setup(x, y, boundary=[0, 1, 0, 1])
        delaunay.construct()
        delaunay.get_area()
        delaunay.get_dens()
        assert delaunay.point_dens is not None

    def test_delaunay2d_run(self):
        x = np.random.rand(10)
        y = np.random.rand(10)
        delaunay = Delaunay2D()
        delaunay.setup(x, y, boundary=[0, 1, 0, 1])
        delaunay.run()
        assert delaunay.simplices is not None
        assert delaunay.point_type is not None

    def test_delaunay2d_estimate(self):
        x = np.random.rand(10) * 0.6 + 0.2
        y = np.random.rand(10) * 0.6 + 0.2
        f = np.random.rand(10)
        delaunay = Delaunay2D()
        delaunay.setup(x, y, boundary=[0, 1, 0, 1])
        delaunay.run()
        delaunay.get_area()
        delaunay.get_dens()
        delaunay.set_field(f)
        est = delaunay.estimate(x[:5], y[:5], allow_border=True)
        assert len(est) == 5
        assert not np.all(np.isnan(est))

    def test_delaunay2d_set_field(self):
        x = np.random.rand(10)
        y = np.random.rand(10)
        f = np.random.rand(10)
        delaunay = Delaunay2D()
        delaunay.setup(x, y, boundary=[0, 1, 0, 1])
        delaunay.run()
        delaunay.set_field(f)
        assert delaunay.f0 is not None
        assert delaunay.dtfe_mode == "field"

    def test_delaunay2d_set_field2dens(self):
        x = np.random.rand(10)
        y = np.random.rand(10)
        delaunay = Delaunay2D()
        delaunay.setup(x, y, boundary=[0, 1, 0, 1])
        delaunay.run()
        delaunay.get_area()
        delaunay.get_dens()
        delaunay.set_field2dens()
        assert delaunay.dtfe_mode == "density"

    def test_delaunay3d_setup(self):
        x = np.random.rand(10)
        y = np.random.rand(10)
        z = np.random.rand(10)
        delaunay = Delaunay3D()
        delaunay.setup(x, y, z, boundary=[0, 1, 0, 1, 0, 1])
        assert delaunay.x is not None
        assert delaunay.y is not None
        assert delaunay.z is not None

    def test_delaunay3d_construct(self):
        x = np.random.rand(10)
        y = np.random.rand(10)
        z = np.random.rand(10)
        delaunay = Delaunay3D()
        delaunay.setup(x, y, z, boundary=[0, 1, 0, 1, 0, 1])
        delaunay.construct()
        # After construct, should have triangulation
        assert delaunay.simplices is not None

    def test_delaunay3d_run(self):
        x = np.random.rand(10)
        y = np.random.rand(10)
        z = np.random.rand(10)
        delaunay = Delaunay3D()
        delaunay.setup(x, y, z, boundary=[0, 1, 0, 1, 0, 1])
        delaunay.run()
        assert delaunay.simplices is not None
        assert delaunay.point_type is not None

    def test_delaunay3d_get_volume(self):
        x = np.random.rand(10)
        y = np.random.rand(10)
        z = np.random.rand(10)
        delaunay = Delaunay3D()
        delaunay.setup(x, y, z, boundary=[0, 1, 0, 1, 0, 1])
        delaunay.construct()
        delaunay.get_volume()
        assert delaunay.delaunay_volume is not None

    def test_delaunay3d_get_dens(self):
        x = np.random.rand(10)
        y = np.random.rand(10)
        z = np.random.rand(10)
        delaunay = Delaunay3D()
        delaunay.setup(x, y, z, boundary=[0, 1, 0, 1, 0, 1])
        delaunay.construct()
        delaunay.get_volume()
        delaunay.get_dens()
        assert delaunay.point_dens is not None

    def test_delaunay3d_estimate(self):
        x = np.random.rand(10) * 0.6 + 0.2
        y = np.random.rand(10) * 0.6 + 0.2
        z = np.random.rand(10) * 0.6 + 0.2
        f = np.random.rand(10)
        delaunay = Delaunay3D()
        delaunay.setup(x, y, z, boundary=[0, 1, 0, 1, 0, 1])
        delaunay.run()
        delaunay.get_volume()
        delaunay.get_dens()
        delaunay.set_field(f)
        valid = np.where(delaunay.simptypes != -1.0)[0]
        assert len(valid) > 0
        simplices = delaunay.simplices[valid[:5]]
        xq = np.mean(delaunay.x[simplices], axis=1)
        yq = np.mean(delaunay.y[simplices], axis=1)
        zq = np.mean(delaunay.z[simplices], axis=1)
        est = delaunay.estimate(xq, yq, zq, allow_border=True)
        assert len(est) == len(xq)
        assert not np.all(np.isnan(est))

    def test_delaunay_merger2d(self):
        x = np.random.rand(20)
        y = np.random.rand(20)
        delaunay1 = Delaunay2D()
        delaunay1.setup(x[:10], y[:10], boundary=[0, 0.5, 0, 1])
        delaunay1.run()
        delaunay1.get_area()
        delaunay1.get_dens()
        delaunay1.set_field2dens()

        delaunay2 = Delaunay2D()
        delaunay2.setup(x[10:], y[10:], boundary=[0.5, 1, 0, 1])
        delaunay2.run()
        delaunay2.get_area()
        delaunay2.get_dens()
        delaunay2.set_field2dens()

        merger = DelaunayMerger2D()
        merger.add_border(delaunay1.get_border())
        merger.add_border(delaunay2.get_border())

        merger.run([0, 2, 0, 2], apply_filter=True)
        est = merger.estimate(x[:5], y[:5])
        assert len(est) == 5

        delaunay1.f = None
        delaunay1.get_border()

    def test_delaunay_merger3d(self):
        x = np.random.rand(20)
        y = np.random.rand(20)
        z = np.random.rand(20)
        delaunay1 = Delaunay3D()
        delaunay1.setup(x[:10], y[:10], z[:10], boundary=[0, 0.5, 0, 1, 0, 1])
        delaunay1.run()
        delaunay1.get_volume()
        delaunay1.get_dens()
        delaunay1.set_field2dens()

        delaunay2 = Delaunay3D()
        delaunay2.setup(x[10:], y[10:], z[10:], boundary=[0.5, 1, 0, 1, 0, 1])
        delaunay2.run()
        delaunay2.get_volume()
        delaunay2.get_dens()
        delaunay2.set_field2dens()

        merger = DelaunayMerger3D()
        merger.add_border(delaunay1.get_border())
        merger.add_border(delaunay2.get_border())
        merger.run([0, 1, 0, 1, 0, 1])
        est = merger.estimate(x[:5], y[:5], z[:5])
        assert len(est) == 5

        delaunay1.f = None
        delaunay1.get_border()

    def test_delaunay_density4grid2d(self):
        x = np.random.rand(50)
        y = np.random.rand(50)
        dens = delaunay_density4grid2D(x, y, boxsize=1.0, ngrid=10)
        assert dens.shape == (10, 10)

    def test_delaunay_field4grid2d(self):
        x = np.random.rand(50)
        y = np.random.rand(50)
        f = np.random.rand(50)
        field = delaunay_field4grid2D(x, y, f, boxsize=1.0, ngrid=10)
        assert field.shape == (10, 10)

    def test_delaunay_density4grid3d(self):
        x = np.random.rand(50)
        y = np.random.rand(50)
        z = np.random.rand(50)
        dens = delaunay_density4grid3D(x, y, z, boxsize=1.0, ngrid=5)
        assert dens.shape == (5, 5, 5)

    def test_delaunay_field4grid3d(self):
        x = np.random.rand(50)
        y = np.random.rand(50)
        z = np.random.rand(50)
        f = np.random.rand(50)
        field = delaunay_field4grid3D(x, y, z, f, boxsize=1.0, ngrid=5)
        assert field.shape == (5, 5, 5)

    # Tests for partition > 1 coverage (2D)
    def test_delaunay_density4grid2d_partition(self):
        """Test 2D density with partition > 1."""
        x = np.random.rand(100)
        y = np.random.rand(100)
        dens = delaunay_density4grid2D(x, y, boxsize=1.0, ngrid=10, partition=2)
        assert dens.shape == (10, 10)
        assert np.all(np.isfinite(dens))

    def test_delaunay_density4grid2d_partition_2x2(self):
        """Test 2D density with 2x2 partition."""
        x = np.random.rand(150)
        y = np.random.rand(150)
        dens = delaunay_density4grid2D(x, y, boxsize=1.0, ngrid=10, partition=[2, 2])
        assert dens.shape == (10, 10)
        assert np.all(np.isfinite(dens))

    def test_delaunay_field4grid2d_partition(self):
        """Test 2D field with partition > 1."""
        x = np.random.rand(100)
        y = np.random.rand(100)
        f = np.random.rand(100)
        field = delaunay_field4grid2D(x, y, f, boxsize=1.0, ngrid=10, partition=2)
        assert field.shape == (10, 10)
        assert np.all(np.isfinite(field))

    # Tests for periodic boundaries (2D)
    def test_delaunay_density4grid2d_periodic(self):
        """Test 2D density with periodic boundaries."""
        x = np.random.rand(100)
        y = np.random.rand(100)
        dens = delaunay_density4grid2D(x, y, boxsize=1.0, ngrid=10, periodic=True)
        assert dens.shape == (10, 10)
        assert np.all(np.isfinite(dens))

    def test_delaunay_density4grid2d_periodic_list(self):
        """Test 2D density with list-based periodic."""
        x = np.random.rand(100)
        y = np.random.rand(100)
        dens = delaunay_density4grid2D(
            x, y, boxsize=1.0, ngrid=10, periodic=[True, False]
        )
        assert dens.shape == (10, 10)
        assert np.all(np.isfinite(dens))

    def test_delaunay_field4grid2d_periodic(self):
        """Test 2D field with periodic boundaries."""
        x = np.random.rand(100)
        y = np.random.rand(100)
        f = np.random.rand(100)
        field = delaunay_field4grid2D(x, y, f, boxsize=1.0, ngrid=10, periodic=True)
        field = delaunay_field4grid2D(
            x, y, f, boxsize=1.0, ngrid=10, periodic=[True, False]
        )
        field = delaunay_field4grid2D(
            x, y, f, boxsize=1.0, ngrid=10, periodic=[False, True]
        )
        assert field.shape == (10, 10)
        assert np.all(np.isfinite(field))
        field = delaunay_field4grid2D(
            x,
            y,
            f,
            boxsize=1.0,
            ngrid=10,
            periodic=[False, True],
            outputexterior=True,
            outputgrid=True,
        )
        field = delaunay_field4grid2D(
            x,
            y,
            f,
            boxsize=1.0,
            ngrid=10,
            periodic=[False, True],
            outputexterior=False,
            outputgrid=True,
        )

    # Tests for output flags (2D)
    def test_delaunay_density4grid2d_outputgrid(self):
        """Test 2D density with outputgrid=True."""
        x = np.random.rand(100)
        y = np.random.rand(100)
        result = delaunay_density4grid2D(x, y, boxsize=1.0, ngrid=10, outputgrid=True)
        assert len(result) == 3
        dens, x2d, y2d = result
        assert dens.shape == (10, 10)
        assert x2d.shape == (100,)
        assert y2d.shape == (100,)

    def test_delaunay_density4grid2d_outputexterior(self):
        """Test 2D density with outputexterior=True."""
        x = np.random.rand(100)
        y = np.random.rand(100)
        result = delaunay_density4grid2D(
            x, y, boxsize=1.0, ngrid=10, outputexterior=True
        )
        assert len(result) == 4
        dens, exterior_border, pixID, count = result
        assert dens.shape == (10, 10)
        assert exterior_border is not None
        assert len(pixID) >= 0
        assert len(count) >= 0

    def test_delaunay_density4grid2d_both_outputs(self):
        """Test 2D density with both outputgrid and outputexterior."""
        x = np.random.rand(100)
        y = np.random.rand(100)
        result = delaunay_density4grid2D(
            x, y, boxsize=1.0, ngrid=10, outputgrid=True, outputexterior=True
        )
        assert len(result) == 6
        dens, x2d, y2d, exterior_border, pixID, count = result
        assert dens.shape == (10, 10)

    # Tests for partition and periodic combination (2D)
    def test_delaunay_density4grid2d_partition_periodic(self):
        """Test 2D density with partition and periodic boundaries."""
        x = np.random.rand(150)
        y = np.random.rand(150)
        dens = delaunay_density4grid2D(
            x, y, boxsize=1.0, ngrid=10, partition=2, periodic=[True, False]
        )
        dens = delaunay_density4grid2D(
            x, y, boxsize=1.0, ngrid=10, partition=2, periodic=[False, True]
        )
        assert dens.shape == (10, 10)
        assert np.all(np.isfinite(dens))

    # Tests for 3D partition
    def test_delaunay_density4grid3d_partition(self):
        """Test 3D density with partition > 1."""
        x = np.random.rand(500)
        y = np.random.rand(500)
        z = np.random.rand(500)
        dens = delaunay_density4grid3D(x, y, z, boxsize=1.0, ngrid=5, partition=2)
        assert dens.shape == (5, 5, 5)
        assert np.all(np.isfinite(dens))

    def test_delaunay_density4grid3d_partition_2x2x2(self):
        """Test 3D density with 2x2x2 partition."""
        x = np.random.rand(200)
        y = np.random.rand(200)
        z = np.random.rand(200)
        dens = delaunay_density4grid3D(
            x, y, z, boxsize=1.0, ngrid=5, partition=[2, 2, 2]
        )
        assert dens.shape == (5, 5, 5)
        assert np.all(np.isfinite(dens))

    def test_delaunay_field4grid3d_partition(self):
        """Test 3D field with partition > 1."""
        x = np.random.rand(150)
        y = np.random.rand(150)
        z = np.random.rand(150)
        f = np.random.rand(150)
        field = delaunay_field4grid3D(x, y, z, f, boxsize=1.0, ngrid=5, partition=2)
        assert field.shape == (5, 5, 5)
        assert np.all(np.isfinite(field))

    # Tests for 3D periodic boundaries
    def test_delaunay_density4grid3d_periodic(self):
        """Test 3D density with periodic boundaries."""
        x = np.random.rand(150)
        y = np.random.rand(150)
        z = np.random.rand(150)
        dens = delaunay_density4grid3D(x, y, z, boxsize=1.0, ngrid=5, periodic=True)
        assert dens.shape == (5, 5, 5)
        assert np.all(np.isfinite(dens))

    def test_delaunay_density4grid3d_periodic_list(self):
        """Test 3D density with list-based periodic."""
        x = np.random.rand(150)
        y = np.random.rand(150)
        z = np.random.rand(150)
        dens = delaunay_density4grid3D(
            x, y, z, boxsize=1.0, ngrid=5, periodic=[True, False, True]
        )
        assert dens.shape == (5, 5, 5)
        assert np.all(np.isfinite(dens))

    def test_delaunay_field4grid3d_periodic(self):
        """Test 3D field with periodic boundaries."""
        x = np.random.rand(150)
        y = np.random.rand(150)
        z = np.random.rand(150)
        f = np.random.rand(150)
        field = delaunay_field4grid3D(x, y, z, f, boxsize=1.0, ngrid=5, periodic=True)
        assert field.shape == (5, 5, 5)
        assert np.all(np.isfinite(field))

    # Tests for 3D output flags
    def test_delaunay_density4grid3d_outputgrid(self):
        """Test 3D density with outputgrid=True."""
        x = np.random.rand(150)
        y = np.random.rand(150)
        z = np.random.rand(150)
        result = delaunay_density4grid3D(x, y, z, boxsize=1.0, ngrid=5, outputgrid=True)
        assert len(result) == 4
        dens, x3d, y3d, z3d = result
        assert dens.shape == (5, 5, 5)
        assert x3d.shape == (125,)

    def test_delaunay_density4grid3d_outputexterior(self):
        """Test 3D density with outputexterior=True."""
        x = np.random.rand(150)
        y = np.random.rand(150)
        z = np.random.rand(150)
        result = delaunay_density4grid3D(
            x, y, z, boxsize=1.0, ngrid=5, outputexterior=True
        )
        assert len(result) == 4
        dens, exterior_border, pixID, count = result
        assert dens.shape == (5, 5, 5)

    def test_delaunay_density4grid3d_both_outputs(self):
        """Test 3D density with both outputgrid and outputexterior."""
        x = np.random.rand(150)
        y = np.random.rand(150)
        z = np.random.rand(150)
        result = delaunay_density4grid3D(
            x, y, z, boxsize=1.0, ngrid=5, outputgrid=True, outputexterior=True
        )
        assert len(result) == 7
        dens, x3d, y3d, z3d, exterior_border, pixID, count = result
        assert dens.shape == (5, 5, 5)

    def test_delaunay_field4grid3d_outputexterior(self):
        """Test 3D field with outputexterior=True."""
        x = np.random.rand(150)
        y = np.random.rand(150)
        z = np.random.rand(150)
        f = np.random.rand(150)
        result = delaunay_field4grid3D(
            x, y, z, f, boxsize=1.0, ngrid=5, outputexterior=True
        )
        assert len(result) == 4
        field, exterior_border, pixID, count = result
        assert field.shape == (5, 5, 5)

    def test_delaunay_field4grid3d_both_outputs(self):
        """Test 3D field with both outputgrid and outputexterior."""
        x = np.random.rand(150)
        y = np.random.rand(150)
        z = np.random.rand(150)
        f = np.random.rand(150)
        result = delaunay_field4grid3D(
            x, y, z, f, boxsize=1.0, ngrid=5, outputgrid=True, outputexterior=True
        )
        assert len(result) == 7
        x = np.random.rand(150)
        y = np.random.rand(150)
        z = np.random.rand(150)
        f = np.random.rand(150)
        result = delaunay_field4grid3D(
            x, y, z, f, boxsize=1.0, ngrid=5, outputgrid=True, outputexterior=False
        )
        field, x3d, y3d, z3d = result
        assert field.shape == (5, 5, 5)

    # Tests for 3D partition and periodic combination
    def test_delaunay_density4grid3d_partition_periodic(self):
        """Test 3D density with partition and periodic boundaries."""
        x = np.random.rand(200)
        y = np.random.rand(200)
        z = np.random.rand(200)
        dens = delaunay_density4grid3D(
            x, y, z, boxsize=1.0, ngrid=5, partition=2, periodic=[True, False, True]
        )
        dens = delaunay_density4grid3D(
            x, y, z, boxsize=1.0, ngrid=5, partition=2, periodic=[False, True, False]
        )
        assert dens.shape == (5, 5, 5)
        assert np.all(np.isfinite(dens))

    def test_delaunay_field4grid3d_partition_periodic(self):
        """Test 3D field with partition and periodic boundaries."""
        x = np.random.rand(200)
        y = np.random.rand(200)
        z = np.random.rand(200)
        f = np.random.rand(200)
        field = delaunay_field4grid3D(
            x, y, z, f, boxsize=1.0, ngrid=5, partition=2, periodic=[True, False, True]
        )
        field = delaunay_field4grid3D(
            x, y, z, f, boxsize=1.0, ngrid=5, partition=2, periodic=[False, True, False]
        )
        assert field.shape == (5, 5, 5)
        assert np.all(np.isfinite(field))

    # Tests with list parameters for boxsize, ngrid, origin
    def test_delaunay_density4grid2d_list_params(self):
        """Test 2D density with list boxsize and ngrid."""
        x = np.random.rand(100) * 1.0
        y = np.random.rand(100) * 0.8
        dens = delaunay_density4grid2D(
            x, y, boxsize=[1.0, 0.8], ngrid=[10, 8], partition=2, periodic=True
        )
        assert dens.shape == (10, 8)
        assert np.all(np.isfinite(dens))
        dens = delaunay_density4grid2D(
            x, y, boxsize=[1.0, 0.8], ngrid=[10, 8], partition=2, periodic=False
        )
        assert dens.shape == (10, 8)
        assert np.all(np.isfinite(dens))

    def test_delaunay_field4grid2d_list_params(self):
        """Test 2D field with list boxsize and ngrid."""
        x = np.random.rand(100) * 1.0
        y = np.random.rand(100) * 0.8
        f = np.random.rand(100)
        field = delaunay_field4grid2D(
            x, y, f, boxsize=[1.0, 0.8], ngrid=[10, 8], partition=2
        )
        assert field.shape == (10, 8)
        assert np.all(np.isfinite(field))
        field = delaunay_field4grid2D(
            x, y, f, boxsize=[1.0, 0.8], ngrid=[10, 8], partition=[2, 2]
        )
        assert field.shape == (10, 8)
        assert np.all(np.isfinite(field))

    def test_delaunay_density4grid3d_list_params(self):
        """Test 3D density with list parameters."""
        x = np.random.rand(150) * 1.0
        y = np.random.rand(150) * 0.8
        z = np.random.rand(150) * 0.6
        dens = delaunay_density4grid3D(
            x, y, z, boxsize=[1.0, 0.8, 0.6], ngrid=[5, 4, 3], partition=2
        )
        assert dens.shape == (5, 4, 3)
        assert np.all(np.isfinite(dens))

    def test_delaunay_field4grid3d_list_params(self):
        """Test 3D field with list parameters."""
        x = np.random.rand(500) * 1.0
        y = np.random.rand(500) * 0.8
        z = np.random.rand(500) * 0.6
        f = np.random.rand(500)
        field = delaunay_field4grid3D(
            x, y, z, f, boxsize=[1.0, 0.8, 0.6], ngrid=[5, 4, 3], partition=2
        )
        assert field.shape == (5, 4, 3)
        assert np.all(np.isfinite(field))

    # Tests for list-based subsampling
    def test_delaunay_density4grid2d_list_subsampling(self):
        """Test 2D density with list-based subsampling."""
        x = np.random.rand(100)
        y = np.random.rand(100)
        dens = delaunay_density4grid2D(
            x, y, boxsize=[1.0, 1.0], ngrid=10, subsampling=[2, 2]
        )
        assert dens.shape == (10, 10)
        assert np.all(np.isfinite(dens))

    def test_delaunay_field4grid2d_list_subsampling(self):
        """Test 2D field with list-based subsampling."""
        x = np.random.rand(100)
        y = np.random.rand(100)
        f = np.random.rand(100)
        field = delaunay_field4grid2D(
            x, y, f, boxsize=[1.0, 1.0], ngrid=10, subsampling=[2, 2]
        )
        assert field.shape == (10, 10)
        assert np.all(np.isfinite(field))

    def test_delaunay_density4grid3d_list_subsampling(self):
        """Test 3D density with list-based subsampling."""
        x = np.random.rand(150)
        y = np.random.rand(150)
        z = np.random.rand(150)
        dens = delaunay_density4grid3D(
            x, y, z, boxsize=[1.0, 1.0, 1.0], ngrid=5, subsampling=[2, 2, 2]
        )
        assert dens.shape == (5, 5, 5)
        assert np.all(np.isfinite(dens))

    def test_delaunay_field4grid3d_list_subsampling(self):
        """Test 3D field with list-based subsampling."""
        x = np.random.rand(150)
        y = np.random.rand(150)
        z = np.random.rand(150)
        f = np.random.rand(150)
        field = delaunay_field4grid3D(
            x,
            y,
            z,
            f,
            boxsize=[1.0, 1.0, 1.0],
            ngrid=5,
            subsampling=[2, 2, 2],
            partition=[2, 2, 2],
        )
        assert field.shape == (5, 5, 5)
        assert np.all(np.isfinite(field))

    # Tests with list-based origin
    def test_delaunay_density4grid2d_list_origin(self):
        """Test 2D density with list-based origin."""
        x = np.random.rand(100) + 0.1
        y = np.random.rand(100) + 0.1
        dens = delaunay_density4grid2D(x, y, boxsize=1.0, ngrid=10, origin=[0.1, 0.1])
        assert dens.shape == (10, 10)
        assert np.all(np.isfinite(dens))

    def test_delaunay_field4grid2d_list_origin(self):
        """Test 2D field with list-based origin."""
        x = np.random.rand(100) + 0.1
        y = np.random.rand(100) + 0.1
        f = np.random.rand(100)
        field = delaunay_field4grid2D(x, y, f, boxsize=1.0, ngrid=10, origin=[0.1, 0.1])
        assert field.shape == (10, 10)
        assert np.all(np.isfinite(field))

    def test_delaunay_density4grid3d_list_origin(self):
        """Test 3D density with list-based origin."""
        x = np.random.rand(150) + 0.1
        y = np.random.rand(150) + 0.1
        z = np.random.rand(150) + 0.1
        dens = delaunay_density4grid3D(
            x, y, z, boxsize=1.0, ngrid=5, origin=[0.1, 0.1, 0.1]
        )
        assert dens.shape == (5, 5, 5)
        assert np.all(np.isfinite(dens))

    def test_delaunay_field4grid3d_list_origin(self):
        """Test 3D field with list-based origin."""
        x = np.random.rand(150) + 0.1
        y = np.random.rand(150) + 0.1
        z = np.random.rand(150) + 0.1
        f = np.random.rand(150)
        field = delaunay_field4grid3D(
            x, y, z, f, boxsize=1.0, ngrid=5, origin=[0.1, 0.1, 0.1]
        )
        assert field.shape == (5, 5, 5)
        assert np.all(np.isfinite(field))
