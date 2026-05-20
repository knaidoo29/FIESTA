================
DTFE onto a Grid
================

This section provides tutorials and usage examples for the DTFE 
grid interpolation routines.

Overview
========

The following functions are available:

* ``delaunay_density4grid2D`` -- 2D density estimation on a grid
* ``delaunay_field4grid2D`` -- 2D field interpolation on a grid
* ``delaunay_density4grid3D`` -- 3D density estimation on a grid
* ``delaunay_field4grid3D`` -- 3D field interpolation on a grid

These routines:

* Construct Delaunay tessellations from particle data
* Estimate densities or interpolate scalar fields
* Support periodic boundaries
* Support partitioned tessellations for large datasets
* Support pixel/voxel subsampling
* Return structured grid outputs

2D Density Estimation
=====================

Basic Example
-------------

Generate a 2D density field from particle coordinates.

.. code-block:: python

    import numpy as np
    import fiesta

    # Random particle positions
    npts = 1000

    x = np.random.uniform(0.0, 100.0, npts)
    y = np.random.uniform(0.0, 100.0, npts)

    dens = fiesta.dtfe.delaunay_density4grid2D(
        x,
        y,
        boxsize=100.0,
        ngrid=128,
    )

The returned array has shape:

.. code-block:: python

    print(dens.shape)

::

    (128, 128)

Using Mass Weights
------------------

Mass weights can be supplied for weighted density estimation.

.. code-block:: python

    mass = np.random.uniform(1.0, 10.0, npts)

    dens = fiesta.dtfe.delaunay_density4grid2D(
        x,
        y,
        boxsize=100.0,
        ngrid=128,
        mass=mass,
    )

Non-Square Grids
----------------

Different box sizes and grid dimensions may be used.

.. code-block:: python

    dens = fiesta.dtfe.delaunay_density4grid2D(
        x,
        y,
        boxsize=[200.0, 100.0],
        ngrid=[256, 128],
    )

Periodic Boundaries
-------------------

Enable periodic boundaries independently per axis.

.. code-block:: python

    dens = fiesta.dtfe.delaunay_density4grid2D(
        x,
        y,
        boxsize=100.0,
        ngrid=128,
        periodic=[True, False],
    )

Subsampling
-----------

Subsampling improves pixel averaging by evaluating multiple
points per pixel.

.. code-block:: python

    dens = fiesta.dtfe.delaunay_density4grid2D(
        x,
        y,
        boxsize=100.0,
        ngrid=128,
        subsampling=4,
    )

This evaluates ``4 x 4 = 16`` samples per pixel.

Returning Grid Coordinates
--------------------------

.. code-block:: python

    dens, xgrid, ygrid = fiesta.dtfe.delaunay_density4grid2D(
        x,
        y,
        boxsize=100.0,
        ngrid=128,
        outputgrid=True,
    )

Flattened Output
----------------

.. code-block:: python

    dens = fiesta.dtfe.delaunay_density4grid2D(
        x,
        y,
        boxsize=100.0,
        ngrid=128,
        flatten=True,
    )

The output shape becomes:

.. code-block:: python

    print(dens.shape)

::

    (16384,)

Large Dataset Partitioning
--------------------------

Partitioning reduces memory usage for large datasets.

.. code-block:: python

    dens = fiesta.dtfe.delaunay_density4grid2D(
        x,
        y,
        boxsize=100.0,
        ngrid=256,
        partition=4,
    )

This creates:

::

    4 x 4 = 16

internal tessellations.

2D Field Interpolation
======================

Basic Example
-------------

Interpolate scalar field values onto a grid.

.. code-block:: python

    f = np.sin(x / 10.0) + np.cos(y / 10.0)

    field = fiesta.dtfe.delaunay_field4grid2D(
        x,
        y,
        f,
        boxsize=100.0,
        ngrid=128,
    )

Returning Grid Coordinates
--------------------------

.. code-block:: python

    field, xgrid, ygrid = fiesta.dtfe.delaunay_field4grid2D(
        x,
        y,
        f,
        boxsize=100.0,
        ngrid=128,
        outputgrid=True,
    )

Using Partitions
----------------

.. code-block:: python

    field = fiesta.dtfe.delaunay_field4grid2D(
        x,
        y,
        f,
        boxsize=100.0,
        ngrid=256,
        partition=[2, 2],
    )

3D Density Estimation
=====================

Basic Example
-------------

.. code-block:: python

    npts = 5000

    x = np.random.uniform(0.0, 100.0, npts)
    y = np.random.uniform(0.0, 100.0, npts)
    z = np.random.uniform(0.0, 100.0, npts)

    dens = fiesta.dtfe.delaunay_density4grid3D(
        x,
        y,
        z,
        boxsize=100.0,
        ngrid=64,
    )

Output shape:

.. code-block:: python

    print(dens.shape)

::

    (64, 64, 64)

Anisotropic Volumes
-------------------

.. code-block:: python

    dens = fiesta.dtfe.delaunay_density4grid3D(
        x,
        y,
        z,
        boxsize=[200.0, 100.0, 50.0],
        ngrid=[128, 64, 32],
    )

Periodic Boundaries
-------------------

.. code-block:: python

    dens = fiesta.dtfe.delaunay_density4grid3D(
        x,
        y,
        z,
        boxsize=100.0,
        ngrid=64,
        periodic=[True, True, False],
    )

3D Partitioning
---------------

Large 3D tessellations are computationally expensive.

Partitioning is strongly recommended.

.. code-block:: python

    dens = fiesta.dtfe.delaunay_density4grid3D(
        x,
        y,
        z,
        boxsize=100.0,
        ngrid=128,
        partition=[2, 2, 2],
    )

This generates:

::

    2 x 2 x 2 = 8

internal tessellations.

3D Field Interpolation
======================

Basic Example
-------------

.. code-block:: python

    f = np.sin(x / 10.0) + np.cos(y / 10.0) + z / 50.0

    field = fiesta.dtfe.delaunay_field4grid3D(
        x,
        y,
        z,
        f,
        boxsize=100.0,
        ngrid=64,
    )

Returning Coordinates
---------------------

.. code-block:: python

    field, xgrid, ygrid, zgrid = fiesta.dtfe.delaunay_field4grid3D(
        x,
        y,
        z,
        f,
        boxsize=100.0,
        ngrid=64,
        outputgrid=True,
    )

Advanced Usage
==============

Output Exterior Information
---------------------------

Exterior information is useful for debugging tessellation coverage.

.. code-block:: python

    dens, exterior_border, pixID, count = fiesta.dtfe.delaunay_density4grid2D(
        x,
        y,
        boxsize=100.0,
        ngrid=128,
        outputexterior=True,
    )

Returned values:

* ``exterior_border`` -- tessellation border information
* ``pixID`` -- undersampled pixel indices
* ``count`` -- sample counts per pixel

Disabling Normalisation
-----------------------

.. code-block:: python

    dens = fiesta.dtfe.delaunay_density4grid2D(
        x,
        y,
        boxsize=100.0,
        ngrid=128,
        normalise=False,
    )

Buffer Regions
--------------

The ``fbuffer`` parameter controls periodic boundary buffer size.

.. code-block:: python

    dens = fiesta.dtfe.delaunay_density4grid3D(
        x,
        y,
        z,
        boxsize=100.0,
        ngrid=64,
        fbuffer=0.3,
    )

Performance Recommendations
===========================

2D Recommendations
------------------

* Use ``partition > 1`` for datasets larger than ``10^5`` particles
* Use ``subsampling = 1`` for fastest performance
* Use ``subsampling >= 2`` for smoother fields

3D Recommendations
------------------

* Always partition large datasets
* Memory usage scales rapidly with particle count
* ``partition=[2,2,2]`` or larger is recommended for dense volumes
* Use lower grid resolutions when prototyping

Example Workflow
================

Complete 2D workflow example:

.. code-block:: python

    import numpy as np

    # Generate particle distribution
    npts = 10000

    x = np.random.uniform(0.0, 200.0, npts)
    y = np.random.uniform(0.0, 200.0, npts)

    # Generate density field
    dens = fiesta.dtfe.delaunay_density4grid2D(
        x,
        y,
        boxsize=200.0,
        ngrid=256,
        partition=4,
        periodic=True,
        subsampling=2,
    )

    print(dens.shape)

Expected output:

::

    (256, 256)

API Reference
=============

.. autofunction:: fiesta.dtfe.delaunay_density4grid2D
    :no-index:

.. autofunction:: fiesta.dtfe.delaunay_field4grid2D
    :no-index:

.. autofunction:: fiesta.dtfe.delaunay_density4grid3D
    :no-index:

.. autofunction:: fiesta.dtfe.delaunay_field4grid3D
    :no-index: