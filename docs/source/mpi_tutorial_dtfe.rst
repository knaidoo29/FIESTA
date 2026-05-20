======================
MPI DTFE Grid Routines
======================

This tutorial explains the MPI-enabled DTFE routines for constructing density and interpolated fields on 2D and 3D grids.

The routines covered are:

* ``mpi_delaunay_density4grid2D``
* ``mpi_delaunay_field4grid2D``
* ``mpi_delaunay_density4grid3D``
* ``mpi_delaunay_field4grid3D``

These functions parallelise DTFE calculations by decomposing the simulation domain into slabs distributed across MPI ranks.

MPI Parallelisation Strategy
============================

The MPI versions divide the simulation box along the x-axis.

Each MPI rank:

1. Receives a slab of particles.
2. Computes a local DTFE.
3. Exchanges boundary simplices with neighbouring ranks.
4. Merges border tessellations.
5. Evaluates grid values inside its local domain.

This enables very large DTFE calculations that would not fit into memory on a single process.

mpi_delaunay_density4grid2D
===========================

Computes a 2D DTFE density field on a regular grid using MPI.

Minimal Example: Density
------------------------

.. code-block:: python

    # Assuming correct initialisation of MPI object, see tutorial.
    import numpy as np
    import fiesta

    N = 100000

    x = np.random.uniform(0, 1000, N)
    y = np.random.uniform(0, 1000, N)

    rho = fiesta.dtfe.mpi_delaunay_density4grid2D(
        x=x,
        y=y,
        boxsize=1000.0,
        ngrid=256,
        MPI=MPI,
        periodic=True
    )

Each rank owns only its slab.

Mass Weighting
--------------

.. code-block:: python

    # Assuming correct initialisation of MPI object, see tutorial.
    import fiesta

    field = fiesta.dtfe.mpi_delaunay_field4grid2D(
        x=x,
        y=y,
        f=temperature,
        mass=particle_mass,
        boxsize=1000,
        ngrid=256,
        MPI=MPI
    )

With Coordinate Output
----------------------

.. code-block:: python

    # Assuming correct initialisation of MPI object, see tutorial.
    import fiesta

    rho, x2d, y2d = fiesta.dtfe.mpi_delaunay_density4grid2D(
        x,
        y,
        boxsize=1000,
        ngrid=256,
        MPI=MPI,
        outputgrid=True
    )

mpi_delaunay_field4grid2D
=========================

Interpolates an arbitrary scalar field onto a 2D grid. Unlike the density routine, 
this function uses particle field values:

.. code-block:: python
    
    f

Example: Temperature Interpolation
----------------------------------

.. code-block:: python

    # Assuming correct initialisation of MPI object, see tutorial.
    import numpy as np
    import fiesta

    temperature = np.sin(x / 100.0) * np.cos(y / 100.0)

    field = fiesta.dtfe.mpi_delaunay_field4grid2D(
        x=x,
        y=y,
        f=temperature,
        boxsize=1000,
        ngrid=256,
        MPI=MPI,
        periodic=True
    )

mpi_delaunay_density4grid3D
===========================

Computes a full 3D DTFE density field. This is the MPI parallel version of 
volumetric density estimation.

Example
-------

.. code-block:: python

    # Assuming correct initialisation of MPI object, see tutorial.
    import numpy as np
    import fiesta

    N = 2_000_000

    x = np.random.uniform(0, 500, N)
    y = np.random.uniform(0, 500, N)
    z = np.random.uniform(0, 500, N)

    rho = fiesta.dtfe.mpi_delaunay_density4grid3D(
        x=x,
        y=y,
        z=z,
        boxsize=500.0,
        ngrid=256,
        MPI=MPI,
        periodic=True,
        subsampling=2
    )

mpi_delaunay_field4grid3D
=========================

Interpolates scalar fields in full 3D.

Examples:

* temperature
* velocity divergence
* metallicity
* pressure
* potential

Example
-------

.. code-block:: python
    
    # Assuming correct initialisation of MPI object, see tutorial.
    import numpy as np
    import fiesta

    velocity_divergence = np.random.normal(size=N)

    field = fiesta.dtfe.mpi_delaunay_field4grid3D(
        x=x,
        y=y,
        z=z,
        f=velocity_divergence,
        boxsize=500,
        ngrid=256,
        MPI=MPI,
        periodic=True
    )

Performance Tuning
==================

Increase MPI Ranks
------------------

Large DTFE calculations scale well with:

* more particles
* larger grids
* higher memory demand

Adjust partition
----------------

Higher values:

.. code-block:: python

    partition = 4

reduce local tessellation memory.

Reduce subsampling
------------------

Subsampling is expensive.

Fastest:

.. code-block:: python

    subsampling = 1

Higher accuracy:

.. code-block:: python
    
    subsampling = 2 # 4 or more, but larger values increase accuracy but increase computation time.

Common Pitfalls
===============

Mismatched Periodicity
----------------------

If your simulation is periodic:

.. code-block:: python
    
    periodic=True

must be enabled. Otherwise discontinuities appear at boundaries.

Memory Explosion in 3D
----------------------

3D Delaunay tessellations can become extremely large.

Recommendations:

* increase MPI ranks
* use smaller partitions
* reduce grid size
* lower subsampling

API Reference
=============

.. autofunction:: fiesta.dtfe.mpi_delaunay_density4grid2D
    :no-index:

.. autofunction:: fiesta.dtfe.mpi_delaunay_field4grid2D
    :no-index:

.. autofunction:: fiesta.dtfe.mpi_delaunay_density4grid3D
    :no-index:

.. autofunction:: fiesta.dtfe.mpi_delaunay_field4grid2D
    :no-index: