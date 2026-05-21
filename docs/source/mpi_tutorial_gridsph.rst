===================
MPI Grid-SPH Method
===================

This module provides MPI-parallel implementations of the Grid-SPH algorithm.
The method combines:

- particle-to-grid mapping,
- integral-image acceleration,
- and MPI slab decomposition,

to efficiently compute SPH-like density and field estimates on very large
datasets.

mpi_gridSPH2D
=============

MPI-parallel 2D Grid-SPH density and field estimation. This function distributes 
the computational grid along the x-axis across MPI ranks. Each rank computes a 
local Grid-SPH estimate using ghost-cell padding (buffer regions) to preserve 
SPH consistency near slab boundaries.

Algorithm
---------

1. Split the x-axis grid across MPI ranks.
2. Deposit particles onto local grids.
3. Exchange ghost-cell boundary regions.
4. Compute 2D integral images.
5. Estimate adaptive SPH areas.
6. Convert areas into densities or field estimates.
7. Remove ghost-cell padding.

Example
-------

.. code-block:: python

    # Assuming correct initialisation of MPI object, see tutorial.
    import fiesta

    boxsize=100
    ngrid=512
    
    rho = fiesta.gridsph.mpi_gridSPH2D(x, y, boxsize, ngrid, MPI, minpart=8, buffer_size=4)

Field Estimation Example
------------------------

.. code-block:: python

    # Assuming correct initialisation of MPI object, see tutorial.
    import fiesta

    temp = fiesta.gridsph.mpi_gridSPH2D(x, y, [100, 100], 256, MPI, f=temperature, minpart=16)

mpi_gridSPH3D
=============

MPI-parallel 3D Grid-SPH density and field estimation.

Overview
--------

This function extends the Grid-SPH algorithm into three dimensions with
distributed-memory parallelism. The computational domain is divided into 
x-axis slabs, with buffer regions used to maintain SPH consistency across 
MPI boundaries.

Algorithm
---------

1. Split the x-axis grid across MPI ranks.
2. Deposit particles onto local grids.
3. Exchange ghost-cell boundary regions.
4. Compute 3D integral images.
5. Estimate adaptive SPH volumes.
6. Convert volumes into densities or field estimates.
7. Remove ghost-cell padding.

Example
-------

.. code-block:: python

    # Assuming correct initialisation of MPI object, see tutorial.
    import fiesta

    boxsize=200
    ngrid=256

    rho = fiesta.gridsph.mpi_gridSPH3D(x, y, z, boxsize, ngrid, MPI, minpart=32, buffer_size=6)

3D Field Example
----------------

.. code-block:: python

    # Assuming correct initialisation of MPI object, see tutorial.
    import fiesta

    temp = fiesta.gridsph.mpi_gridSPH3D(x, y, z, [200, 200, 200], 128, MPI, f=temperature, minpart=20)

Performance Notes
=================

Advantages
----------

- Efficient for extremely large particle datasets.
- Integral-image operations scale very well.
- Reduced memory compared to full-domain SPH methods.
- Better scalability than KDTree-based SPH.

Choosing ``buffer_size``
========================

The ``buffer_size`` parameter is critical for accurate MPI-SPH estimation.

Guidelines
----------

- Too small:
    SPH kernels become truncated near slab boundaries.

- Too large:
    Increased MPI communication and memory overhead.

Limitations
===========

- Slab decomposition occurs only along the x-axis.
- Load imbalance may appear for clustered particle distributions.
- Accuracy depends on grid resolution and ``minpart``.