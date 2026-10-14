=======================
MPI SPH Grid Estimation
=======================

This module provides distributed-memory (MPI) implementations of SPH
grid interpolation. Data is partitioned across MPI ranks along the x-axis,
with optional buffer regions to ensure correct neighbor evaluation.

mpi_sph4grid2D
==============

Distributed SPH interpolation in 2D using MPI.

Algorithm
---------

1. Pack particle data into MPI-distributed structure.
2. Sort and distribute particles along x-axis.
3. Add buffer regions for neighbor completeness.
4. Build local SPH grid via ``sph4grid2D``.

Example
-------

.. code-block:: python

    # Assuming correct initialisation of MPI object, see tutorial.
    import fiesta

    boxsize=100
    ngrid=256
    
    # compute density field
    dgrid = fiesta.sph.mpi_sph4grid2D(
        x, y,
        boxsize,
        ngrid,
        MPI,
        k=20
    )

    # compute field
    fgrid = fiesta.sph.mpi_sph4grid2D(
        x, y, 
        boxsize,
        ngrid,
        MPI,
        f=f,
        k=20
    )

mpi_sph4grid3D
===============

MPI-parallel SPH interpolation in 3D.

Algorithm
---------

1. Pack particle data into MPI format.
2. Distribute particles along x-axis slabs.
3. Apply buffer regions for neighbor consistency.
4. Compute local SPH grid using ``sph4grid3D``.

Important Notes
---------------

- ``buffer_length`` must be smaller than slab width.
- Periodic boundaries require consistent global box size.
- Load imbalance may occur for clustered particle distributions.

Example
-------

.. code-block:: python

    # Assuming correct initialisation of MPI object, see tutorial.
    import fiesta

    boxsize=[200, 200, 200]
    ngrid=128

    # compute density field
    dgrid = fiesta.sph.mpi_sph4grid3D(
        x, y, z,
        boxsize,
        ngrid,
        MPI,
        k=32
    )

    # compute field
    fgrid = fiesta.sph.mpi_sph4grid3D(
        x, y, z,
        boxsize,
        ngrid,
        MPI,
        f=f,
        k=32
    )
