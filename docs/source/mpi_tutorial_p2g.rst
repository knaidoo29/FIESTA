====================
MPI Particle-to-Grid
====================

Introduction
============

FIESTA provides MPI-enabled particle-to-grid interpolation routines for
distributed-memory parallel computations. These functions allow very large
particle datasets to be assigned onto regular Cartesian grids using
multiple MPI processes.

The MPI routines currently support:

- ``mpi_part2grid2D``
- ``mpi_part2grid3D``

Parallelisation Strategy
========================

The MPI routines distribute particles across processes along the ``x`` axis.
Each MPI rank computes a local portion of the grid and communicates overlap
regions between neighbouring processes when required by higher-order
assignment schemes such as:

- ``CIC``
- ``TSC``
- ``PCS``

The ``NGP`` method does not require overlap communication.


2D Example
==========

The following example distributes particles across MPI processes and assigns
them onto a 2D grid using the ``CIC`` interpolation scheme.

.. code-block:: python

   # Assuming correct initialisation of MPI object, see tutorial.
   import numpy as np
   import fiesta

   # Generate particles on root rank
   if MPI.rank == 0:

       np.random.seed(0)

       npart = 100000

       x = np.random.uniform(0.0, 1.0, npart)
       y = np.random.uniform(0.0, 1.0, npart)

       # Particle weights
       f = np.ones(npart)

   else:

       x = None
       y = None
       f = None

   # Interpolate particles onto the grid
   fgrid = fiesta.mpi_part2grid2D(
       x,
       y,
       f,
       boxsize=1.0,
       ngrid=256,
       MPI=MPI,
       method='CIC',
       periodic=True
   )

   print(MPI.rank, fgrid.shape)

Run using:

.. code-block:: bash

   mpirun -np 4 python example_2d.py

Notes
-----

- Only rank 0 initially contains the particle data.
- The particles are automatically distributed across MPI ranks.
- Each rank returns its local grid slab.

3D Example
==========

The following example performs distributed 3D particle-to-grid interpolation.

.. code-block:: python

   # Assuming correct initialisation of MPI object, see tutorial.
   import numpy as np
   import fiesta

   if MPI.rank == 0:

       np.random.seed(0)

       npart = 500000

       x = np.random.uniform(0.0, 1.0, npart)
       y = np.random.uniform(0.0, 1.0, npart)
       z = np.random.uniform(0.0, 1.0, npart)

       f = np.ones(npart)

   else:

       x = None
       y = None
       z = None
       f = None

   fgrid = fiesta.mpi_part2grid3D(
       x,
       y,
       z,
       f,
       boxsize=1.0,
       ngrid=128,
       MPI=MPI,
       method='TSC',
       periodic=True
   )

   print(MPI.rank, fgrid.shape)

Run using:

.. code-block:: bash

   mpirun -np 8 python example_3d.py

Assignment Schemes
==================

The MPI routines support the following interpolation schemes:

+----------+----------------------------+
| Method   | Description                |
+==========+============================+
| ``NGP``  | Nearest Grid Point         |
+----------+----------------------------+
| ``CIC``  | Cloud-In-Cell              |
+----------+----------------------------+
| ``TSC``  | Triangular Shaped Cloud    |
+----------+----------------------------+
| ``PCS``  | Piecewise Cubic Spline     |
+----------+----------------------------+

Higher-order schemes produce smoother fields but require additional
inter-process communication.

Periodic Boundaries
===================

Periodic boundary conditions can be enabled using:

.. code-block:: python

   periodic=True

When periodic boundaries are enabled, overlap regions are exchanged across
domain boundaries.

Grid Distribution
=================

The returned grids are distributed across MPI ranks along the ``x`` axis.

For example:

.. code-block:: text

   Rank 0 -> x = [0 : nx0]
   Rank 1 -> x = [nx0 : nx1]
   Rank 2 -> x = [nx1 : nx2]

Each rank stores only its local section of the full grid.

Gathering the Full Grid
=======================

To reconstruct the complete grid on a single rank, use MPI gather operations.

Example:

.. code-block:: python

   grids = MPI.collect(fgrid)

   if MPI.rank == 0:

       full_grid = np.concatenate(grids, axis=0)

Memory Considerations
=====================

The MPI routines are intended for large distributed-memory simulations where
single-node memory is insufficient. These remain the fastest schemes to assign 
particles/fields to a grid.

API Reference
=============

.. autofunction:: fiesta.p2g.mpi_part2grid2D
    :no-index:

.. autofunction:: fiesta.p2g.mpi_part2grid3D
    :no-index: