=====================
Distributing with MPI
=====================

``FIESTA`` can compute distributed job using MPI. This requires installing ``mpi4py``
and following the instructions carefully so that the distribution is correctly installed
with the MPI distribution on the user's computer. Typically this doesn't matter too much 
unless you are installing ``mpi4py`` on a HPC.

Running MPI Jobs
================

MPI scripts should be launched using ``mpirun`` or ``mpiexec``.

Example:

.. code-block:: bash

   mpirun -np 4 python example.py

This launches the script using 4 MPI processes.

Creating the MPI Object
=======================

FIESTA uses the MPI utilities from ``shift-fft``.

Create the MPI object as follows:

.. code-block:: python

   import sys
   from os import environ

   # Set thread environment variables FIRST
   N_THREADS = '1'
   environ['OMP_NUM_THREADS'] = N_THREADS
   environ['OPENBLAS_NUM_THREADS'] = N_THREADS
   environ['MKL_NUM_THREADS'] = N_THREADS
   environ['VECLIB_MAXIMUM_THREADS'] = N_THREADS
   environ['NUMEXPR_NUM_THREADS'] = N_THREADS

   from shift.mpiutils import MPI

   MPI = MPI()

For more details follow instructions `here <https://shift.readthedocs.io/en/latest/mpiutils.html>`_.

Distributing Points
===================

Points can be distributed across MPI ranks using the ``fiesta.coords.distribute_points_by_x`` function. 
By distributing the points based on their x coordinates, we can ensure that each rank only processes the 
points that fall within its assigned slab of the simulation box.

Points loaded in at rank 0 can be distributed to the correct ranks using:

.. code-block:: python

   # Let's for example only define points on rank 0
   if MPI.rank == 0:
       x = np.random.uniform(0, 1000, N)
       y = np.random.uniform(0, 1000, N)
       z = np.random.uniform(0, 1000, N)
   else:
       x = None
       y = None
       z = None

   if MPI.rank == 0:
      # Now construct the data array
      # For 2D points:
      data = fiesta.coords.xyz2points(x, y)

      # For 3D points:
      data = fiesta.coords.xyz2points(x, y, z)

      # For ND points, or 2D/3D points with additional fields (e.g. f1, f2, ...):
      data = fiesta.coords.coord2points([x, y, z, f1, f2, ...])

      # Now distribute the points across ranks based on their x coordinates
      data = fiesta.coords.distribute_points_by_x(data, boxsize=1000.0, ngrid=256, origin=0.0, MPI=MPI)
   else:
      # create an empty data array on other ranks to receive the distributed points
      data = None
   
   x = data[:, 0]
   y = data[:, 1]
   z = data[:, 2] # For 3D points
   f1 = data[:, 3] # For additional fields ... etc.

Points loaded in on all ranks, can be distributed to the correct ranks using the same function, however
in practice this is not recommended as it will involve processor^2 communications. It is better to load 
points in on rank 0 and distribute them to the correct ranks using the above method or load the correct
points for each slab. If you do load points in on all ranks, ensure they are not duplicated across ranks.

.. code-block:: python

   # Let's for example define points on all ranks
   x = np.random.uniform(0, 1000, N)
   y = np.random.uniform(0, 1000, N)
   z = np.random.uniform(0, 1000, N)

   # Now construct the data array
   # For 2D points:
   data = fiesta.coords.xyz2points(x, y)

   # For 3D points:
   data = fiesta.coords.xyz2points(x, y, z)

   # For ND points, or 2D/3D points with additional fields (e.g. f1, f2, ...):
   data = fiesta.coords.coord2points([x, y, z, f1, f2, ...])

   # Now distribute the points across ranks based on their x coordinates
   data = fiesta.coords.distribute_points_by_x(data, boxsize=1000.0, ngrid=256, origin=0.0, MPI=MPI)

   x = data[:, 0]
   y = data[:, 1]
   z = data[:, 2] # For 3D points
   f1 = data[:, 3] # For additional fields ... etc.

Distributed Grids
=================

Grids are distributed across MPI ranks using the MPI decomposition. Each rank owns a slab of the grid in 
the x-direction. The shape of the local grid on each rank should match the expected shape based on the MPI 
decomposition. This can be constructed using the ``shift`` grid utilities, for example:

.. code-block:: python

   import shift

   BOXSIZE = 1.0
   NGRID = 128

   # Local slab, should follow the shape of the MPI decomposition.
   # This defines the expected local x and y grid for each rank.
   x2D, y2D = shift.cart.mpi_grid2D(BOXSIZE, NGRID, MPI)

   # For 3D grids:
   x3D, y3D, z3D = shift.cart.mpi_grid3D(BOXSIZE, NGRID, MPI)

   # For ND grids, you can access the local grid for the xaxis as:
   xedges, x1D = shift.cart.mpi_grid1D(BOXSIZE, NGRID, MPI)

   # The minimum and maximum x coordinates for the local slab can be obtained from xedges:
   x_min = xedges[0]
   x_max = xedges[-1]
   # This can be used to determine which points fall within the local slab.

MPI Tutorials
=============

.. toctree::
   :maxdepth: 1
   
   mpi_tutorial_p2g
   mpi_tutorial_dtfe
   mpi_tutorial_sph_grid
   mpi_tutorial_gridsph
   mpi_tutorial_diff
   mpi_tutorial_interp

