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

MPI Tutorials
=============

.. toctree::
  :maxdepth: 1

  mpi_tutorial_p2g
  mpi_tutorial_dtfe
  mpi_tutorial_sph_grid
  mpi_tutorial_gridsph
  mpi_tutorial_diff

