Boundary
========

These functions are responsible for providing boundary particles outside of a box, either using particle
repetition when the data has periodic boundary conditions or random buffer particles. The MPI implementation,
using slab decomposition, will pass boundary and internal buffer particles between slabs.

Randoms
-------

.. autofunction:: fiesta.boundary.buffer_random_2D
.. autofunction:: fiesta.boundary.buffer_random_3D

.. autofunction:: fiesta.boundary.mpi_buffer_random_2D
.. autofunction:: fiesta.boundary.mpi_buffer_random_3D
.. autofunction:: fiesta.boundary.mpi_buffer_random_utils

Periodic
--------

.. autofunction:: fiesta.boundary.buffer_periodic
.. autofunction:: fiesta.boundary.subbox_buffer_periodic
.. autofunction:: fiesta.boundary.buffer_periodic_2D
.. autofunction:: fiesta.boundary.subbox_buffer_periodic_2D
.. autofunction:: fiesta.boundary.buffer_periodic_3D
.. autofunction:: fiesta.boundary.subbox_buffer_periodic_3D

.. autofunction:: fiesta.boundary.mpi_buffer_internal_2D
.. autofunction:: fiesta.boundary.mpi_buffer_internal_3D
.. autofunction:: fiesta.boundary.mpi_buffer_periodic_2D
.. autofunction:: fiesta.boundary.mpi_buffer_periodic_3D
