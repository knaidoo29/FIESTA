========================================
MPI Bilinear and Trilinear Interpolation
========================================

The MPI interpolation functions provide distributed interpolation from
MPI-decomposed grids onto arbitrary particle or coordinate positions.

The interpolation is performed locally on each MPI slab after:

* distributing query points to the appropriate MPI rank,
* exchanging ghost grid cells between neighboring ranks,
* extending slab boundaries for interpolation consistency.

The module includes:

* ``mpi_bilinear`` — distributed bilinear interpolation for 2D grids
* ``mpi_trilinear`` — distributed trilinear interpolation for 3D grids

mpi_bilinear
=============

Overview
--------

``mpi_bilinear`` performs distributed bilinear interpolation on a 2D grid.

The grid is decomposed along the x-axis across MPI ranks.

Point Redistribution
--------------------

Interpolation coordinates are redistributed to the MPI slab containing
their x position.

This minimizes unnecessary communication and ensures interpolation is
performed locally.

You should either define these points only on rank 0 or ensure they are 
correctly distributed across ranks. If points are only defined on rank 0, 
they will be redistributed to the correct ranks for interpolation. If points 
are defined on all ranks, they should already be distributed according to 
the MPI decomposition. The function ``fiesta.coords.distribute_points_by_x`` 
can be used to distribute points based on their x coordinates.

Example: 2D Distributed Interpolation
-------------------------------------

.. code-block:: python

    # Assuming correct initialisation of MPI object, see tutorial.
    import numpy as np
    import shift
    import fiesta

    BOXSIZE = 1.0
    NGRID = 128

    # Local slab, should follow the shape of the MPI decomposition.
    # This defines the expected local x and y grid for each rank.
    x2D, y2D = shift.cart.mpi_grid2D(BOXSIZE, NGRID, MPI) 

    fgrid = np.random.random(np.shape(x2D))

    if MPI.rank == 0:
        # Query points, defined only on rank 0 for redistribution.
        x = np.random.random(1000)
        y = np.random.random(1000)
    else:
        x, y = None, None

    # Perform distributed bilinear interpolation.
    x, y, f = mpiinterp.mpi_bilinear(fgrid, NGRID, BOXSIZE, x, y, MPI, periodic=True)

    # x, y, f now contain the interpolated values for points in the local slab of each rank.

Shape Validation
----------------

The function validates whether the local slab shape matches the expected
MPI decomposition.

If not:

.. code-block:: text

    ERROR: Shape of fgrid does not match expectation for distributed array

is printed on rank 0.

mpi_trilinear
=============

Overview
--------

``mpi_trilinear`` performs distributed trilinear interpolation on 3D grids.

The domain decomposition occurs along x only.

Example: 3D Distributed Interpolation
-------------------------------------

.. code-block:: python

    # Assuming correct initialisation of MPI object, see tutorial.
    import numpy as np
    import shift
    import fiesta

    BOXSIZE = 1.0
    NGRID = 128

    # Local slab, should follow the shape of the MPI decomposition.
    # This defines the expected local x and y grid for each rank.
    x3D, y3D, z3D = shift.cart.mpi_grid3D(BOXSIZE, NGRID, MPI) 

    fgrid = np.random.random(np.shape(x3D))

    if MPI.rank == 0:
        # Query points, defined only on rank 0 for redistribution.
        x = np.random.random(1000)
        y = np.random.random(1000)
        z = np.random.random(1000)
    else:
        x, y, z = None, None, None

    # Perform distributed trilinear interpolation.
    x, y, z, f = mpiinterp.mpi_trilinear(fgrid, NGRID, BOXSIZE, x, y, z, MPI, periodic=True)

    # x, y, z, f now contain the interpolated values for points in the local slab of each rank

API Reference
=============

.. autofunction:: fiesta.interp.mpi_bilinear
    :noindex:

.. autofunction:: fiesta.interp.mpi_trilinear
    :noindex: