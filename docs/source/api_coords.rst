Coords
======

This manages data points in 1D/2D/3D dimensions, with particular attention in how to distribute
particles to the correct slab when using MPI.

Points
------

.. autofunction:: fiesta.coords.x2points
.. autofunction:: fiesta.coords.points2x
.. autofunction:: fiesta.coords.xy2points
.. autofunction:: fiesta.coords.points2xy
.. autofunction:: fiesta.coords.xyz2points
.. autofunction:: fiesta.coords.points2xyz
.. autofunction:: fiesta.coords.coord2points

.. autofunction:: fiesta.coords.split_limits_by_grid
.. autofunction:: fiesta.coords.mpi_find_range_2D
.. autofunction:: fiesta.coords.mpi_find_range_3D
.. autofunction:: fiesta.coords.mpi_open_2D
.. autofunction:: fiesta.coords.mpi_open_3D
.. autofunction:: fiesta.coords.check_coords_at_MPI_0
.. autofunction:: fiesta.coords.distribute_points_by_x

.. autoclass:: fiesta.coords.MPI_SortByX
    :members: 
    :undoc-members: