DTFE
====

DTFE related functions for computing the Delaunay tesselation and interpolation densities and 
field values. MPI related classes take specific care to only consider points within complete
tesselations and distribute border points to adjacent slabs to complete DTFE computations
accurately without approximations.

DTFE Objects
------------

.. autoclass:: fiesta.dtfe.Delaunay2D
    :members: 
    :undoc-members:

.. autoclass:: fiesta.dtfe.DelaunayMerger2D
    :members: 
    :undoc-members:

.. autoclass:: fiesta.dtfe.Delaunay3D
    :members: 
    :undoc-members:

.. autoclass:: fiesta.dtfe.DelaunayMerger3D
    :members: 
    :undoc-members:

DTFE to Grid
------------

.. autofunction:: fiesta.dtfe.delaunay_density4grid2D
.. autofunction:: fiesta.dtfe.delaunay_field4grid2D
.. autofunction:: fiesta.dtfe.delaunay_density4grid3D
.. autofunction:: fiesta.dtfe.delaunay_field4grid3D


MPI DTFE to Grid
----------------

.. autofunction:: fiesta.dtfe.mpi_delaunay_density4grid2D
.. autofunction:: fiesta.dtfe.mpi_delaunay_field4grid2D
.. autofunction:: fiesta.dtfe.mpi_delaunay_density4grid3D
.. autofunction:: fiesta.dtfe.mpi_delaunay_field4grid3D