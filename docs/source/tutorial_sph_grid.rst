===============
SPH onto a Grid
===============

Overview
========

The functions ``sph4grid2D`` and ``sph4grid3D`` compute Smooth Particle Hydrodynamics (SPH)
estimates of either:

- density fields (if ``f is None``)
- scalar fields (if ``f`` is provided)

on a regular Cartesian grid.

They are built on top of the ``SPH2D`` and ``SPH3D`` KDTree-based SPH implementations
and use nearest-neighbour smoothing with a cubic kernel.

sph4grid2D
==========

Description
-----------

Estimates SPH density or field values on a **2D Cartesian grid**.

Algorithm
---------

1. Build KDTree using SPH2D
2. Assign particle masses
3. Construct regular grid
4. Define subsampling offsets within each pixel
5. For each offset:
   - shift grid positions
   - compute SPH estimate
6. Average over all subsamples
7. Reshape into (nxgrid, nygrid)


Example
-------

.. code-block:: python

    # To compute the density field
    dgrid = fiesta.sph.sph4grid2D(x, y, boxsize=100, ngrid=256)

    # To compute the field
    fgrid = fiesta.sph.sph4grid2D(x, y, boxsize=100, ngrid=256, f=f)

    # Compute the field and output the x and y grid
    fgrid, xg, yg = fiesta.sph.sph4grid2D(
        x, y,
        boxsize=100,
        ngrid=256,
        outputgrid=True
    )


sph4grid3D
==========

Same as ``sph4grid2D`` but extended to **3D SPH grid estimation**.

Algorithm
---------

1. Build KDTree using SPH3D
2. Assign masses
3. Create 3D grid
4. Generate subsampling offsets inside each voxel
5. For each offset:
   - shift grid positions
   - evaluate SPH density or field
6. Average over subsamples
7. Reshape into (nxgrid, nygrid, nzgrid)

Performance Notes
-----------------

- Computational cost scales as:
  - ``O(Ngrid × subsampling^d)``
- 3D is significantly more expensive than 2D
- Nearest-neighbour search dominates runtime
- Recommended:
  - keep ``k`` modest (20–50)
  - avoid large subsampling unless necessary

Example
-------

.. code-block:: python

    # To compute the density field
    dgrid = fiesta.sph.sph4grid3D(x, y, z, boxsize=100, ngrid=128)

    # To compute the field
    fgrid = fiesta.sph.sph4grid3D(x, y, z, boxsize=100, ngrid=128, f=f)

    # Compute the field and output the x and y grid
    fgrid, xg, yg, zg = fiesta.sph.sph4grid3D(
        x, y, z,
        boxsize=100,
        ngrid=128,
        outputgrid=True
    )

Notes
-----

- If ``mass`` is not provided, all particles are assumed to have unit mass.
- Periodic boundaries are handled internally in the KDTree construction.
- Subsampling improves accuracy by reducing pixel-center bias.

API Reference
=============

.. autofunction:: fiesta.sph.sph4grid2D
    :no-index:

.. autofunction:: fiesta.sph.sph4grid3D
    :no-index: