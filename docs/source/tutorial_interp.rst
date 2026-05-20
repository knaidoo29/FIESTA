==============
Interpolations
==============

Bilinear and Trilinear Interpolation
====================================

The ``bilinear`` and ``trilinear`` functions interpolate values from
regular Cartesian grids onto arbitrary coordinates.

These routines are useful for:

- Sampling scalar fields on regular grids.
- Interpolating simulation outputs.
- Evaluating density or fields at particle positions.
- Mapping gridded data back onto particles or trajectories.

The interpolation supports:

- Periodic boundaries
- Non-periodic boundaries
- Axis-dependent periodicity
- Custom coordinate origins
- Fill values for points outside the domain

Bilinear Interpolation
======================

The ``bilinear`` function interpolates values from a 2D grid.

Basic Example
--------------

Create a 2D grid and interpolate values at arbitrary coordinates.

.. code-block:: python

    import numpy as np
    import fiesta

    # Grid dimensions
    ngrid = 128
    boxsize = 10.0

    # Create coordinates
    xgrid = np.linspace(0, boxsize, ngrid)
    ygrid = np.linspace(0, boxsize, ngrid)

    X, Y = np.meshgrid(xgrid, ygrid, indexing="ij")

    # Example scalar field
    fgrid = np.sin(X) * np.cos(Y)

    # Query coordinates
    x = np.array([1.2, 3.7, 5.5])
    y = np.array([0.5, 4.2, 8.1])

    # Interpolate
    f = fiesta.interp.bilinear(fgrid, boxsize, x, y)

    print(f)

Periodic Interpolation
----------------------

For simulations with periodic boundaries:

.. code-block:: python

    f = fiesta.interp.bilinear(
        fgrid,
        boxsize,
        x,
        y,
        periodic=True,
    )

Coordinates wrapping around the domain edges are handled automatically.

Non-Periodic Interpolation
--------------------------

Disable periodic boundaries:

.. code-block:: python

    f = fiesta.interp.bilinear(
        fgrid,
        boxsize,
        x,
        y,
        periodic=False,
        fill_value=-1.0,
    )

Points outside the box receive ``fill_value``.

Axis-Dependent Periodicity
--------------------------

Use different boundary conditions along each axis.

Example:

- periodic along x
- non-periodic along y

.. code-block:: python

    f = fiesta.interp.bilinear(
        fgrid,
        boxsize=[10.0, 20.0],
        x=x,
        y=y,
        periodic=[True, False],
    )

Using a Custom Origin
---------------------

Shift the coordinate system origin.

.. code-block:: python

    f = fiesta.interp.bilinear(
        fgrid,
        boxsize=10.0,
        x=x,
        y=y,
        origin=-5.0,
    )

The valid domain becomes:

.. code-block:: text

    x ∈ [-5, 5]
    y ∈ [-5, 5]

Outside Boundary Handling
-------------------------

Coordinates outside the interpolation domain are assigned
``fill_value``.

Example:

.. code-block:: python

    x = np.array([-1.0, 5.0, 12.0])

    f = fiesta.interp.bilinear(
        fgrid,
        boxsize=10.0,
        x=x,
        y=y,
        periodic=False,
        fill_value=np.nan,
    )

Performance Notes
-----------------

The interpolation routines are optimized for large datasets.

For best performance:

- Use NumPy arrays.
- Avoid Python loops.
- Interpolate many coordinates simultaneously.

Example:

.. code-block:: python

    f = fiesta.interp.bilinear(fgrid, boxsize, x_array, y_array)

instead of:

.. code-block:: python

    for i in range(len(x_array)):
        fiesta.interp.bilinear(...)

Trilinear Interpolation
=======================

The ``trilinear`` function performs interpolation on 3D regular grids.

Basic Example
-------------

.. code-block:: python

    import numpy as np
    import fiesta

    ngrid = 64
    boxsize = 20.0

    xgrid = np.linspace(0, boxsize, ngrid)
    ygrid = np.linspace(0, boxsize, ngrid)
    zgrid = np.linspace(0, boxsize, ngrid)

    X, Y, Z = np.meshgrid(
        xgrid,
        ygrid,
        zgrid,
        indexing="ij",
    )

    # Example 3D field
    fgrid = np.sin(X) + np.cos(Y) + Z**2

    # Query points
    x = np.array([1.0, 5.0, 10.0])
    y = np.array([2.0, 6.0, 11.0])
    z = np.array([3.0, 7.0, 15.0])

    # Interpolate
    f = fiesta.interp.trilinear(fgrid, boxsize, x, y, z)

    print(f)

Periodic Boundaries
-------------------

Enable periodic interpolation:

.. code-block:: python

    f = fiesta.interp.trilinear(
        fgrid,
        boxsize,
        x,
        y,
        z,
        periodic=True,
    )

Non-Periodic Boundaries
-----------------------

.. code-block:: python

    f = fiesta.interp.trilinear(
        fgrid,
        boxsize,
        x,
        y,
        z,
        periodic=False,
        fill_value=0.0,
    )

Mixed Boundary Conditions
-------------------------

Specify periodicity independently per axis.

Example:

- periodic along x
- periodic along y
- non-periodic along z

.. code-block:: python

    f = fiesta.interp.trilinear(
        fgrid,
        boxsize=[10.0, 10.0, 50.0],
        x=x,
        y=y,
        z=z,
        periodic=[True, True, False],
    )

Using Custom Origins
--------------------

.. code-block:: python

    f = fiesta.interp.trilinear(
        fgrid,
        boxsize=10.0,
        x=x,
        y=y,
        z=z,
        origin=-5.0,
    )

Valid coordinate range becomes:

.. code-block:: text

    [-5, 5]

along all axes.

Implementation Details
----------------------

The interpolation internally:

1. Locates surrounding grid cells.
2. Computes interpolation weights.
3. Applies weighted averaging.
4. Handles periodic wrapping if enabled.
5. Assigns fill values outside the domain.

The routines automatically flatten the input grid for optimized
low-level evaluation.

API Reference
=============

.. autofunction:: fiesta.interp.bilinear
    :no-index:

.. autofunction:: fiesta.interp.trilinear
    :no-index: