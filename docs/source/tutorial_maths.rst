=====
Maths
=====

Numerical Differentiation
=========================

The differentiation routines compute finite-difference derivatives using:

- symmetric two-point derivatives
- non-symmetric three-point estimators near boundaries

These functions operate on regularly sampled grids.

dfdx
====

Computes the derivative along the x-axis.

Basic Example
-------------

Differentiate a sine function.

.. code-block:: python

    import numpy as np
    import fiesta

    xgrid = np.linspace(0, 2*np.pi, 100)

    f = np.sin(xgrid)

    dfdx_values = fiesta.maths.dfdx(xgrid, f)

Periodic Derivatives
--------------------

For periodic domains:

.. code-block:: python

    dfdx_values = fiesta.maths.dfdx(
        xgrid,
        f,
        periodic=True,
    )

This wraps the boundary values automatically.

Differentiating Multi-Dimensional Arrays
----------------------------------------

``dfdx`` differentiates along the first axis.

Example:

.. code-block:: python

    xgrid = np.linspace(0, 1, 128)

    f = np.random.random((128, 64))

    gx = fiesta.maths.dfdx(xgrid, f)

The derivative is computed along axis 0.

Boundary Handling
-----------------

When ``periodic=False``:

- symmetric derivatives are used internally
- edge cells use non-symmetric finite differences

When ``periodic=True``:

- ghost cells are added
- boundaries wrap around automatically

dfdy
====

Computes the derivative along the y-axis.

Basic Example
-------------

.. code-block:: python

    import numpy as np

    ygrid = np.linspace(0, 2*np.pi, 128)

    Y = np.meshgrid(ygrid, indexing="ij")[0]

    f = np.sin(Y)

    gy = fiesta.maths.dfdy(ygrid, f)

2D Field Example
----------------

.. code-block:: python

    nx = 128
    ny = 128

    x = np.linspace(0, 1, nx)
    y = np.linspace(0, 1, ny)

    X, Y = np.meshgrid(x, y, indexing="ij")

    f = np.sin(X) * np.cos(Y)

    gy = fiesta.maths.dfdy(y, f)

Periodic Boundaries
-------------------

.. code-block:: python

    gy = fiesta.maths.dfdy(
        ygrid,
        f,
        periodic=True,
    )

This is useful for periodic simulation domains.

dfdz
====

Computes the derivative along the z-axis.

Basic Example
-------------

.. code-block:: python

    import numpy as np

    zgrid = np.linspace(0, 2*np.pi, 64)

    f = np.sin(zgrid)

    gz = fiesta.maths.dfdz(zgrid, f)

3D Field Example
----------------

.. code-block:: python

    nx = ny = nz = 64

    x = np.linspace(0, 1, nx)
    y = np.linspace(0, 1, ny)
    z = np.linspace(0, 1, nz)

    X, Y, Z = np.meshgrid(
        x, y, z,
        indexing="ij"
    )

    f = X**2 + Y**2 + Z**2

    gz = fiesta.maths.dfdz(z, f)

Periodic 3D Simulations
-----------------------

.. code-block:: python

    gz = fiesta.maths.dfdz(
        zgrid,
        f,
        periodic=True,
    )

Finite Difference Notes
=======================

The derivative functions assume:

- Regularly spaced grids
- Uniform coordinate spacing

Grid spacing is inferred automatically:

.. code-block:: python

    dx = xgrid[1] - xgrid[0]

The implementation internally:

1. Extends arrays if periodic
2. Applies finite differences
3. Removes ghost zones
4. Returns derivatives on original grids

3x3 Determinant
===============

The ``det3by3`` function computes determinants for 3x3 matrices.

Single Matrix Example
---------------------

.. code-block:: python

    det = fiesta.maths.det3by3(
        1, 2, 3,
        0, 1, 4,
        5, 6, 0,
    )

Vectorized Matrix Evaluation
----------------------------

The routine supports arrays directly.

.. code-block:: python

    import numpy as np

    n = 1000

    m00 = np.random.random(n)
    m01 = np.random.random(n)
    m02 = np.random.random(n)

    m10 = np.random.random(n)
    m11 = np.random.random(n)
    m12 = np.random.random(n)

    m20 = np.random.random(n)
    m21 = np.random.random(n)
    m22 = np.random.random(n)

    det = fiesta.maths.det3by3(
        m00, m01, m02,
        m10, m11, m12,
        m20, m21, m22,
    )

This computes determinants for 1000 matrices simultaneously.

4x4 Determinant
===============

The ``det4by4`` function computes determinants for 4x4 matrices.

Basic Example
-------------

.. code-block:: python

    det = fiesta.maths.det4by4(
        1, 2, 3, 4,
        5, 6, 7, 8,
        2, 6, 4, 8,
        3, 1, 1, 2,
    )

Batch Matrix Determinants
-------------------------

The implementation is fully vectorized.

.. code-block:: python

    n = 10000

    matrices = np.random.random((n, 4, 4))

    det = fiesta.maths.det4by4(
        matrices[:,0,0],
        matrices[:,0,1],
        matrices[:,0,2],
        matrices[:,0,3],

        matrices[:,1,0],
        matrices[:,1,1],
        matrices[:,1,2],
        matrices[:,1,3],

        matrices[:,2,0],
        matrices[:,2,1],
        matrices[:,2,2],
        matrices[:,2,3],

        matrices[:,3,0],
        matrices[:,3,1],
        matrices[:,3,2],
        matrices[:,3,3],
    )

Implementation Details
----------------------

The determinant is computed via cofactor expansion:

.. math::

    \det(M) =
    m_{00} C_{00}
    - m_{01} C_{01}
    + m_{02} C_{02}
    - m_{03} C_{03}

where each cofactor uses the ``det3by3`` routine internally.

---------------------------------------------------------------------

Performance Notes
-----------------

The determinant routines:

- avoid Python loops
- operate directly on NumPy arrays
- support broadcasting
- efficiently evaluate large batches

API Reference
=============

Differentiation Functions
-------------------------

.. autofunction:: fiesta.maths.dfdx
    :no-index:

.. autofunction:: fiesta.maths.dfdy
    :no-index:

.. autofunction:: fiesta.maths.dfdz
    :no-index:

Determinant Functions
---------------------

.. autofunction:: fiesta.maths.det3by3
    :no-index:

.. autofunction:: fiesta.maths.det4by4
    :no-index:
