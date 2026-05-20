=============================
MPI Numerical Differentiation
=============================

The MPI differentiation functions provide distributed numerical derivatives
for data decomposed across MPI ranks. These functions are designed for
large datasets where the computational domain is split between processes.

The implementation supports:

* Distributed finite-difference differentiation
* Periodic and non-periodic boundaries
* 1D, 2D, and 3D scalar fields

Overview
========

The MPI differentiation module contains:

* ``mpi_dfdx`` — MPI derivative along the x-axis
* ``mpi_dfdy`` — MPI derivative along the y-axis
* ``mpi_dfdz`` — MPI derivative along the z-axis

The x-axis derivative performs explicit MPI communication because the data
domain is decomposed along the x direction. The y and z derivatives are
performed locally on each MPI slab.

MPI Domain Decomposition
------------------------

The data is assumed to be distributed across MPI ranks along the x-axis.

For example, a global 3D grid:

.. code-block:: text

    Global grid:
    +-----------------------------+
    | Rank 0 | Rank 1 | Rank 2    |
    +-----------------------------+

Each MPI rank owns a slab of the x-domain.

To compute derivatives near slab boundaries, neighboring values are exchanged
between MPI ranks using halo communication.

mpi_dfdx
========

Overview
--------

``mpi_dfdx`` computes the derivative along the x-axis for distributed data.

Because the x-axis is distributed across MPI processes, edge cells are
communicated between neighboring ranks before differentiation.

Periodic Boundaries
-------------------

For periodic domains:

* First and last MPI ranks communicate with each other
* Derivatives wrap continuously across the domain

For non-periodic domains:

* Edge ranks use one-sided finite differences

Example: 1D Distributed Derivative
----------------------------------

.. code-block:: python
    
    # Assuming correct initialisation of MPI object, see tutorial.
    import numpy as np
    import fiesta

    # Local grid
    x_local = np.linspace(0, 1, 32)

    # Local field
    f_local = np.sin(2 * np.pi * x_local)

    # Distributed derivative
    dfdx_local = fiesta.maths.mpi_dfdx(
        x_local,
        f_local,
        MPI,
        periodic=True
    )

Example: 3D Scalar Field
------------------------

.. code-block:: python

    # Assuming correct initialisation of MPI object, see tutorial.
    import numpy as np
    import fiesta

    nx, ny, nz = 32, 64, 64

    x_local = np.linspace(0, 1, nx)

    field = np.random.random((nx, ny, nz))

    dfdx = fiesta.maths.mpi_dfdx(
        x_local,
        field,
        MPI,
        periodic=True
    )

mpi_dfdy
========

Overview
--------

``mpi_dfdy`` computes derivatives along the y-axis.

Since the domain decomposition is along x, all y-values required for the
derivative already exist locally. No MPI communication is necessary.

Example
-------

.. code-block:: python

    # Assuming correct initialisation of MPI object, see tutorial.
    import fiesta

    dfdy = fiesta.maths.mpi_dfdy(
        ygrid,
        field,
        MPI,
        periodic=True
    )

mpi_dfdz
========

Overview
--------

``mpi_dfdz`` computes derivatives along the z-axis.

As with ``mpi_dfdy``, no MPI communication is required because the z-axis
is fully local to each rank.

Example
-------

.. code-block:: python

    # Assuming correct initialisation of MPI object, see tutorial.
    import fiesta

    dfdz = fiesta.maths.mpi_dfdz(
        zgrid,
        field,
        MPI,
        periodic=True
    )

Finite Difference Scheme
========================

The derivatives use:

* Symmetric two-point finite differences
* Three-point one-sided differences near boundaries

For interior points:

.. math::

    \frac{df}{dx} \approx
    \frac{f(x+\Delta x) - f(x-\Delta x)}
    {2\Delta x}

Boundary points use one-sided estimators.

Periodic Domains
================

When ``periodic=True``:

* Boundary values wrap around
* MPI edge ranks exchange data cyclically
* Derivatives remain continuous across domain edges

Example:

.. code-block:: python

    # Assuming correct initialisation of MPI object, see tutorial.
    import fiesta

    dfdx = fiesta.maths.mpi_dfdx(
        xgrid,
        field,
        MPI,
        periodic=True
    )

Non-Periodic Domains
====================

When ``periodic=False``:

* One-sided derivatives are applied near edges
* No wrapping occurs
* Physical domain boundaries are preserved

Example:

.. code-block:: python

    # Assuming correct initialisation of MPI object, see tutorial.
    import fiesta

    dfdx = fiesta.maths.mpi_dfdx(
        xgrid,
        field,
        MPI,
        periodic=False
    )

API References
==============

.. autofunction:: fiesta.maths.mpi_dfdx
    :no-index:

.. autofunction:: fiesta.maths.mpi_dfdy
    :no-index:

.. autofunction:: fiesta.maths.mpi_dfdz
    :no-index: