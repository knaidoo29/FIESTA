===========
SPH Classes
===========

Overview
========

The SPH classes implement Smooth Particle Hydrodynamics (SPH) density estimation
and field interpolation using KDTree-based nearest-neighbour searches.

Supported dimensions:

- ``SPH1D`` — 1D particle systems
- ``SPH2D`` — 2D particle systems
- ``SPH3D`` — 3D particle systems

All classes support:

- Adaptive smoothing lengths based on nearest neighbours
- Density estimation
- Field interpolation
- Optional particle masses
- Cubic spline kernel

The general SPH estimator follows:

.. math::

    \rho(\mathbf{x}) = \sum_i m_i W(|\mathbf{x} - \mathbf{x}_i|, h)

    f(\mathbf{x}) = \frac{\sum_i m_i f_i W(|\mathbf{x} - \mathbf{x}_i|, h)}
                          {\sum_i m_i W(|\mathbf{x} - \mathbf{x}_i|, h)}

where ``W`` is the cubic SPH kernel and ``h`` is the smoothing length.

SPH1D
=====

The ``SPH1D`` class performs SPH estimation in one dimension.

Initialization
--------------

.. code-block:: python

    sph = fiesta.sph.SPH1D()

By default, the cubic kernel is used:

.. code-block:: python

    sph.kernel_type = "cubic"

Setup
-----

.. code-block:: python

    sph.setup(k=20, mass=None)

Parameters:

- ``k`` : int
    Number of nearest neighbours used for smoothing
- ``mass`` : array, optional
    Particle masses (defaults to 1 if None)

Assigning Mass
--------------

.. code-block:: python

    sph.assign_mass(mass)

If no mass is provided:

.. code-block:: python

    sph.mass = np.ones(len(points))

Density Estimation
------------------

.. code-block:: python

    dens = sph.get_density(x)

Returns SPH density evaluated at positions ``x``.

Field Interpolation
-------------------

.. code-block:: python

    sph.set_field(f)
    f_est = sph.estimate(x)

- ``f`` must be defined at particle positions
- Returns SPH-smoothed field values

Core Method
-----------

.. code-block:: python

    sph.sph_estimate(x, f=None, dens=None)

If ``f`` is None → returns density  
If ``f`` is provided → returns interpolated field

SPH2D
=====

The ``SPH2D`` class performs SPH estimation in two dimensions.

Initialization
--------------

.. code-block:: python

    sph = fiesta.sph.SPH2D()

Setup
-----

.. code-block:: python

    sph.setup(k=20, mass=None)

Arguments:

- ``k`` : number of nearest neighbours
- ``mass`` : optional particle masses

Kernel
------

The SPH kernel is evaluated in 2D:

.. code-block:: python

    sph.kernel(r, h)

Uses cubic spline kernel:

.. math::

    W(r, h) = W_{\mathrm{cubic}}^{2D}(r, h)

Density Estimation
------------------

.. code-block:: python

    dens = sph.get_density(x, y)

Field Interpolation
-------------------

.. code-block:: python

    sph.set_field(f)
    f_est = sph.estimate(x, y)

Core Method
-----------

.. code-block:: python

    sph.sph_estimate(x, y, f=None, dens=None)

Returns:

- density if ``f is None``
- interpolated field otherwise

SPH3D
=====

The ``SPH3D`` class performs SPH estimation in three dimensions.

Initialization
--------------

.. code-block:: python

    sph = fiesta.sph.SPH3D()

Setup
-----

.. code-block:: python

    sph.setup(k=50, mass=None)

Parameters:

- ``k`` : number of nearest neighbours (default: 50)
- ``mass`` : optional particle mass array

Kernel
------

.. code-block:: python

    sph.kernel(r, h)

Uses cubic kernel in 3D:

.. math::

    W(r, h) = W_{\mathrm{cubic}}^{3D}(r, h)

Density Estimation
------------------

.. code-block:: python

    dens = sph.get_density(x, y, z)

Field Interpolation
-------------------

.. code-block:: python

    sph.set_field(f)
    f_est = sph.estimate(x, y, z)

Core Method
-----------

.. code-block:: python

    sph.sph_estimate(x, y, z, f=None, dens=None)

Behavior:

- If ``f is None`` → returns density
- If ``f is provided`` → returns interpolated field

Algorithm Details
=================

Nearest Neighbour Search
------------------------

All SPH classes rely on KDTree queries:

.. code-block:: python

    nind, ndist = sph.nearest(..., k=self.k, return_dist=True)

- ``nind``: indices of nearest neighbours
- ``ndist``: distances to neighbours

Smoothing Length
----------------

For each query point:

.. code-block:: python

    h = max(ndist)

This defines the adaptive smoothing radius.

Kernel Weighting
----------------

Each contribution is weighted:

.. code-block:: python

    w_i = m_i * W(r_i, h)

Density or field is computed via summation.

Performance Notes
=================

- Complexity: O(N log N) for KDTree queries
- Best performance with vectorized mass arrays

Usage Example
-------------

.. code-block:: python

    sph = fiesta.sph.SPH2D()
    sph.setup(k=30)

    sph.points = np.column_stack([x, y])
    sph.assign_mass()

    dens = sph.get_density(x, y)

    sph.set_field(f)
    f_interp = sph.estimate(x, y)
