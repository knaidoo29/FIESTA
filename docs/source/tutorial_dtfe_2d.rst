DTFE in 2D
==========

Introduction
------------

FIESTA provides Delaunay Tessellation Field Estimation (DTFE) tools for
interpolating irregularly sampled data onto continuous fields using
Delaunay triangulations.

The DTFE implementation supports:

- density estimation
- scalar field interpolation
- boundary detection
- simplex querying
- distributed tessellation merging

Main classes:

- ``Delaunay2D``
- ``DelaunayMerger2D``

Basic Workflow
--------------

The typical DTFE workflow is:

1. create the tessellation
2. compute densities or assign fields
3. estimate values at arbitrary coordinates

Example workflow:

.. code-block:: python

   dtfe = fiesta.dtfe.Delaunay2D()

   dtfe.setup(x, y, boundary)

   dtfe.run()

   dtfe.get_dens()

   dtfe.set_field2dens()

   f_est = dtfe.estimate(xq, yq)

Creating a Delaunay Tessellation
--------------------------------

Generate random points:

.. code-block:: python

   import numpy as np
   import fiesta

   np.random.seed(0)

   npoints = 1000

   x = np.random.uniform(0.0, 1.0, npoints)
   y = np.random.uniform(0.0, 1.0, npoints)

Define the domain boundary:

.. code-block:: python

   boundary = [0.0, 1.0, 0.0, 1.0]

Construct the tessellation:

.. code-block:: python

   dtfe = fiesta.dtfe.Delaunay2D()

   dtfe.setup(x, y, boundary)

   dtfe.run()

The tessellation is now available through:

.. code-block:: python

   dtfe.simplices

Density Estimation
------------------

DTFE naturally estimates densities from the inverse area surrounding each
particle.

Compute point densities:

.. code-block:: python

   dtfe.get_dens()

The densities are stored in:

.. code-block:: python

   dtfe.point_dens

Convert the density into the interpolated field:

.. code-block:: python

   dtfe.set_field2dens()

Field Interpolation
-------------------

Instead of density estimation, arbitrary scalar fields can be assigned to
the particles.

Example:

.. code-block:: python

   f = np.sin(2*np.pi*x) * np.cos(2*np.pi*y)

Assign the field:

.. code-block:: python

   dtfe.set_field(f)

The field can now be interpolated continuously.

Estimating the Field
--------------------

Evaluate the field at arbitrary coordinates.

Generate query points:

.. code-block:: python

   xq = np.random.uniform(0.0, 1.0, 1000)
   yq = np.random.uniform(0.0, 1.0, 1000)

Estimate the field:

.. code-block:: python

   f_est = dtfe.estimate(xq, yq)

The returned values contain interpolated estimates inside valid simplices.

Points outside the tessellation return:

.. code-block:: python

   np.nan

Finding Simplices
-----------------

Determine which simplex contains a coordinate:

.. code-block:: python

   simplex_ids = dtfe.find_simplex(xq, yq)

Returned values:

- ``-1`` indicates the point lies outside the tessellation
- non-negative integers identify simplex indices

Boundary Handling
-----------------

DTFE interpolation near boundaries can become unreliable because simplices
may extend beyond the sampled domain.

FIESTA automatically classifies:

- interior simplices
- boundary simplices
- exterior simplices

The simplex classifications are stored in:

.. code-block:: python

   dtfe.simptypes

Point classifications are stored in:

.. code-block:: python

   dtfe.point_type

Classification meanings:

+--------+--------------------------------+
| Value  | Meaning                        |
+========+================================+
| 1      | interior                       |
+--------+--------------------------------+
| 0      | boundary                       |
+--------+--------------------------------+
| -1     | invalid / exterior             |
+--------+--------------------------------+

Allowing Border Interpolation
-----------------------------

By default, interpolation excludes boundary simplices.

To allow interpolation near borders:

.. code-block:: python

   f_est = dtfe.estimate(
       xq,
       yq,
       allow_border=True
   )

This option is only valid for non-density field interpolation.

Circumcircles
-------------

The circumcircle of each simplex can be computed:

.. code-block:: python

   dtfe.get_circumcircles()

Results are stored in:

.. code-block:: python

   dtfe.circumcircles

Each row contains:

.. code-block:: text

   [x_center, y_center, radius]

Triangle Areas
--------------

Compute simplex areas:

.. code-block:: python

   dtfe.get_area()

Results:

.. code-block:: python

   dtfe.delaunay_area

Border Extraction
-----------------

Extract boundary information from a tessellation:

.. code-block:: python

   border = dtfe.get_border()

The returned dictionary contains:

- border particles
- simplex information
- field values
- point classifications

This is useful for distributed or tiled DTFE calculations.

Distributed DTFE Merging
------------------------

``DelaunayMerger2D`` merges border regions from multiple DTFE calculations.

This is intended for:

- MPI workflows
- tiled datasets
- domain decomposition methods

Creating a Merger
-----------------

.. code-block:: python

   merger = fiesta.dtfe.DelaunayMerger2D()

Adding Borders
--------------

Suppose multiple tessellations were computed independently:

.. code-block:: python

   border1 = dtfe1.get_border()

   border2 = dtfe2.get_border()

Add them to the merger:

.. code-block:: python

   merger.add_border(border1)

   merger.add_border(border2)

Construct the merged tessellation:

.. code-block:: python

   merger.run(boundary)

Estimating from the Merged DTFE
-------------------------------

Interpolate using the merged tessellation:

.. code-block:: python

   f_est = merger.estimate(xq, yq)

Filtering Boundary Regions
--------------------------

The merger can optionally remove simplices outside the desired domain:

.. code-block:: python

   merger.run(
       boundary,
       apply_filter=True
   )

Cleaning Objects
----------------

Reset a DTFE object:

.. code-block:: python

   dtfe.clean()

Reset a merger object:

.. code-block:: python

   merger.clean()

Performance Notes
-----------------

DTFE computations involve:

- Delaunay triangulation construction
- simplex neighbour analysis
- geometric interpolation

For large datasets:

- triangulation can become memory intensive
- interpolation scales approximately linearly after triangulation
- boundary merging adds overhead

The implementation is best suited for:

- irregular spatial data
- particle simulations
- adaptive density estimation
- mesh-free interpolation

API Reference
-------------

.. autoclass:: fiesta.dtfe.Delaunay2D
   :members:
   :no-index:

.. autoclass:: fiesta.dtfe.DelaunayMerger2D
   :members:
   :no-index: