SRC
===

Functions written with numba, used under-the-hood by other functions in the library.

Particle to Grid
----------------

.. autofunction:: fiesta.src.xgrid
.. autofunction:: fiesta.src.xgrids
.. autofunction:: fiesta.src.which_pix
.. autofunction:: fiesta.src.which_pixs
.. autofunction:: fiesta.src.pix1dto2d
.. autofunction:: fiesta.src.pix1dto3d
.. autofunction:: fiesta.src.periodic_pix
.. autofunction:: fiesta.src.ngp_pix
.. autofunction:: fiesta.src.cic_pix
.. autofunction:: fiesta.src.tsc_pix
.. autofunction:: fiesta.src.pcs_pix
.. autofunction:: fiesta.src.weight_cic
.. autofunction:: fiesta.src.weight_tsc
.. autofunction:: fiesta.src.weight_pcs
.. autofunction:: fiesta.src.part2grid_ngp_2d
.. autofunction:: fiesta.src.part2grid_cic_2d
.. autofunction:: fiesta.src.part2grid_tsc_2d
.. autofunction:: fiesta.src.part2grid_pcs_2d
.. autofunction:: fiesta.src.part2grid_ngp_3d
.. autofunction:: fiesta.src.part2grid_cic_3d
.. autofunction:: fiesta.src.part2grid_tsc_3d
.. autofunction:: fiesta.src.part2grid_pcs_3d

Bi/Trilinear Interpolation
--------------------------

.. autofunction:: fiesta.src.bilinear_periodic
.. autofunction:: fiesta.src.bilinear_nonperiodic
.. autofunction:: fiesta.src.bilinear_axisperiodic
.. autofunction:: fiesta.src.trilinear_periodic
.. autofunction:: fiesta.src.trilinear_nonperiodic
.. autofunction:: fiesta.src.trilinear_axisperiodic

Delaunay Functions
------------------

.. autofunction:: fiesta.src.triangle_area
.. autofunction:: fiesta.src.sum_triangle_area
.. autofunction:: fiesta.src.tetrahedron_volume
.. autofunction:: fiesta.src.voronoi_2d_area
.. autofunction:: fiesta.src.voronoi_3d_volume
.. autofunction:: fiesta.src.delaunay_area_2d
.. autofunction:: fiesta.src.sum_delaunay_area_4_points_2d
.. autofunction:: fiesta.src.get_delf0_2d
.. autofunction:: fiesta.src.delaunay_estimate_2d
.. autofunction:: fiesta.src.delaunay_volume_3d
.. autofunction:: fiesta.src.sum_delaunay_vol_4_points_3d
.. autofunction:: fiesta.src.get_delf0_3d
.. autofunction:: fiesta.src.delaunay_estimate_3d

Matrices
--------

.. autofunction:: fiesta.src.inv2by2
.. autofunction:: fiesta.src.inv3by3
.. autofunction:: fiesta.src.eig2by2
.. autofunction:: fiesta.src.symeig3by3

Differentiation
---------------

.. autofunction:: fiesta.src.dfdx_1d_periodic
.. autofunction:: fiesta.src.dfdx_2d_periodic
.. autofunction:: fiesta.src.dfdy_2d_periodic
.. autofunction:: fiesta.src.dfdx_3d_periodic
.. autofunction:: fiesta.src.dfdy_3d_periodic
.. autofunction:: fiesta.src.dfdz_3d_periodic

GridSPH
-------

.. autofunction:: fiesta.src.sum_from_integral_image_2D
.. autofunction:: fiesta.src.sum_from_integral_image_3D
.. autofunction:: fiesta.src.get_volume_enclosing_box_2D
.. autofunction:: fiesta.src.get_volume_enclosing_box_3D