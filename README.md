# FIeld ESTimAtor

Estimates a field from values based on known input points. The module provides
field assignment and interpolation functions.

To-do:

* kDTree for fast nearest neighbour search.
* Dimensions:
  - 2D
  - 3D
  - Unit Sphere
* Assignment schemes:
  - Nearest neighbour
  - Grid based: NGP, CIC, TSC
  - Delaunay Tesselation
* Interpolation:
  - Delaunay Tesselation
  - Bilinear and Trilinear interpolation from a grid.
* Density estimation:
  - Grid based: NGP, CIC, TSC density estimation
  - Voronoi estimation
* Boundary management:
  - periodic conditions
  - buffer particles
  - mask handling for buffer particles.
* Serial implementation for data management.

## Python module

### Depencies

* numpy
* scipy
* f2py

### Functions
