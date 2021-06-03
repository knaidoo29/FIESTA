# FIESTA : FIeld ESTimAtor

Author:         Krishna Naidoo                    
Version:        0.0.0                               
Homepage:       https://github.com/knaidoo29/FIESTA     
Documentation:  tba                                  


Routines for Estimating a field from values based on known input points.

## Development

FIESTA is currently in development. The plan is to implement the bulk of the code in python with Fortran subroutines for quick computations on large datasets. At some point much of this will have to be ported to C/C++ for very large datasets and the efficient use of MPI.

### TO DO

* Dimensions:
  * Unit Sphere
  * Polar Coordinates
* Assignment schemes:
  * Nearest neighbour
* Interpolation:
  * Delaunay Tesselation
* Density estimates:
  * Grid
  * Voronoi
* Field estimator general
  * Nearest neighbour
  * Grid
  * Delaunay
* Boundary management:
 * Buffer particles outside mask
* Serial implementation for very large datasets

### Implemented in Python

* kDTree for fast nearest neighbour search.
* Dimensions:
  * 2D
  * 3D
* Assignment schemes:
  * Grid based: NGP, CIC, TSC
* Interpolation:
  * Bilinear
  * Trilinear
* Boundary particles:
  * Random buffer
  * Periodic

## Python module

### Installation

Clone or download FIESTA. Navigate to the python directory in a terminal and install using the following:

```
    python setup.py build
    python setup.py install
```

You will then be able to load FIESTA in python using:

```
  import fiesta
```

### Depencies

* python3
* numpy
* scipy
* f2py

### Functions

* `boundary` : Boundary particles.
  * `boundary.buffer_random_particles_2d` : Random buffer particles in 2D.
  * `boundary.buffer_random_particles_3d` : Random buffer particles in 3D.
  * `boundary.buffer_periodic_particles_2d` : Periodic particles in 2D.
  * `boundary.buffer_periodic_particles_3d` : Periodic particles in 3D.

* `grid` : Grid based functions.
  * `grid.grid2d` : Generates a 2D grid.
  * `grid.grid3d` : Generates a 3D grid.
  * `grid.part2grid2d` : Particle to grid assignment in 2D.
  * `grid.part2grid3d` : Particle to grid assignment in 3D.

* `interp` : Interpolation functions.
  * `interp.bilinear` : Bilinear interpolation from a grid.
  * `interp.trilinear` : Trilinear interpolation from a grid.

* `kdtree` : kDTree class built on `scipy`.
  * `kdtree.KDTree2D` : Class implemented in 2D.
  * `kdtree.KDTree3D` : Class implemented in 3D.

* `randoms` : Random generation.
  * `randoms.random_uniform` : Uniform randoms in 1D.
  * `randoms.random_box` : Uniform randoms in 2D.
  * `randoms.random_cube` : Uniform randoms in 3D.

* `src` : Direct access to lower level Fortran source code.
  * `src.xgrid` : Computes grid point based on the grid index and grid size.
  * `src.weight_cic` : Cloud-in-cell weight.
  * `src.weight_tsc` : Triangular-shaped-cloud weight.
  * `src.p2g_ngp_2d` : Nearest-grid-point particle assignment in 2D.
  * `src.p2g_cic_2d_periodic` : Cloud-in-cell particle assignment in 2D with periodic boundaries.
  * `src.p2g_cic_2d_nonperiodic` : Cloud-in-cell particle assignment in 2D with non-periodic boundaries.
  * `src.p2g_tsc_2d_periodic` : Triangular-shaped-cloud particle assignment in 2D with periodic boundaries.
  * `src.p2g_tsc_2d_nonperiodic` : Triangular-shaped-cloud particle assignment in 2D with non-periodic boundaries.
  * `src.p2g_ngp_3d` : Nearest-grid-point particle assignment in 3D.
  * `src.p2g_cic_3d_periodic` : Cloud-in-cell particle assignment in 3D with periodic boundaries.
  * `src.p2g_cic_3d_nonperiodic` : Cloud-in-cell particle assignment in 3D with non-periodic boundaries.
  * `src.p2g_tsc_3d_periodic` : Triangular-shaped-cloud particle assignment in 3D with periodic boundaries.
  * `src.p2g_tsc_3d_nonperiodic` : Triangular-shaped-cloud particle assignment in 3D with non-periodic boundaries.
  * `src.bilinear_periodic` : Bilinear interpolation from a grid with periodic boundaries.
  * `src.bilinear_nonperiodic` : Bilinear interpolation from a grid with non-periodic boundaries.
  * `src.trilinear_periodic` : Trilinear interpolation from a grid with periodic boundaries.
  * `src.trilinear_nonperiodic` : Trilinear interpolation from a grid with non-periodic boundaries.
