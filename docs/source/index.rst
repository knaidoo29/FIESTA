FIESTA : FIeld ESTimAtor
========================

| | |
|---------|----------------|
| Author  | Krishna Naidoo |
| Version | 0.0.1          |
| Homepage | https://github.com/knaidoo29/FIESTA |
| Documentation | TBA |

FIESTA is a python library for general interpolation from uniform and non-uniform
input points. The library is predominantly written in python with a backend of Fortran
for speed.

MPI functionality can be enabled through the installation of the python library
mpi4py but will require the additional installation of MPIutils which handles
all of the MPI enabled functions. The class MPI is passed as an additional argument
for parallelisation.
