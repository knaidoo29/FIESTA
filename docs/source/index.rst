FIeld ESTimAtor
===============


+---------------+-----------------------------------------+
| Author        | Krishna Naidoo                          |
+---------------+-----------------------------------------+
| Version       | 0.0.1                                   |
+---------------+-----------------------------------------+
| Homepage      | https://github.com/knaidoo29/FIESTA     |
+---------------+-----------------------------------------+
| Documentation | https://fiesta-docs.readthedocs.io/     |
+---------------+-----------------------------------------+


Introduction
------------

FIESTA is a python library for general interpolation from uniform and non-uniform
input points. The library is predominantly written in python with a backend of
Fortran for speed.

MPI functionality can be enabled through the installation of the python library
mpi4py but will require the additional installation of MPIutils which handles
all of the MPI enabled functions. The class MPI is passed as an additional argument
for parallelisation.


Dependencies
------------

FIESTA is being developed on Python 3.8 but should work on all versions >3.4. Due
to the Fortran source code you will additionally require a fortran compiler usually
gfortran (which will come with gcc). The following Python modules are required.

* `numpy <http://www.numpy.org/>`_
* `scipy <https://scipy.org/>`_
* `SHIFT <https://github.com/knaidoo29/SHIFT>`_

For testing you will require `nose <https://nose.readthedocs.io/en/latest/>`_ or
`pytest <http://pytest.org/en/latest/>`_ .
