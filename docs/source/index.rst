.. image:: _static/FiEsta_logo_large_white.jpg
   :align: center
   :class: only-light

.. image:: _static/FiEsta_logo_large_black.jpg
   :align: center
   :class: only-dark


FIeld ESTimAtor
===============


+---------------+-----------------------------------------+
| Author        | Krishna Naidoo                          |
+---------------+-----------------------------------------+
| Version       | 0.2.0                                   |
+---------------+-----------------------------------------+
| Repository    | https://github.com/knaidoo29/FIESTA     |
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

* `numba <https://numba.pydata.org/>`_
* `numpy <http://www.numpy.org/>`_
* `scipy <https://scipy.org/>`_
* `SHIFT <https://github.com/knaidoo29/SHIFT>`_
* `mpi4py <https://mpi4py.readthedocs.io/en/stable/>`_ [Optional: enables MPI parallelism]

Installation
------------

FIESTA can be installed by cloning the github repository::

    git clone https://github.com/knaidoo29/FIESTA.git
    cd FIESTA
    python setup.py build
    python setup.py install

Once this is done you should be able to call FIESTA from python:

.. code-block:: python

    import fiesta

Support
-------

If you have any issues with the code or want to suggest ways to improve it please
open a new issue (`here <https://github.com/knaidoo29/FIESTA/issues>`_)
or (if you don't have a github account) email krishna.naidoo.11@ucl.ac.uk.

Contents
--------

Version History
---------------

**Version 0.0**:

  * Interpolation and field assignments for cartesian coordinates up to 3 dimensions.

  * Grid based methods:

    - Field assignments: Nearest Grid Point, Cloud In Cell and Triangular Shaped Cloud.

    - Interpolation: bilinear and trilinear interpolation.

  * Non-uniform based methods:

    - Field assignments and interpolation estimation via voronoi tesselation, delaunay tesselation and smooth particle hydrodynamics.
