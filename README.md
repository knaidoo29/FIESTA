![biglogo](
https://raw.githubusercontent.com/knaidoo29/FIESTA/main/docs/source/_static/FiEsta_logo_large_github.jpg)

# FIESTA: Field Interpolation and Estimation using Spatial Techniques and Algorithms

<p align="center">
    <a href="https://github.com/knaidoo29/FIESTA/actions/workflows/tests.yml">
    <img src="https://github.com/knaidoo29/FIESTA/actions/workflows/tests.yml/badge.svg" alt="Python Tests">
    </a>
    <a href="https://codecov.io/gh/knaidoo29/FIESTA" > 
    <img src="https://codecov.io/gh/knaidoo29/FIESTA/graph/badge.svg?token=RFSEAKTTKG"/> 
    </a>
    <a href="https://img.shields.io/badge/Python-3.8%20|%203.9%20|%203.10%20|%203.11%20|%203.12-blue">
    <img src="https://img.shields.io/badge/Python-3.8%20|%203.9%20|%203.10%20|%203.11%20|%203.12-blue" alt="Python Version Support">
    </a>
    <a href="https://img.shields.io/github/v/release/knaidoo29/FIESTA">  
    <img src="https://img.shields.io/github/v/release/knaidoo29/FIESTA" alt="Version">
    </a>
    <a href="https://pypi.org/project/fiesta-pkg/">
    <img src="https://img.shields.io/pypi/v/fiesta-pkg.svg" alt="PyPI version">
    </a>
    <a href="https://fiesta-docs.readthedocs.io/">
    <img src="https://readthedocs.org/projects/fiesta-docs/badge/?version=latest" alt="Documentation Status">
    </a>
    <a href="https://github.com/knaidoo29/FIESTA">
    <img src="https://img.shields.io/badge/GitHub-repo-blue?logo=github" alt="GitHub repository">
    </a>
    <a href="https://img.shields.io/github/stars/knaidoo29/fiesta">
    <img src="https://img.shields.io/github/stars/knaidoo29/fiesta" alt="github: stars">
    </a>
    <a href="https://img.shields.io/github/stars/knaidoo29/fiesta">
    <img src="https://img.shields.io/github/forks/knaidoo29/fiesta" alt="github: forks">
    </a>
    <a href="https://opensource.org/licenses/MIT">
    <img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT">
    </a>
    <a href="https://github.com/psf/black">
    <img src="https://img.shields.io/badge/code%20style-black-000000.svg" alt="Code style: black">
    </a>
    <!-- <a href="https://doi.org/10.5281/zenodo.17093446">
    <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.17093446.svg" alt="zenodo: DOI">
    </a> -->
</p>

## Introduction

``FIESTA`` is a python library for general interpolation from uniform and non-uniform input points. The library is written in python with numba acceleration for speed. The library has the *optional* capability to be used on large data-sets by distributing jobs across multiple processes via ``MPI``. This relies on the ``mpi4py`` library and the ``MPI`` object from the ``shift-fft`` package (see [here](https://shift.readthedocs.io/en/latest/mpiutils.html)) which is passed as an additional object into ``MPI`` related functions.


## Dependencies

* [numba](https://numba.pydata.org/)
* [numpy](http://www.numpy.org/)
* [scipy](http://scipy.org/)
* [shift-fft](https://shift.readthedocs.io/en/latest/index.html)
* [mpi4py](https://mpi4py.readthedocs.io/en/stable/) [Optional: enables MPI distributed processes]

## Installation

FIESTA can be installed via ``pip``:

```
  pip install fiesta-pkg
```

or by cloning the repository:

```
    git clone https://github.com/knaidoo29/fiesta.git
    cd FIESTA
    pip install .
```

Once this is done you should be able to call `fiesta` from python:

```
    import fiesta
```

To use the ``MPI`` functionality please take a look at the documentation in FIESTA which instructs users how to use the ``shift-fft`` ``MPI`` object and how to run these distributed jobs successfully without errors or MPI related hanging.

## Documentation

In depth documentation and tutorials are provided [here](https://fiesta-docs.readthedocs.io/).

## Citation

If you use **FIESTA** in your work, please cite:

<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17093446.svg)](https://doi.org/10.5281/zenodo.17093446)

```bibtex
  @software{naidoo_shift_2025,
    author       = {Naidoo, Krishna},
    title        = {SHIFT: a scalable MPI library for computing fast Fourier transforms in python},
    year         = 2025,
    publisher    = {Zenodo},
    doi          = {10.5281/zenodo.17093446},
    url          = {https://doi.org/10.5281/zenodo.17093446}}
``` -->

## Support

If you have any issues with the code or want to suggest ways to improve it please open a new issue ([here](https://github.com/knaidoo29/FIESTA/issues)) or (if you don't have a github account) email _krishna.naidoo.11@ucl.ac.uk_.
