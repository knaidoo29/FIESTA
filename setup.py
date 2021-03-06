import setuptools
from numpy.distutils.core import setup, Extension

# read the contents of your README file
from os import path
this_directory = path.abspath(path.dirname(__file__))

with open(path.join(this_directory, 'README.md')) as f:
    long_description = f.read()

ext1 = Extension(name='fiesta.src.grid', sources=['fiesta/src/grid.f90'])
ext2 = Extension(name='fiesta.src.part2grid_pix', sources=['fiesta/src/part2grid_pix.f90'])
ext3 = Extension(name='fiesta.src.part2grid_wei', sources=['fiesta/src/part2grid_wei.f90'])
ext4 = Extension(name='fiesta.src.part2grid_2d', sources=['fiesta/src/part2grid_2d.f90'])
ext5 = Extension(name='fiesta.src.part2grid_3d', sources=['fiesta/src/part2grid_3d.f90'])
ext6 = Extension(name='fiesta.src.bilinear', sources=['fiesta/src/bilinear.f90'])
ext7 = Extension(name='fiesta.src.trilinear', sources=['fiesta/src/trilinear.f90'])
ext8 = Extension(name='fiesta.src.polygon', sources=['fiesta/src/polygon.f90'])
ext9 = Extension(name='fiesta.src.voronoi2d', sources=['fiesta/src/voronoi2d.f90'])
ext10 = Extension(name='fiesta.src.polyhedron', sources=['fiesta/src/polyhedron.f90'])
ext11 = Extension(name='fiesta.src.voronoi3d', sources=['fiesta/src/voronoi3d.f90'])
ext12 = Extension(name='fiesta.src.matrix', sources=['fiesta/src/matrix.f90'])
ext13 = Extension(name='fiesta.src.delaunay2d', sources=['fiesta/src/delaunay2d.f90'])
ext14 = Extension(name='fiesta.src.delaunay3d', sources=['fiesta/src/delaunay3d.f90'])
ext15 = Extension(name='fiesta.src.differentiate', sources=['fiesta/src/differentiate.f90'])

exts = [ext1, ext2, ext3, ext4, ext5, ext6, ext7, ext8, ext9, ext10,
        ext11, ext12, ext13, ext14, ext15]

setup(name = 'fiesta',
      version = '0.0.1',
      description       = "FIeld ESTimAtor.",
      long_description  = long_description,
      long_description_content_type = 'text/markdown',
      url               = 'https://github.com/knaidoo29/fiesta',
      author            = "Krishna Naidoo",
      author_email      = "krishna.naidoo.11@ucl.ac.uk",
      license='MIT',
      packages=setuptools.find_packages(),
      install_requires=['numpy', 'scipy'],
      ext_modules = exts,
      python_requires = '>=3',
      classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Scientific/Engineering :: Simulation',
      ],
      test_suite='nose.collector',
      tests_require=['nose'],
      )
