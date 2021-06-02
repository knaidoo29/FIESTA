import setuptools
from numpy.distutils.core import setup, Extension

# read the contents of your README file
from os import path
this_directory = path.abspath(path.dirname(__file__))
print(this_directory[:-6])
with open(path.join(this_directory[:-6], 'README.md')) as f:
    long_description = f.read()

ext1 = Extension(name = 'fiesta.grid._part2grid', sources = ['fiesta/grid/_part2grid.f90'])
ext2 = Extension(name = 'fiesta.interp._bilinear', sources = ['fiesta/interp/_bilinear.f90'])
ext3 = Extension(name = 'fiesta.interp._trilinear', sources = ['fiesta/interp/_trilinear.f90'])

setup(name = 'fiesta',
      version = '0.0.0',
      description       = "FIeld ESTimAtor.",
      long_description  = long_description,
      long_description_content_type = 'text/markdown',
      url               = 'https://github.com/knaidoo29/fiesta',
      author            = "Krishna Naidoo",
      author_email      = "krishna.naidoo.11@ucl.ac.uk",
      license='MIT',
      packages=setuptools.find_packages(),
      install_requires=['numpy', 'scipy'],
      ext_modules = [ext1, ext2, ext3],
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
