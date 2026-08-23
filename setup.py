#!/usr/bin/env python

from setuptools import setup

with open('README.rst') as fh:
    long_description = fh.read()

setup(name='molecular_hydrogen',
      version='0.2',
      description='Molecular Parameters for H2',
      long_description=long_description,
      long_description_content_type='text/x-rst',
      author='Adam Ginsburg',
      author_email='adam.g.ginsburg@gmail.com',
      url='https://github.com/keflavich/molecular_hydrogen/',
      packages=['molecular_hydrogen', 'molecular_hydrogen.tests'],
      package_data={'molecular_hydrogen': ['*.pyx']},
      install_requires=['numpy', 'astropy'],
      extras_require={'hitran': ['astroquery'],
                      'speed': ['cython'],
                      'test': ['pytest', 'astroquery']},
      license='MIT',
      classifiers=[
                   "Development Status :: 4 - Beta",
                   "Intended Audience :: Science/Research",
                   "Programming Language :: Python :: 3",
                   "Topic :: Scientific/Engineering :: Astronomy",
                   "License :: OSI Approved :: MIT License",
                  ],
     )
