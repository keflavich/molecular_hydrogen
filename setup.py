#!/usr/bin/env python

import os
import shutil
import sys
from distutils.core import setup, Command

if 'develop' in sys.argv:
    # use setuptools for develop, but nothing else
    from setuptools import setup

setup(name='molecular_hydrogen',
      version='0.1',
      description='Molecular Parameters for H2',
      author=['Adam Ginsburg'],
      author_email=['adam.g.ginsburg@gmail.com'], 
      url='https://github.com/keflavich/molecular_hydrogen/',
      packages=['molecular_hydrogen'],
      requires=['numpy'],
      classifiers=[
                   "Development Status :: 4 - Beta",
                   "Programming Language :: Python",
                   "License :: OSI Approved :: MIT License",
                  ],
      
     )

