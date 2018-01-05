#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 26 12:56:47 2017

@author: marcel
"""

from setuptools import setup

README = 'README.md'
VERSION = 0.11
DESCRIPTION = 'Phantom Toolkit for analyzing SimpleITK images'
NAME = 'SimplePhantomToolkit'




def readme():
    with open(README) as f:
        return f.read()


setup(name=NAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description=readme(),
      classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering :: Medical Science Apps.',
        'Intended Audience :: Science/Research',
        'Natural Language :: English'
      ],
      keywords='phantom analysis simpleitk medical dicom',
      url='https://github.com/heydude1337/SimplePhantomToolkit',
      author='HeyDude',
      author_email='heydude1337@gmail.com',
      license='MIT',
      packages=['SimplePhantomToolkit', 'SimplePhantomToolkit.phantoms'],
      install_requires=[
          'SimpleITK', 'matplotlib', 'pandas', 'numpy', 'scipy', 'seaborn'
      ],
      include_package_data=True,
      zip_safe=False)