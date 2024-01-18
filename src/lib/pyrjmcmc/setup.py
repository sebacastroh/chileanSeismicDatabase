# -*- coding: utf-8 -*-
"""
Created on Sat May 15 18:33:42 2021

@author: sebac
"""

try:
    from setuptools import setup
    from setuptools import Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension

from Cython.Distutils import build_ext

ext_modules = [Extension("utils",["utils.pyx"])]

setup(name= 'Generic model class',
      cmdclass = {'build_ext': build_ext},
      ext_modules = ext_modules)
