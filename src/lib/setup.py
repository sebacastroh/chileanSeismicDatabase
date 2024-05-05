try:
    from setuptools import setup
    from setuptools import Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension

from Cython.Distutils import build_ext
import numpy as np

ext_modules = [
    Extension("seismic", ["seismic.pyx"], extra_compile_args=["-fopenmp"], extra_link_args=["-fopenmp"])
]

setup(
    name= 'Generic model class',
    cmdclass = {'build_ext': build_ext},
    include_dirs = [np.get_include()],
    ext_modules = ext_modules
)
