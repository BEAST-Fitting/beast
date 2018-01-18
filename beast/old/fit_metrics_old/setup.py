#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy.distutils.misc_util


c_ext = [
    Extension("c_likelihood", ["likelihood.pyx"],
              extra_compile_args=['-fopenmp'],
              extra_link_args=['-fopenmp'],
              ),
    Extension("c_common", ["common.pyx"])
]


setup(
    name="proba_ext",
    cmdclass={"build_ext": build_ext},
    ext_modules=c_ext,
    include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs(),
)

# python setup.py build_ext --inplace
