from distutils.core import Extension
from distutils.command.build_ext import build_ext
import sys
import subprocess
import numpy as np
import os

#from Cython.Build import cythonize
#from numpy.distutils.core import setup as numpy_setup, Extension as numpy_Extension
from distutils.core import setup


setup(name="delamo2",
      description="De-la-mo v2",
      author="Adarsh Krisnamurthy, Stephen D. Holland",
      url="http://thermal.cnde.iastate.edu/de-la-mo.xhtml",
      #ext_modules=[],
      packages=["delamo2"],
      scripts=["scripts/delamo_process","scripts/delamo_test_processor"],
      cmdclass = {"build_ext": build_ext_openmp},
)
