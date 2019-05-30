# Copyright 2016-2018 Iowa State University Research Foundation, Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from distutils.core import Extension
from distutils.command.build_ext import build_ext
import sys
import subprocess
import numpy as np
import os

#from Cython.Build import cythonize
#from numpy.distutils.core import setup as numpy_setup, Extension as numpy_Extension
from distutils.core import setup


setup(name="delamo",
      description="De-la-mo v2",
      author="Adarsh Krisnamurthy, Stephen D. Holland",
      url="http://thermal.cnde.iastate.edu/de-la-mo.xhtml",
      #ext_modules=[],
      package_dir={
            'autofiber': 'autofiber/autofiber'
      },
      packages=["delamo", "autofiber"],
      scripts=["scripts/delamo_process","scripts/delamo_test_processor"],
      #cmdclass = {"build_ext": build_ext_openmp},
)
