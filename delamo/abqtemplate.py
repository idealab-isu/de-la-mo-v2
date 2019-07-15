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


import sys
import re
import copy
import collections
import os
import math

###PATHLINE

import abaqus as abq
import abaqusConstants as abqC
import regionToolset
import mesh
import section
import assembly
import connector
#import connectorBehavior
#import job
#import load
import optimization
#import part
#import sketch
#import step
#import visualization
#import xyPlot

import numpy as np

if sys.version_info > (3,):
    long = int
    pass


# run initialization script (initintrs)

###RUN_INITINSTRS


# run model assembly script (assemblyinstrs)

###RUN_ASSEMBLYINSTRS

# run boundary condition script (bcinstrs)

###RUN_BCINSTRS

# run mesh generation script (meshinstrs)

###RUN_MESHINSTRS

# run fiber orientation script (fiberinstrs)

###RUN_FIBERINSTRS

# run model-execution script (runinstrs)

###RUN_RUNINSTRS



 
