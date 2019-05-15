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



 
