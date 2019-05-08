import numpy as np

from delamo.api import DelamoModeler 
from delamo.api import Layer
from delamo.api import bond_layers
from delamo.api import SimpleCoordSys
from delamo import process
from delamo.layer import LayerMold

import os

# Front matter
# ------------


# Initialize the DeLaMo model
DM=DelamoModeler.Initialize(globals(),
                            pointtolerancefactor=100.0,
                            normaltolerance=100e-4)

# This script then generates both a CAD file and a Python script.
# The Python script can be run from Abaqus. It includes the 
# initialization script referenced above, and also opens 
# the CAD file and builds the model. 


# The name of the script file to generate and
# the name of the CAD file to write are returned
# by process.output_filenames()

# The first parameter to output_filenames
# should always match the name of the original script
# with the ".py" stripped

# In manually generated scripts, always specify phase
# to be "ORIGINAL"
(script_to_generate,
 cad_file_path_from_script,
 layer_boundary_template) = process.output_filenames("99_Test",phase="ORIGINAL")

# When writing a DeLaMo script, you start by creating a 
# finite element initialization script. This is a 
# Python script for ABAQUS that defines your various parameters
# -- material properties, etc. as Python variables.  
# In this case they are stored in the "abqparams_CFRP.py" file
DM.abaqus_init_script("abqparams_CFRP.py",globals())

# The above call automatically inserts wrapped copies of variables 
# defined in those scripts into the global variable space. Then you 
# can reference those variables in this script 

# (you can have as many init scripts as you like)


# The Delamo model contains generates several sets of instructions
# for different phases of the finite element modeling process:
# DM.initinstrs  (initialization)
# DM.assemblyinstrs (model assembly)
# DM.bcinstrs (boundary conditions)
# DM.meshinstrs (meshing)

# All methods called from those variables will go generally be executed
# in the assemblyinstrs pool unless otherwise overridden. You can
# use e.g. DM.meshinstrs.rewrapobj() to get a reference to 
# one of these variables that will execute in an alternate context. 
#
# For example, 
LaminateAssemblyMeshing=DM.meshinstrs.rewrapobj(LaminateAssembly)
# Creates a reference to the LaminateAssembly, for which method calls
# execute in the meshing context



# Basic parameters

# Set layer thickness we are planning on using
thickness = 0.199

# Load a NURBS mold surface from a file
Mold = LayerMold.FromFile(os.path.join("..","data","CurvedMold1.STEP"))

# Define a coordinate system
# This example defines +x direction along 0 deg. fibers,
# +y direction across 0 deg fibers, equivalent to
# the default (when coordsys is not specified)
coordsys=SimpleCoordSys((1.0,0.0,0.0),(0.0,1.0,0.0)) 


# Create 1st layer by moving the distance specified by thickness
# in the OFFSET_DIRECTION
layer1 = Layer.CreateFromMold(DM,Mold,"OFFSET",thickness,"Layer_1",LaminaSection,0,coordsys=coordsys)

# Once any breaks, etc. of a given layer are complete, it must be
# finalized. 
layer1.Finalize(DM)

# The MeshSimple method is a shortcut over the underlying ABAQUS routines
# It loops over each part in the layer and calls setElementType() with
# the specified MeshElemTypes, setMeshControls() with the given shape
# and technique, and seedPart() with the given mesh size, deviation factor,
# and minsizefactor. and refines the mesh near any given refined_edges

# Note that ABAQUS constants must be referenced as part of abqC
# rather than used directly 
layer1.MeshSimple(MeshElemTypes,meshsize,abqC.HEX_DOMINATED,abqC.SYSTEM_ASSIGN)

# Create and add point marker for fixed faced boundary condition
# There is a surface at y=-25 mm  from z= 0...0.2 mm
# This point identifies it
FixedPoint=[-40.0,-50.0,0.1] 

# Define a fixed boundary condition based on that point.
# EncastreBC is an ABAQUS function that was found by
# using the ABAQUS/CAE interface and then looking at the
# replay (.rpy) file. 
FEModel.EncastreBC(name="FixedFace_%d" % (DM.get_unique()),
                   createStepName=ApplyForceStep.name,
                   region=layer1.singlepart.GetInstanceFaceRegion(FixedPoint,0.07))


# Create 2nd layer
layer2 = Layer.CreateFromMold(DM,layer1.gk_layer.OffsetMold(),"OFFSET",thickness,"Layer_2", LaminaSection,-45,coordsys=coordsys)
#layer2.Split(os.path.join("..","data","SplitLine2.csv"),DM.abqpointtolerance)
#layer2.Split(os.path.join("..","data","nasa-delam12-1.csv"),DM.abqpointtolerance)

layer2.Finalize(DM)
layer2.MeshSimple(MeshElemTypes,meshsize/1.8,abqC.HEX_DOMINATED,abqC.SYSTEM_ASSIGN)

# Bond layers 1 and 2. With no other parameters, the layers are attached
# with a TIE boundary condition
#bond_layers(DM,layer1, layer2)
bond_layers(DM,layer1, layer2, defaultBC="COHESIVE",
            CohesiveInteraction=CohesiveInteraction,
            ContactInteraction=ContactInteraction,
            delaminationlist= [ "../data/nasa-delam12-2.csv" ])


# Create 3rd layer
#layer3 = Layer.CreateFromMold(DM,layer2.gk_layer.OffsetMold(),"OFFSET",thickness,"Layer_3",LaminaSection,45,coordsys=coordsys)

#layer3.Finalize(DM)
#layer3.MeshSimple(MeshElemTypes,meshsize/1.8,abqC.HEX_DOMINATED,abqC.SYSTEM_ASSIGN)

#bond_layers(DM,layer2, layer3)



# Finalization generates the output script and CAD model. 
DM.Finalize(script_to_generate,cad_file_path_from_script)
