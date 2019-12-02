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

import numpy as np

from delamo.api import DelamoModeler
from delamo.api import Layer
from delamo.api import bond_layers
from delamo.api import SimpleCoordSys
from delamo.api import solid_solid_coupling
from delamo import process
from delamo.layer import LayerMold

import os

# Front matter
# ------------


# Initialize the DeLaMo model
DM = DelamoModeler.Initialize(globals(),
                              pointtolerancefactor=100.0,
                              normaltolerance=100e-4,
                              GapWidth=0.0)

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
 layer_boundary_template) = process.output_filenames("07_SolidSolidCoupling", process="DEFECT_INSERTION",
                                                     phase="ORIGINAL",
                                                     apply_damage_script="07_SolidSolidCoupling_add_damage.py", )

# When writing a DeLaMo script, you start by creating a 
# finite element initialization script. This is a 
# Python script for ABAQUS that defines your various parameters
# -- material properties, etc. as Python variables.  
# In this case they are stored in the "abqparams_CFRP.py" file
DM.abaqus_init_script("abqparams_CFRP.py", globals())

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
LaminateAssemblyMeshing = DM.meshinstrs.rewrapobj(LaminateAssembly)
# Creates a reference to the LaminateAssembly, for which method calls
# execute in the meshing context


# Basic parameters

# Set layer thickness for lamina
# *** MUST BE KEPT IN SYNC WITH 04_Delam_plate_add_damage.py ***
thickness1 = 2.19456 / 8.0
thickness2 = (4.57197 - 2.19456)/ 8.0

(OrigMold, SolidSolidCoupling) = solid_solid_coupling.from_solid_and_tool(DM,
                                                                          os.path.join("..", "data", "NASAShellOverwrap.STEP"),
                                                                          os.path.join("..", "data",
                                                                                       "CuttingTool2.STEP"),
                                                                          OrigDirPoint=np.array((0.0, 60.0, 0.0)),
                                                                          OrigDirNormal=np.array((0.0, 0.0, 1.0)))

# Define a coordinate system
# This example defines +x direction along 0 deg. fibers,
# +y direction across 0 deg fibers, equivalent to
# the default (when coordsys is not specified)
coordsys = SimpleCoordSys((1.0, 0.0, 0.0), (0.0, 1.0, 0.0))

layup = [0, 45, -45, 90, 90, -45, 45, 0, 0, 45, -45, 90, 90, -45, 45, 0]  # assumed layup

# Create and add point marker for fixed faced boundary condition
# There is a surface at y=-25 mm  from z= 0...0.2 mm
# This point identifies it
FixedPoint = [0, -107.95, 1.]
ForcePoint = [0, +107.95, 1.]

Mold = OrigMold
previouslayer = None
layers = []

# Create the flat region
for layernum in range(16):

    # Set the thickness for the 2 zones
    if (layernum < 8):
        thickness = thickness1
    else:
        thickness = thickness2

    layer = Layer.CreateFromMold(DM, Mold, "OFFSET", thickness, "Layer_%d" % (layernum + 1), LaminaSection,
                                 layup[layernum], coordsys=coordsys)


    # If it is the 9th layer, then cut the layer
    if (layernum == 8):
        layer.Split(os.path.join("..", "data", "SplitLineNASA.csv"), DM.abqpointtolerance)
        layer.gk_layer.RemoveLayerBody(1)
        pass

    layers.append(layer)

    # Once any breaks, etc. of a given layer are complete, it must be
    # finalized. 
    layer.Finalize(DM)

    # The MeshSimple method is a shortcut over the underlying ABAQUS routines
    # It loops over each part in the layer and calls setElementType() with
    # the specified MeshElemTypes, setMeshControls() with the given shape
    # and technique, and seedPart() with the given mesh size, deviation factor,
    # and minsizefactor. and refines the mesh near any given refined_edges

    # Note that ABAQUS constants must be referenced as part of abqC
    # rather than used directly 
    layer.MeshSimple(MeshElemTypes, meshsize, abqC.HEX_DOMINATED, abqC.SYSTEM_ASSIGN)

    # Bond layer to previous layer
    if previouslayer is not None:
        bond_layers(DM, previouslayer, layer, CohesiveInteraction=CohesiveInteraction,
                    ContactInteraction=ContactInteraction)
        pass

    # Embed layer in solid
    SolidSolidCoupling.embed_layer(DM, layer)

    # mold for next layer
    Mold = layer.gk_layer.OffsetMold()

    # previouslayer for next layer
    previouslayer = layer

    # loop back

    pass

# Define a fixed boundary condition based on FixedPoint.
# EncastreBC is an ABAQUS function that was found by
# using the ABAQUS/CAE interface and then looking at the
# replay (.rpy) file. 
# FEModel.EncastreBC(name="FixedFace_%d" % (DM.get_unique()),
#                    createStepName=ApplyForceStep.name,
#                    region=layer.singlepart.GetInstanceFaceRegion(FixedPoint, 0.02))
#
# ForceVector = [0.0, 0.0, -5e-2]  # Units of MPa
#
# # Call ABAQUS SurfaceTraction method
# # Again, this came from looking at ABAQUS replay (.rpy) output
# # Observe again that all ABAQUS symbolic constants need the "abqC"
# # prefix.
# FEModel.SurfaceTraction(name="SurfaceTraction_%d" % (DM.get_unique()),
#                         createStepName=ApplyForceStep.name,
#                         region=layers[0].singlepart.GetInstanceFaceRegionSurface(ForcePoint, 0.1),
#                         distributionType=abqC.UNIFORM,
#                         field='',
#                         localCsys=None,
#                         traction=abqC.GENERAL,
#                         follower=abqC.OFF,
#                         resultant=abqC.ON,
#                         magnitude=np.linalg.norm(ForceVector),
#                         directionVector=((0.0, 0.0, 0.0), tuple(ForceVector / np.linalg.norm(ForceVector))),
#                         amplitude=abqC.UNSET)

# You can have the job auto-start when the Python script is run
# DM.RunJob(BendingJob)

# Finalization generates the output script and CAD model. 
DM.Finalize(script_to_generate, cad_file_path_from_script)
