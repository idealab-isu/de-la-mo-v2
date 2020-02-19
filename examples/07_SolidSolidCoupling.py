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
# *** MUST BE KEPT IN SYNC WITH 07_SolidSolidCoupling_add_damage.py ***
numLayers = 8  # If updated must also update loop limit, below
thickness1 = 2.194565 / numLayers
thickness2 = (4.57197 - 2.194565)/ numLayers

(OrigMold, SolidSolidCoupling) = solid_solid_coupling.from_solid_and_tool(DM,
                                                                          os.path.join("..", "data", "NASAShellOverwrap.STEP"),
                                                                          os.path.join("..", "data",
                                                                                       "CuttingTool2.STEP"),
                                                                          OrigDirPoint=np.array((-200.0, 50.0, 0.0)),
                                                                          OrigDirNormal=np.array((0.0, 0.0, 1.0)))


MoldEdgePoints = OrigMold.GetPointsOnOuterEdges()




# Define a coordinate system
# This example defines +x direction along 0 deg. fibers,
# +y direction across 0 deg fibers, equivalent to
# the default (when coordsys is not specified)
coordsys = SimpleCoordSys((1.0, 0.0, 0.0), (0.0, 1.0, 0.0))

layup = [0, 45, -45, 90, 90, -45, 45, 0, 0, 45, -45, 90, 90, -45, 45, 0, 0, 45, -45, 90, 90, -45, 45, 0, 0, 45, -45, 90, 90, -45, 45, 0 ]  # assumed layup



# define properties for solid
solidmeshsize=15.0 # mm

# Define ABAQUS section for surrounding solid model.. Should really make principal axes function of position for stiffener portion...

## Unfortunately ABAQUS can't use a solid section for a complicated geometry
## that needs to be Tet meshed. But if it weren't so complicated and it
## could be meshed, here's how you would do it: 

#sectionlayers = [ section.SectionLayer(material='CFRPLaminaMat', thickness=0.125, orientAngle = layupangle, numIntPts=1,plyName='') for layupangle in layup ]

#SolidSection=FEModel.CompositeSolidSection(name='SolidSection',layupName='SolidSectionLayup',symmetric=False, layup = sectionlayers )

# Instead we model the solid as just a general solid with a single orientation
# For the moment, model the solid as uniaxial... Should probably build a hybrid stiffness model based on laminate theory
SolidSection=FEModel.HomogeneousSolidSection(name='LaminaSection',material=CFRPLaminaMat.name,thickness=None)

SolidSolidCoupling.solidpart.MeshSimple(MeshElemTypes, solidmeshsize, ElemShape=abqC.TET, ElemTechnique=abqC.FREE,refined_edges = MoldEdgePoints,pointtolerance=DM.abqpointtolerance,refinedmeshsize=meshsize)
SolidSolidCoupling.solidpart.AssignSection(SolidSection)
SolidSolidCoupling.solidpart.ApplyLayup(coordsys,0.0) # orientation of 0 means that 0 degrees as defined in the SolidSection layers lines up with the first axis (fiber direction) of the coordsys. 




# Create and add point marker for fixed faced boundary condition
# This point identifies it
FixedPoint = [0, +107.95, 1.]
ForcePoint = [0, -107.95, 1.]

Mold = OrigMold
previouslayer = None
layers = []

# Create the flat region
for layernum in range(16):  # Iteration limit should match numLayers*2 but must be a constant so that loop unwrapping can work

    # Set the thickness for the 2 zones
    if (layernum < numLayers):
        thickness = thickness1
        pass
    
    if layernum >= numLayers: # Avoid using else clause because that triggers a redbaron bug
        thickness = thickness2
        pass

    layer = Layer.CreateFromMold(DM, Mold, "OFFSET", thickness, "Layer_%d" % (layernum + 1), LaminaSection,
                                 layup[layernum], coordsys=coordsys)


    # If it is the 9th layer, then cut the layer
    if (layernum == numLayers):
        layer.Split(os.path.join("..", "data", "SplitLineNASA.csv"), DM.abqpointtolerance)
        layer.gk_layer.RemoveLayerBodyByPointInFace(np.array((-200.0, 60.0, 2.19456)), DM.abqpointtolerance)
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

    # Here layer #7 (0-based) has delaminations on both sides so it must be Tet-meshed
    if layernum==7:
        layer.MeshSimple(MeshElemTypes, meshsize, abqC.TET, abqC.SYSTEM_ASSIGN)
        pass
    else:
        layer.MeshSimple(MeshElemTypes, meshsize, abqC.HEX_DOMINATED, abqC.SYSTEM_ASSIGN)
        pass
    
    # Bond layer to previous layer
    if previouslayer is not None:
        bond_layers(DM, previouslayer, layer, CohesiveInteraction=CohesiveInteraction,ContactInteraction=ContactInteraction)
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
FEModel.EncastreBC(name="FixedFace_%d" % (DM.get_unique()),
                   createStepName=ApplyForceStep.name,
                   region=SolidSolidCoupling.solidpart.GetInstanceFaceRegion(FixedPoint, 0.02))


ForceVector = [0.0, 0.0, -5e-2]  # Units of MPa
#
# Call ABAQUS SurfaceTraction method
# Again, this came from looking at ABAQUS replay (.rpy) output
# Observe again that all ABAQUS symbolic constants need the "abqC"
# prefix.
FEModel.SurfaceTraction(name="SurfaceTraction_%d" % (DM.get_unique()),
                        createStepName=ApplyForceStep.name,
                        region=SolidSolidCoupling.solidpart.GetInstanceFaceRegionSurface(ForcePoint, 0.1),
                        distributionType=abqC.UNIFORM,
                        field='',
                        localCsys=None,
                        traction=abqC.GENERAL,
                        follower=abqC.OFF,
                        resultant=abqC.ON,
                        magnitude=np.linalg.norm(ForceVector),
                        directionVector=((0.0, 0.0, 0.0), tuple(ForceVector / np.linalg.norm(ForceVector))),
                        amplitude=abqC.UNSET)

# You can have the job auto-start when the Python script is run
# DM.RunJob(BendingJob)

# Finalization generates the output script and CAD model. 
DM.Finalize(script_to_generate, cad_file_path_from_script)
