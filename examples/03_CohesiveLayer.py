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
 layer_boundary_template) = process.output_filenames("03_CohesiveLayer",phase="ORIGINAL")

# When writing a DeLaMo script, you start by creating a 
# finite element initialization script. This is a 
# Python script for ABAQUS that defines your various parameters
# -- material properties, etc. as Python variables.  
# In this case they are stored in the "abqparams_CFRP.py" file
DM.abaqus_init_script("abqparams_CFRP.py",globals())
DM.abaqus_init_script("abqparams_CFRP_cohesive.py",globals())

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

# Set cohesive layer thickness
cohesivethickness = 0.001


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
layer2.Finalize(DM)
layer2.MeshSimple(MeshElemTypes,meshsize/1.8,abqC.HEX_DOMINATED,abqC.SYSTEM_ASSIGN)

# Bond layers 1 and 2. With no other parameters, the layers are attached
# with a TIE boundary condition
bond_layers(DM,layer1, layer2)

# Update and add point marker for fixed faced boundary condition
FixedPoint[2]+=thickness
FEModel.EncastreBC(name="FixedFace_%d" % (DM.get_unique()),
                   createStepName=ApplyForceStep.name,
                   region=layer2.singlepart.GetInstanceFaceRegion(FixedPoint,0.07))


# Create cohesive layer
layer23cohesive = Layer.CreateFromMold(DM,layer2.gk_layer.OffsetMold(),"OFFSET",cohesivethickness,"Layer23cohesive", CohesiveSection,0) # Orientation doesn't really matter since we have made the cohesive layer isotropic

#layer23cohesive.Finalize(DM)  # bond_layers() will do the finalize step on the cohesive layer so we don't have to

# Update and add point marker for fixed faced boundary condition
FixedPoint[2]+=cohesivethickness


# Create 3rd layer
layer3 = Layer.CreateFromMold(DM,layer23cohesive.gk_layer.OffsetMold(),"OFFSET",thickness,"Layer_3",LaminaSection,45,coordsys=coordsys)

layer3.Finalize(DM)
layer3.MeshSimple(MeshElemTypes,meshsize/1.8,abqC.HEX_DOMINATED,abqC.SYSTEM_ASSIGN)

# The line below performs a bonding operation with a delamination
# and a contact zone inside the delamination surrounded by a
# cohesive zone into which the delamination may grow
bond_layers(DM,layer2, layer3,
            cohesive_layer=layer23cohesive,
            defaultBC="COHESIVE_LAYER",
            delaminationlist=[os.path.join("..","data","nasa-delam12-1.csv")],
            ContactInteraction=ContactInteraction)

layer23cohesive.MeshCohesive(meshsize/1.8,abqC.HEX_DOMINATED)

            
# Update and add point marker for fixed faced boundary condition
FixedPoint[2]+=thickness
FEModel.EncastreBC(name="FixedFace_%d" % (DM.get_unique()),
                   createStepName=ApplyForceStep.name,
                   region=layer3.singlepart.GetInstanceFaceRegion(FixedPoint,0.07))

# Create 4th layer 
layer4 = Layer.CreateFromMold(DM,layer3.gk_layer.OffsetMold(),"OFFSET",thickness,"Layer_4",LaminaSection,90,coordsys=coordsys)
layer4.Finalize(DM)
layer4.MeshSimple(MeshElemTypes,meshsize/2.0,abqC.HEX_DOMINATED,abqC.SYSTEM_ASSIGN)

bond_layers(DM,layer3, layer4)

# Update and add point marker for fixed faced boundary condition
FixedPoint[2]+=thickness
FEModel.EncastreBC(name="FixedFace_%d" % (DM.get_unique()),
                   createStepName=ApplyForceStep.name,
                   region=layer4.singlepart.GetInstanceFaceRegion(FixedPoint,0.07))


# Create 5th layer over the layer 4 or the stiffener contour, if present
# ... for we just tell it to follow the layer 4 contour, which
# the stiffener automagically expanded
layer5 = Layer.CreateFromMold(DM,layer4.gk_layer.OffsetMold(),"OFFSET",thickness,"Layer_5",LaminaSection,90,coordsys=coordsys)
layer5.Finalize(DM)
layer5.MeshSimple(MeshElemTypes,meshsize/2.0,abqC.HEX_DOMINATED,abqC.SYSTEM_ASSIGN)

bond_layers(DM,layer4, layer5) 
    
# Update and add point marker for fixed faced boundary condition
FixedPoint[2]+=thickness
FixedPoint[1]-=.07 # accommodate outward shift as we go up
FEModel.EncastreBC(name="FixedFace_%d" % (DM.get_unique()),
                   createStepName=ApplyForceStep.name,
                   region=layer5.singlepart.GetInstanceFaceRegion(FixedPoint,0.07))


# Create 6th layer
layer6 = Layer.CreateFromMold(DM,layer5.gk_layer.OffsetMold(),"OFFSET",thickness,"Layer_6",LaminaSection,45,coordsys=coordsys)
layer6.Finalize(DM)
layer6.MeshSimple(MeshElemTypes,meshsize/2.0,abqC.HEX_DOMINATED,abqC.SYSTEM_ASSIGN)

bond_layers(DM,layer5, layer6)

# Update and add point marker for fixed faced boundary condition
FixedPoint[2]+=thickness
FEModel.EncastreBC(name="FixedFace_%d" % (DM.get_unique()),
                   createStepName=ApplyForceStep.name,
                   region=layer6.singlepart.GetInstanceFaceRegion(FixedPoint,0.07))


# Create 7th layer
layer7 = Layer.CreateFromMold(DM,layer6.gk_layer.OffsetMold(),"OFFSET",thickness,"Layer_7",LaminaSection,-45,coordsys=coordsys)
layer7.Finalize(DM)
layer7.MeshSimple(MeshElemTypes,meshsize/2.0,abqC.HEX_DOMINATED,abqC.SYSTEM_ASSIGN)

bond_layers(DM,layer6, layer7)

# Update and add point marker for fixed faced boundary condition
FixedPoint[2]+=thickness
FEModel.EncastreBC(name="FixedFace_%d" % (DM.get_unique()),
                   createStepName=ApplyForceStep.name,
                   region=layer7.singlepart.GetInstanceFaceRegion(FixedPoint,0.07))


# Create 8th layer
layer8 = Layer.CreateFromMold(DM,layer7.gk_layer.OffsetMold(),"OFFSET",thickness,"Layer_8",LaminaSection,0,coordsys=coordsys)
layer8.Finalize(DM)
layer8.MeshSimple(MeshElemTypes,meshsize/2.0,abqC.HEX_DOMINATED,abqC.SYSTEM_ASSIGN)

bond_layers(DM,layer7, layer8)

# Update and add point marker for fixed faced boundary condition
FixedPoint[2]+=thickness
FEModel.EncastreBC(name="FixedFace_%d" % (DM.get_unique()),
                   createStepName=ApplyForceStep.name,
                   region=layer8.singlepart.GetInstanceFaceRegion(FixedPoint,0.07))

# Can define a "Surface" that is visible in the Abaqus output database
# This is a direct ABAQUS call on the part object
# within layer1 (assumes layer1 is not split due to fiber/matrix breakage)
layer1.singlepart.fe_part.Surface(name="ForceSurface",
                                  side1Faces=layer1.singlepart.GetPartFace((-49.0,-49.0,thickness*0),0.1))


ForceVector=[ 0.0, 0.0, -5e-2 ] # Units of MPa 

# Call ABAQUS SurfaceTraction method
# Again, this came from looking at ABAQUS replay (.rpy) output
# Observe again that all ABAQUS symbolic constants need the "abqC"
# prefix. 
FEModel.SurfaceTraction(name="SurfaceTraction_%d" % (DM.get_unique()),
                        createStepName=ApplyForceStep.name,
                        region=layer1.singlepart.GetInstanceFaceRegionSurface((-49.0,-49.0,thickness*0.0),0.1),
                        distributionType=abqC.UNIFORM,
                        field='',
                        localCsys=None,
                        traction=abqC.GENERAL,
                        follower=abqC.OFF,
                        resultant=abqC.ON,
                        magnitude=np.linalg.norm(ForceVector),
                        directionVector=((0.0,0.0,0.0),tuple(ForceVector/np.linalg.norm(ForceVector))),
                        amplitude=abqC.UNSET)


# You can have the job auto-start when the Python script is run
#DM.RunJob(BendingJob)

# Finalization generates the output script and CAD model. 
DM.Finalize(script_to_generate,cad_file_path_from_script)
