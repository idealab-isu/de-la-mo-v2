import os
import os.path
import glob
import re
import numpy as np


# Set layer thickness for lamina
# *** MUST BE KEPT IN SYNC WITH 04_Delam_plate.py ***
thickness = 0.130


damage_directory = os.path.join("..","data","TRI_delams")
output_directory = "04_Delam_plate_output"

boundary_layers = glob.glob(os.path.join(output_directory,"layerboundary_PREDAMAGE_*.stl"))
if len(boundary_layers) != 7:
    raise ValueError("Did not find exactly seven PREDAMAGE .stl files in %s" % (output_directory))

linenumbers = [ int(re.match(r"""layerboundary_PREDAMAGE_([0-9]+)[.]stl""",os.path.split(boundary_layer_filename)[1]).group(1)) for boundary_layer_filename in boundary_layers ]
linenumbers.sort()

delamcounts = np.ones(len(linenumbers),dtype=np.uint32)

damage_filenames = glob.glob(os.path.join(damage_directory,"*.npz"))
for damage_filename in damage_filenames:
    damage_fh=np.load(damage_filename)
    delam_outline_raw = damage_fh["delam_coords"]

    delam_layernum = delam_outline_raw[0,0]
    assert((delam_outline_raw[0,:]==delam_layernum).all()) # All outline entries should specify the same layer number

    assert(delam_layernum==int(delam_layernum)) # should be an integer
    
    assert(delam_layernum > 0) # layer #1 is the first boundary
    
    assert((delam_outline_raw[3,:]==0.0)) # surface z coordinate assumed to be zero everywhere

    delam_depth = delam_layernum*thickness

    delam_coords = np.array((delam_outline_raw[1,:],delam_outline_raw[2,:],np.ones(delam_outline_raw.shape[1],dtype='d')*delam_depth),dtype='d')
    
    linenumber = linenumbers[int(delam_layernum)-1] # line number of bond_layers call

    
    csvfilename = os.path.join(output_directory,"layerdelam_%5.5d_%2.2d.csv" % (linenumber,delamcounts[int(delam_layernum)-1]))
    
    delamcounts[int(delam_layernum)-1] += 1

    np.savetxt(csvfilename,delam_coords.T,delimiter=", ")
    
    pass
