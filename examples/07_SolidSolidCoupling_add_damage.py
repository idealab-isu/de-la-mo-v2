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

import copy
import os
import os.path
import glob
import re
import numpy as np
import scipy 
import scipy.spatial

import pandas as pd

doplots = False
numLayers = 8

# Set layer thickness for lamina
# *** MUST BE KEPT IN SYNC WITH 07_SolidSolidCoupling.py ***
thickness1 = 2.194565 / numLayers
thickness2 = (4.57197 - 2.194565)/ numLayers

damage_directory = os.path.join("..","data","NASA_delams")
output_directory = "07_SolidSolidCoupling_output"

boundary_layers = glob.glob(os.path.join(output_directory,"layerboundary_PREDAMAGE_*.stl"))
if len(boundary_layers) != 15:
    raise ValueError("Did not find exactly fifteen PREDAMAGE .stl files in %s" % (output_directory))

linenumbers = [ int(re.match(r"""layerboundary_PREDAMAGE_([0-9]+)[.]stl""",os.path.split(boundary_layer_filename)[1]).group(1)) for boundary_layer_filename in boundary_layers ]
linenumbers.sort()

delamcounts = np.ones(len(linenumbers),dtype=np.uint32)

if doplots:
    from matplotlib import pyplot as pl
    pass

legends={}

pixel_size = 0.02*25.4  # pixel size (mm) corresponding to .02 inches

damage_filenames = glob.glob(os.path.join(damage_directory,"*.xlsx"))
for damage_filename in damage_filenames:
    damage_dataframe = pd.read_excel(damage_filename,header=1,index_col=1)
    damage_dataframe.drop(columns=damage_dataframe.keys()[0],inplace=True)
    
    # Remove first column (with label y and nothing else)
    
    damage_matrix = damage_dataframe.values
    x_scan_coordinates = np.array([float(colindex) for colindex in damage_dataframe.keys()])*pixel_size
    y_scan_coordinates = damage_dataframe.index.values.astype('d')*pixel_size
    # Column positions represent -x from our cad model
    x_cad_coordinates = x_scan_coordinates - 256.54
    y_cad_coordinates = y_scan_coordinates - 107.95

    filename_matchobj = re.match(r"""11x8_10MHz_ccw_cscan_export_delam_(\d+)-(\d+)_binary.xlsx""",os.path.split(damage_filename)[1])
    if filename_matchobj is None:
        raise ValueError("Failed to match filename pattern for file \"%s\"" % (filename_matchobj))
    startlayer=int(filename_matchobj.group(1))
    endlayer=int(filename_matchobj.group(2))
    assert(endlayer==startlayer+1)
    delam_layernum = startlayer
    linenumber = linenumbers[delam_layernum-1] # delamination 1-2 would correspond to our first delamination (linenumbers[0])
    
    if delam_layernum < numLayers:
        delam_depth = delam_layernum*thickness1
        pass
    else:
        delam_depth = numLayers*thickness1 + (delam_layernum-numLayers)*thickness2
        pass

    delam_coordindexes = np.where(damage_matrix)
    delam_y_coords = y_cad_coordinates[delam_coordindexes[0]]
    delam_x_coords = x_cad_coordinates[delam_coordindexes[1]]
    
    delam_coords = np.array((delam_x_coords,delam_y_coords),dtype='d') # ,np.ones(delam_outline_raw.shape[1],dtype='d')*delam_depth)

    # This next bit of code converts an array of points inside the delaminated region to the delamination outline. It assumes the delamination is convex, which is true for the data of this example.
    delam_qhull = scipy.spatial.ConvexHull(delam_coords.T)
    delam_hull_outline = np.array((delam_coords[0,delam_qhull.vertices],delam_coords[1,delam_qhull.vertices],np.ones(delam_qhull.vertices.shape[0],dtype='d')*delam_depth),dtype='d')
    
    delam_hull_outline_closed=np.concatenate((delam_hull_outline,delam_hull_outline[:,0:1]),axis=1)

    ## The outlines seem to have x and y swapped relative to the CAD model. This fixes the outlines
    #delam_hull_outline_closed_swapxy=np.concatenate((delam_hull_outline_closed[1:2,:],delam_hull_outline_closed[0:1,:],delam_hull_outline_closed[2:3,:]),axis=0)
    

    linenumber = linenumbers[int(delam_layernum)-1] # line number of bond_layers call

    
    csvfilename = os.path.join(output_directory,"layerdelam_%5.5d_%2.2d.csv" % (linenumber,delamcounts[int(delam_layernum)-1]))
    
    # !!!!*** Temporary ignore all delams after the first in each layer because they all overlap
    if delamcounts[int(delam_layernum)-1] > 1:
        continue
    
    delamcounts[int(delam_layernum)-1] += 1

    

    #np.savetxt(csvfilename,delam_hull_outline_closed_swapxy.T,delimiter=", ")
    np.savetxt(csvfilename,delam_hull_outline_closed.T,delimiter=", ")
    

    if doplots:
        pl.figure(delam_layernum)
        #pl.plot(delam_hull_outline_closed_swapxy[0,:],delam_hull_outline_closed_swapxy[1,:],'-')
        pl.plot(delam_hull_outline_closed[0,:],delam_hull_outline_closed[1,:],'-')
        pl.xlabel('X position (mm)')
        pl.ylabel('Y position (mm)')
        pl.title('Layer boundary #%d' % (delam_layernum))
        oldlegends=[]
        if delam_layernum in legends:
            oldlegends=copy.copy(legends[delam_layernum])
            pass
        oldlegends.append(os.path.split(damage_filename)[1])
        pl.legend(oldlegends)
        legends[delam_layernum]=oldlegends
        
        pass
    
    pass

if doplots:
    pl.show()
    pass
