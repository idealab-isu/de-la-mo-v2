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
import os
import os.path

# Old pythonocc defined core functionality directly in OCC
# New pythonocc defines it in OCC.core
use_OCC_core=False
try:
    import OCC.TopoDS
    pass
except ImportError:
    use_OCC_core=True
    pass

if use_OCC_core:
    from OCC.Core.TopoDS import topods
    from OCC.Core.TopoDS import TopoDS_Face
    from OCC.Core.TopoDS import TopoDS_Shape
    from OCC.Core.BRep import BRep_Builder
    from OCC.Core.BRep import BRep_Tool
    from OCC.Core.BRepTools import breptools_Read
    from OCC.Core.TopExp import TopExp_Explorer
    from OCC.Core.TopAbs import TopAbs_FACE

    from OCC.Core.STEPControl import STEPControl_Reader
    from OCC.Core.IGESControl import IGESControl_Reader
    from OCC.Core.IFSelect import IFSelect_RetDone, IFSelect_ItemsByEntity
    #from OCC.Core.Display.SimpleGui import init_display
    pass
else:
    from OCC.TopoDS import topods
    from OCC.TopoDS import TopoDS_Face
    from OCC.TopoDS import TopoDS_Shape
    from OCC.BRep import BRep_Builder
    from OCC.BRep import BRep_Tool
    from OCC.BRepTools import breptools_Read
    from OCC.TopExp import TopExp_Explorer
    from OCC.TopAbs import TopAbs_FACE

    from OCC.STEPControl import STEPControl_Reader
    from OCC.IGESControl import IGESControl_Reader
    from OCC.IFSelect import IFSelect_RetDone, IFSelect_ItemsByEntity
    #from OCC.Display.SimpleGui import init_display
    pass


def load_brep(BRepFile):
    BRepShape=TopoDS_Shape()
    builder=BRep_Builder()
    breptools_Read(BRepShape,BRepFile,builder)

    return BRepShape

def load_step_singleshape(filename):
    step_reader = STEPControl_Reader()
    status = step_reader.ReadFile(filename)
    if status != IFSelect_RetDone:
        raise IOError("STEP reader status = %s" % (str(status)))

    step_reader.TransferRoot(1)  # Not sure what this does
    
    if step_reader.NbShapes() > 1:
        raise ValueError("STEP file contains multiple shapes")

    return step_reader.Shape(1)
    
    

def load_iges(filename):

    iges_reader = IGESControl_Reader()
    status = iges_reader.ReadFile(filename)
    if status != IFSelect_RetDone:
        raise IOError("IGES reader status = %s" % (str(status)))

    iges_reader.TransferRoots()  # Not sure what this does

    BRepShape=iges_reader.Shape(1)
    
    return BRepShape

def load_byfilename(filename):
    ext = os.path.splitext(filename)[1].lower()

    if ext==".brep" or ext==".brp":
        return load_brep(filename)
    elif ext==".step" or ext==".stp":
        return load_step_singleshape(filename)
    elif ext==".iges" or ext==".igs":
        return load_iges(filename)
    raise ValueError("Unknown filename extension %s" % (ext))
    


def GetFacesSurfaces(BRepShape):
    FaceExplorer=TopExp_Explorer(BRepShape,TopAbs_FACE)

    Faces = []
    while FaceExplorer.More():
        ShapeFace=FaceExplorer.Current() # a TopoDS_Shape
        # Convert TopoDS_Shape to a TopoDS_Face
        Face=topods.Face(ShapeFace)
        
        Faces.append(Face)
        FaceExplorer.Next()
        pass

    Surfaces = [  BRep_Tool().Surface(Face) for Face in Faces ] 
    SurfObjs = [ Surface.GetObject() for Surface in Surfaces ]

    return (Faces,Surfaces,SurfObjs)
