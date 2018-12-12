import sys
import os
import os.path

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
