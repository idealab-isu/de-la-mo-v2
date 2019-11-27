from __future__ import print_function

import sys
import os
import os.path
import csv

import numpy as np
import math

from OCC.TopoDS import topods

from OCC.TopoDS import TopoDS_Face
from OCC.TopoDS import TopoDS_Shape
from OCC.TopoDS import TopoDS_Compound
from OCC.TopoDS import TopoDS_Wire
from OCC.TopoDS import TopoDS_Vertex
from OCC.TopoDS import TopoDS_Shell
from OCC.TopoDS import TopoDS_Solid
from OCC.TopoDS import topods_Shell
from OCC.TopoDS import topods_Face
from OCC.TopoDS import topods_Edge
from OCC.TopoDS import topods_Vertex
from OCC.BRep import BRep_Builder
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeEdge
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeWire
from OCC.BRep import BRep_Tool
from OCC.BRepExtrema import BRepExtrema_DistShapeShape
from OCC import BRepLib
from OCC.BRepLib import BRepLib_FuseEdges
from OCC import BRepOffsetAPI
from OCC import BRepOffset
from OCC import BRepFeat
from OCC import BRepBuilderAPI
#from OCC.BRepClass import BRepClass_FacePassiveClassifier
from OCC.BRepClass import BRepClass_FaceExplorer
from OCC.BRepClass import BRepClass_FClassifier
from OCC.ShapeAnalysis import ShapeAnalysis_FreeBoundsProperties
from OCC.ShapeAnalysis import ShapeAnalysis_FreeBounds
from OCC.ShapeAnalysis import ShapeAnalysis_Surface
from OCC.BRepTools import breptools_Read
from OCC.TopExp import TopExp_Explorer
from OCC.TopAbs import TopAbs_ON
from OCC.TopAbs import TopAbs_IN
from OCC.TopAbs import TopAbs_OUT
from OCC.TopAbs import TopAbs_FACE
from OCC.TopAbs import TopAbs_VERTEX
from OCC.TopAbs import TopAbs_EDGE
from OCC.TopAbs import TopAbs_SHELL
from OCC.TopAbs import TopAbs_SOLID
from OCC.TopAbs import TopAbs_FORWARD
from OCC.TopAbs import TopAbs_REVERSED
from OCC.GeomAbs import GeomAbs_Arc
from OCC.TopTools import TopTools_ListIteratorOfListOfShape
from OCC.TopTools import TopTools_ListOfShape
from OCC.TopoDS import TopoDS_Iterator
from OCC.GeomLProp import GeomLProp_SLProps
from OCC.gp import gp_Pnt2d
from OCC.gp import gp_Vec
from OCC.gp import gp_Dir
from OCC.gp import gp_Pnt
from OCC.GEOMAlgo import GEOMAlgo_Splitter
from OCC import GeomProjLib
from OCC.TColgp import TColgp_Array1OfPnt
from OCC.TColgp import TColgp_HArray1OfPnt
from OCC.GeomAPI import (GeomAPI_Interpolate, GeomAPI_PointsToBSpline)
from OCC.BRepAlgoAPI import BRepAlgoAPI_Section
from OCC.GProp import GProp_GProps
from OCC.BRepGProp import brepgprop_SurfaceProperties
from OCC.BRepGProp import brepgprop_VolumeProperties

from OCC.BRepMesh import BRepMesh_IncrementalMesh
from OCC.Bnd import Bnd_Box
from OCC.BRepBndLib import brepbndlib_Add
from OCC.TopLoc import TopLoc_Location

from OCC.BRepTools import breptools_Read
from OCC.BRepTools import breptools_Write
from OCC.STEPControl import STEPControl_Reader
from OCC.STEPControl import STEPControl_Writer
from OCC.STEPControl import STEPControl_ShellBasedSurfaceModel
from OCC.STEPControl import STEPControl_ManifoldSolidBrep
from OCC.STEPControl import STEPControl_GeometricCurveSet
from OCC.IGESControl import IGESControl_Reader
from OCC.IFSelect import IFSelect_RetDone, IFSelect_ItemsByEntity


from OCC.BRepAlgoAPI import BRepAlgoAPI_BooleanOperation
from OCC.BRepAlgoAPI import BRepAlgoAPI_Cut
from OCC.BOPAlgo import BOPAlgo_CUT


from .tools import ProjectEdgesOntoFace,FindOCCPointNormal,SelectFaceByPointNormal,OCCPointInFace

writecount = 0

class Solid(object):
    """ Represents a solid object. It is mutable, but only in certain ways
    (e.g. by subtraction, etc.). It will be an error if any bonded faces
    fail to exist in a successor object. """
    
    Name = None  # Name by which this solid is identified (for assigning boundary conditions, etc) 
    ImmutableSolid=None # ImmutableSolid instance

    # define .Shape property for consistency with LayerBody to pull
    # out the OCC Shape
    @property
    def Shape(self):
        return self.ImmutableSolid.Solid
    
    def __init__(self,**kwargs):

        for argname in kwargs:
            if not hasattr(self,argname):
                raise ValueError("Invalid attribute")
            setattr(self,argname,kwargs[argname])
            pass
        pass

    @classmethod
    def FromOCC(cls,Name,SolidShape,PointTolerance=1e-5,NormalTolerance=1e-6):
        _ImmutableSolid = ImmutableSolid.FromOCC(SolidShape,PointTolerance=PointTolerance,NormalTolerance=NormalTolerance)
        return cls(Name=Name,ImmutableSolid=_ImmutableSolid)


    def SubtractLayer(self, layer,PointTolerance=1e-5,NormalTolerance=1e-6):
        """Subtract the layer from the inputSolid
        """

        # Cut the layer solid model from the InputSolid

        #BooleanOp = BRepAlgoAPI_BooleanOperation()
        #BooleanOp.SetOperation(BOPAlgo_CUT)

        BooleanOp = BRepAlgoAPI_Cut()

        BooleanOp_ArgumentShapes = TopTools_ListOfShape()
        BooleanOp_ArgumentShapes.Append(self.ImmutableSolid.Solid)
        BooleanOp.SetArguments(BooleanOp_ArgumentShapes)
        
        BooleanOp_ToolShapes = TopTools_ListOfShape()
        
        for layerbody in layer.BodyList:
            BooleanOp_ToolShapes.Append(layerbody.Shape)
            pass

        BooleanOp.SetTools(BooleanOp_ToolShapes)

        BooleanOp.Build()

        if BooleanOp.ErrorStatus() != 0:
            raise ValueError("Error in subtraction boolean operation")

        BooleanResult=BooleanOp.Shape()

        global writecount
        
        step_writer2=STEPControl_Writer()
        step_writer2.Transfer(BooleanResult,STEPControl_ManifoldSolidBrep,True)
        step_writer2.Write("/tmp/BooleanResult%d.step" % (writecount))
        writecount +=1
        
        self.ImmutableSolid = ImmutableSolid.FromOCC(BooleanResult,PointTolerance=PointTolerance,NormalTolerance=NormalTolerance)

        pass
    

    def layer_adjacent_side_faces(self,layer,PointTolerance=1e-5,NormalTolerance=1e-6):
        """ return a face adjacency list indicating the faces within this solid that
        map to (by point/normal identification) faces within the layer. 

        NOTE: If the layer extends beyond the surface of the solid then the face in the layer
        and the face in the solid will not be the same, because the surface of the solid does not
        get imprinted onto the side of the layer. As a result such an adjacency may not be 
        identified by this function depending on the locations of the points and thus such layers 
        may not be successfully bonded (TIE boundary condition) to the solid. A warning messages
        is generated in this case. As long as the number of unbonded layers is much smaller than
        the total number of layers the influence of missing one or (worst-case) two tied layers is 
        minimal then this is not a problem.
        """
        
        # Iterate through layer and identify side faces of the layer that matches with the faces
        # of the layer


        FAL = []
        


        #exp = TopExp_Explorer(self.ImmutableSolid.Solid, TopAbs_FACE)
        #
        ## iterate over all faces
        #
        #
        #while exp.More():
        #    SolidFace=exp.Current()
        #
        #    
        #    for layerbody in layer.BodyList:
        #        for layerbodyface in layerbody.FaceListSide:
        #
        #            print("Checking %s against %s" % (str(layerbodyface.Face),str(SolidFace)))
        #            
        #            if layerbodyface.Face.IsSame(topods_Face(SolidFace)):
        #                # Matches a side face
        #                FAL.append({
        #                    "name1": self.Name,
        #                    "name2": layerbody.Name,
        #                    "bcType": "TIE",
        #                    "point1": layerbodyface.Point,
        #                    "normal1": layerbodyface.Normal,
        #                    })
        #                break
        #            pass
        #        pass
        #    exp.Next()
        #    pass


        for layerbody in layer.BodyList:
            for layerbodyface in layerbody.FaceListSide:
                SolidFace = SelectFaceByPointNormal(self.ImmutableSolid.Solid,layerbodyface.Point,layerbodyface.Normal,PointTolerance,NormalTolerance)

                if SolidFace is not None:
                    # Matches a face in the solid
                    FAL.append({
                        "name1": self.Name,
                        "name2": layerbody.Name,
                        "bcType": "TIE",
                        "point1": layerbodyface.Point,
                        "normal1": layerbodyface.Normal,
                    })

                    pass
                else:
                    print("WARNING: Side face from layer %s not matched in solid %s (if this happens for a large percentage of layers, it is a problem!)" % (layer.Name,self.Name))
                    pass
                pass
            pass
        return FAL

    pass


class ImmutableSolid(object):
    Solid=None # TopoDS_Solid
    Faces=None # List of all ImmutableSolidFaces for the solid

    def __init__(self,**kwargs):

        for argname in kwargs:
            if not hasattr(self,argname):
                raise ValueError("Invalid attribute")
            setattr(self,argname,kwargs[argname])
            pass
        pass

    @classmethod
    def FromOCC(cls,SolidShape,PointTolerance=1e-5,NormalTolerance=1e-6):
        """Create an ImmutableSolid given a topods_Shape (SolidShape)"""

        # Create Face list by exploring the topods_Shape
        Faces=[]        
        FaceExp = TopExp_Explorer(SolidShape,TopAbs_FACE)
        while FaceExp.More():
            
            Faces.append(ImmutableSolidFace.FromOCC(FaceExp.Current(),PointTolerance=PointTolerance,NormalTolerance=NormalTolerance))
            FaceExp.Next()
            pass
        
        return cls(Solid=SolidShape,
                   Faces=Faces)

    pass

class ImmutableSolidFace(object):
    # Note members that largely match definition of LayerBodyFace in layer.py
    
    Face=None # topods_Face
    Point=None
    Normal=None
    ParPoint=None

    def __init__(self,**kwargs):

        for argname in kwargs:
            if not hasattr(self,argname):
                raise ValueError("Invalid attribute")
            setattr(self,argname,kwargs[argname])
            pass
        pass

    def __eq__(self,other):
        # Test equality of faces. Designed to be interoperable
        # with LayerBodyFace from layer.py
        #
        return self.Face.IsSame(other.Face)
    # Implement not equal operator
    def __ne__(self,other):
        return not(self==other)
    

    # Implement hash operator as a constant.
    # This makes us a hashable type (can work in a dictionary),
    # even though there's not a good way to hash this
    # (though with the current semi-broken equality operator
    # we might be able to get away with using the id of the underlying surface).
    #
    # Returning a constant here means that dictionaries indexed by ImmutableSolidFace
    # will perform like lists (but at least you will be able to do direct lookups)
    def __hash__(self):
        return 0

    @classmethod
    def FromOCC(cls,Face,PointTolerance=1e-5,NormalTolerance=1e-6):
        """ Create an ImmutableSolidFace from a TopoDS_Face. 
        IsPointingInside should be None (if this LayerBodyFace is 
        not actually part of a LayerBody, or if the OCC-Evaluated normal
        is guaranteed to be outward-facing) or a function/method/lambda 
        that, given a Normal, returns whether  that normal is pointing inside
        the LayerBody."""

        Face=topods_Face(Face) # In case it is a TopoDS_Shape of type Face, cast it to a TopoDS_Face.

        (Point,Normal,ParPoint)=FindOCCPointNormal(Face,PointTolerance=PointTolerance,NormalTolerance=NormalTolerance)
        
        
        return cls(Face=Face,
                   Point=Point,
                   Normal=Normal,
                   ParPoint=ParPoint)
    
   
    pass



    
