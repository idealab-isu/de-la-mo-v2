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

from .tools import ProjectEdgesOntoFace,FindOCCPointNormal,SelectFaceByPointNormal,OCCPointInFace

class Solid(object):
    """ Represents a solid object. It is mutable, but only in certain ways
    (e.g. by subtraction, etc.). It will be an error if any bonded faces
    fail to exist in a successor object. """
    
    Name = None  # Name by which this solid is identified (for assigning boundary conditions, etc) 
    ImmutableSolid=None # ImmutableSolid instance

    def __init__(self,**kwargs):

        for argname in kwargs:
            if not hasattr(self,argname):
                raise ValueError("Invalid attribute")
            setattr(self,argname,kwargs[argname])
            pass
        pass

    @classmethod
    def FromOCC(cls,SolidShape,PointTolerance=1e-5,NormalTolerance=1e-6):
        _ImmutableSolid = ImmutableSolid.FromOCC(SolidShape,PointTolerance=PointTolerance,NormalTolerance=NormalTolerance)
        return cls(ImmutableSolid=_ImmutableSolid)

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



    
