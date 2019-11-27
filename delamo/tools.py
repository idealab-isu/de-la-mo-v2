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

def ProjectEdgesOntoFace(edge_edges, face):
    edge_curves = [BRep_Tool.Curve(edge) for edge in edge_edges]

    surface = BRep_Tool.Surface(face)

    # Note: edge_curves[i][1] and edge_curves[i][2] appear to be start and end u coordinates for curve

    Projections = [GeomProjLib.geomprojlib_Project(edge_curve[0], surface) for edge_curve in edge_curves]

    # Right here we should be trimming our projection to line up with layerbodyface1 and the unprojected edge (element of edge_edges)
    # But it's probably OK not to, because we are using the projection to make a tool that will be used to cut the face
    # and the extension of the tool beyond the face boundary shouldn't cause any problems, at least so long as thath
    # geometry doesn't get too weird

    ProjectionEdges = [BRepBuilderAPI.BRepBuilderAPI_MakeEdge(Projection).Edge() for Projection in Projections]

    # If we did trimmming, we would need to construct wire from the edge(s) that actually projected to something within the face,
    # with any gaps filled by appropriately trimmed edges from the face.

    # ProjectedWireBuilder = BRepBuilderAPI.BRepBuilderAPI_MakeWire()
    #
    # for ProjectionEdge in ProjectionEdges:
    #    ProjectedWireBuilder.add(ProjectionEdge)
    #    pass
    #
    # ProjectedWire = ProjectedWireBuilder.Wire()

    # Need to take ProjectionEdges, which are located on layerbodysurface1
    # and perform a offset from the surface
    return ProjectionEdges

# Parametric space search start location for FindOCCPointNormal()
FindOCCPointNormal_refParPoint = np.array([0.05,0.1])

def FindOCCPointNormal(Face, PointTolerance, NormalTolerance):
    """ Given a face, find and return a (point, normal, parPoint)
    that uniquely identifies the face. The point must be 
    significantly farther than tolerance from any edge"""

    # breptools_Write(Face,"/tmp/PointFind.brep")

    # Evaluate reference parametric point
    faceSurface = BRep_Tool().Surface(Face)
    refParPoint = FindOCCPointNormal_refParPoint.copy()
    refPointProps = GeomLProp_SLProps(faceSurface, refParPoint[0], refParPoint[1], 1, PointTolerance)

    faceNormal = gp_Vec(refPointProps.Normal())
    if Face.Orientation == TopAbs_REVERSED:
        # Face is reversed from underlying surface -> we need to flip the normal
        faceNormal = -faceNormal
        pass
    facePoint = refPointProps.Value()
    faceParPoint = refParPoint

    # Check if original reference point is inside tha Face
    FaceExplorer=BRepClass_FaceExplorer(Face)
    C=BRepClass_FClassifier()
    C.Perform(FaceExplorer, gp_Pnt2d(refParPoint[0],refParPoint[1]), PointTolerance)

    # If original point is not inside face iterate to find a point inside face
    if (C.State()!=TopAbs_IN):

        origPoint = facePoint
        origPointVertex = BRepBuilderAPI.BRepBuilderAPI_MakeVertex(origPoint).Vertex()
        origParPoint = refParPoint

        # Find the closest point by this method:
        # https://www.opencascade.com/content/closest-point-step-object
        DistanceCalculator = BRepExtrema_DistShapeShape(Face, origPointVertex)
        DistanceCalculator.Perform()
        currentDist = DistanceCalculator.Value()

        if DistanceCalculator.NbSolution() > 0:

            closestDist = currentDist
            # Evaluate (u,v) coordinates on this face of closest point
            # !!!*** ParOnFaceS1 seems to fail sometimes
            # because OCC finds an edge closer.
            # In that case you can call DistanceCalculator.ParOnEdgeS1(1)
            # to get the coordinate along the edge.
            # You can identify which type the support is
            # by calling DistanceCalculator.SupportOnShape(1).ShapeType()
            # and comparing with TopAbs_EDGE, etc.
            # From the documentation it looks like you might get
            # TopAbs_VERTEX, as well

            # Alternate method that might work is calculate the 3D point and
            # then calculate the parametric value on the face
            #(ClosestU, ClosestV) = DistanceCalculator.ParOnFaceS1(1)

            currentCP = DistanceCalculator.PointOnShape1(1)
            SAS = ShapeAnalysis_Surface(faceSurface)
            currentUV = SAS.ValueOfUV(currentCP, PointTolerance)
            ClosestU = currentUV.X()
            ClosestV = currentUV.Y()
            [uMin, uMax, vMin, vMax] = SAS.Bounds()

            angleIncrement = 1
            parIncrement = 0.001 * (math.sqrt((uMax-uMin)*(uMax-uMin) + (vMax-vMin)*(vMax-vMin))/math.sqrt(2.00))
            if (parIncrement > 100):
                print("WARNING: Face bounds of Face %s are incorrect! Reverting to default parametric value"%Face)
                parIncrement = 0.001
            pointFound = False
            for angle in range(0,359,angleIncrement):
                newU = ClosestU + parIncrement * math.cos(angle*math.pi/180.0)
                newV = ClosestV + parIncrement * math.sin(angle*math.pi/180.0)

                newParPoint = np.array([newU, newV])
                C.Perform(FaceExplorer, gp_Pnt2d(newParPoint[0], newParPoint[1]), PointTolerance)
                if (C.State() == TopAbs_IN):
                    # Evaluate reference parametric point
                    newPointProps = GeomLProp_SLProps(faceSurface, newParPoint[0], newParPoint[1], 1,
                                                      PointTolerance)
                    faceNormal = gp_Vec(newPointProps.Normal())
                    facePoint = newPointProps.Value()
                    faceParPoint = newParPoint
                    if Face.Orientation == TopAbs_REVERSED:
                        # Face is reversed from underlying surface -> we need to flip the normal
                        faceNormal = -faceNormal
                        pass
                    pointFound = True
                    break
                pass

            if (not pointFound):
                raise ValueError("Point inside face not found!")

            pass

        pass

    # (vs. we get TopAbs_ON if it is on an edge rather the inside of the face)
        
    #    print('Point in Face')
    # sys.modules["__main__"].__dict__.update(globals())
    # sys.modules["__main__"].__dict__.update(locals())
    # raise ValueError("Break")

    return (np.array((facePoint.X(), facePoint.Y(), facePoint.Z()), dtype='d'),
             np.array((faceNormal.X(), faceNormal.Y(), faceNormal.Z()), dtype='d'),
             np.array((faceParPoint[0], faceParPoint[1]), dtype='d'))

def SelectFaceByPointNormal(InputTopology,Point,Normal,PointTolerance,NormalTolerance):
    """ Select a face from InputTopology that contains Point on its interior within PointTolerance.
    Does not currently use Normal Parameter. Returns a topods_Face or None"""
    # For the moment we don't actually use the normal...

    #import pdb
    #pdb.set_trace()
    FaceExp=TopExp_Explorer(InputTopology,TopAbs_FACE)
    while FaceExp.More():
        # Extract the Surface object (geometry, not topology) underlying this face
        tds_Face=topods_Face(FaceExp.Current())

        if OCCPointInFace(Point,tds_Face,PointTolerance)==TopAbs_IN:
            return tds_Face
        FaceExp.Next()
        pass
    return None
        


def OCCPointInFace(Point, Face, PointTolerance,debug=False):
    """ Given a face and a Point check if point lies on the face"""

    facePoint = gp_Pnt(Point[0],Point[1],Point[2])
    # faceSurface = BRep_Tool().Surface(Face)
    # SAS = ShapeAnalysis_Surface(faceSurface)
    # currentUV = SAS.ValueOfUV(facePoint, PointTolerance)
    # parPointU = currentUV.X()
    # parPointV = currentUV.Y()
    # print(parPointU, parPointV)
    #
    #FaceExplorer = BRepClass_FaceExplorer(Face)
    #C = BRepClass_FClassifier()
    # C.Perform(FaceExplorer, gp_Pnt2d(parPointU, parPointV), PointTolerance)

    #C.Perform(Face, facePoint, PointTolerance )

    #return C.State()

    # Find the closest point by this method
    # https://www.opencascade.com/content/closest-point-step-object


    origPointVertex = BRepBuilderAPI.BRepBuilderAPI_MakeVertex(facePoint).Vertex()

    DistanceCalculator = BRepExtrema_DistShapeShape(Face, origPointVertex)
    DistanceCalculator.Perform()
    if not DistanceCalculator.IsDone():
        raise ValueError("Point to face distance calculation failed")
    currentDist = DistanceCalculator.Value()

    if debug:
        print("currentDist=%f" % (currentDist))
        pass
    
    if (currentDist < PointTolerance):
        return TopAbs_IN
    return TopAbs_OUT

