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

from . import loaders
from .tools import ProjectEdgesOntoFace,FindOCCPointNormal,SelectFaceByPointNormal,OCCPointInFace
from . import solid

def OffsetFaceMoreAlongDirection(face, offsetDir, PointTolerance):
    assert(offsetDir == 1 or offsetDir == -1)
    OffsetDist = 1000.0 * PointTolerance * offsetDir

    mkOffset = BRepOffsetAPI.BRepOffsetAPI_MakeOffsetShape(face, OffsetDist, PointTolerance,
                                                            BRepOffset.BRepOffset_Skin,
                                                            False, False,
                                                            GeomAbs_Arc)
    assert (mkOffset.IsDone())

    OffsetShell = mkOffset.Shape()

    OffsetShellFacesExp = TopExp_Explorer(OffsetShell, TopAbs_FACE)
    OffsetFaces = []

    while OffsetShellFacesExp.More():
        OffsetFaces.append(topods_Face(OffsetShellFacesExp.Current()))
        OffsetShellFacesExp.Next()
        pass

    assert (len(OffsetFaces) == 1)  # Offset of a single face should give a single face

    return OffsetFaces[0]


def VectorNormalize(v):
    return v / np.linalg.norm(v)


def Distance(a, b):
    return np.sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]))


def CreateMidVertex(v1, v2):
    v01 = DMVertex()
    v01.point = (v1.point + v2.point) / 2.0
    v01.normal = (v1.normal + v2.normal) / 2.0
    VectorNormalize(v01.normal)
    return v01


def Tessellate(Object, Shape, MeshSize):
    B = Bnd_Box()
    brepbndlib_Add(Shape, B)
    bBoxMinX, bBoxMinY, bBoxMinZ, bBoxMaxX, bBoxMaxY, bBoxMaxZ = B.Get()

    Object.bBoxMin[0] = bBoxMinX
    Object.bBoxMin[1] = bBoxMinY
    Object.bBoxMin[2] = bBoxMinZ
    Object.bBoxMax[0] = bBoxMaxX
    Object.bBoxMax[1] = bBoxMaxY
    Object.bBoxMax[2] = bBoxMaxZ

    currentModelSize = Distance(np.array([bBoxMinX, bBoxMinY, bBoxMinZ]), np.array([bBoxMaxX, bBoxMaxY, bBoxMaxZ]))
    meshMaxSize = MeshSize * currentModelSize

    mesher = BRepMesh_IncrementalMesh()
    mesher.SetShape(Shape)
    mesher.SetAngle(0.01)
    mesher.SetControlSurfaceDeflection(True)
    mesher.SetDeflection(0.001 * currentModelSize)
    mesher.SetParallel(True)
    mesher.SetInternalVerticesMode(True)
    mesher.SetMinSize(0.005 * currentModelSize)
    mesher.SetRelative(False)
    mesher.Perform()

    return meshMaxSize


class DMObject:

    def __init__(self, **kwargs):
        self.totalNumTriangles = 0
        self.objID = None

        self.bBoxMin = np.zeros((3))
        self.bBoxMax = np.zeros((3))

        self.faces = []

        self.filepath = ""

    def RefineTessellation(self, minEdgeLength):
        DELTA = 1e-3
        numTrianglesAdded = 0
        maxObjectEdgeLen = 0
        # k = 6
        for k in range(0, len(self.faces)):
            continueRefining = True
            face = self.faces[k]
            while continueRefining:
                j = 1
                maxFaceEdgeLen = 0.0
                newMaxEdgeLen = 0.0
                refined = False
                numTriangles = len(face.triangles)
                while j < numTriangles:
                    t = face.triangles[j-1]
                    edgeLen01 = Distance(t.vertices[0].point, t.vertices[1].point)
                    edgeLen12 = Distance(t.vertices[1].point, t.vertices[2].point)
                    edgeLen20 = Distance(t.vertices[2].point, t.vertices[0].point)
                    maxEdgeLen = max(max(edgeLen01, edgeLen12), edgeLen20)
                    if maxEdgeLen > minEdgeLength:
                        t1 = DMTriangle()
                        t2 = DMTriangle()

                        t1.facenormal = t.facenormal
                        t1.visibilityFactor = 0
                        t2.facenormal = t.facenormal
                        t1.visibilityFactor = 0

                        if abs(edgeLen01 - maxEdgeLen) < DELTA:
                            mid01 = CreateMidVertex(t.vertices[0], t.vertices[1])
                            t1[0] = t[0]
                            t1[1] = mid01
                            t1[2] = t[2]

                            t2[0] = mid01
                            t2[1] = t[1]
                            t2[2] = t[2]
                            newMaxEdgeLen = max(max(edgeLen01 / 2.0, edgeLen12), edgeLen20)
                        elif abs(edgeLen12 - maxEdgeLen) < DELTA:
                            mid12 = CreateMidVertex(t.vertices[1], t.vertices[2])

                            t1[0] = t[0]
                            t1[1] = t[1]
                            t1[2] = mid12

                            t2[0] = t[0]
                            t2[1] = mid12
                            t2[2] = t[2]
                            newMaxEdgeLen = max(max(edgeLen01, edgeLen12 / 2.0), edgeLen20)
                        elif abs(edgeLen20 - maxEdgeLen) < DELTA:
                            mid20 = CreateMidVertex(t.vertices[2], t.vertices[0])

                            t1[0] = t[0]
                            t1[1] = t[1]
                            t1[2] = mid20

                            t2[0] = mid20
                            t2[1] = t[1]
                            t2[2] = t[2]
                            newMaxEdgeLen = max(max(edgeLen01, edgeLen12), edgeLen20 / 2.0)

                        face.triangles.append(t1)
                        face.triangles.append(t2)
                        del face.triangles[j-1]
                        numTriangles += 1
                        self.totalNumTriangles += 1
                        numTrianglesAdded += 1
                        j -= 1
                        refined = True

                        if numTrianglesAdded % 100 == 0:
                            print(".", end=" ")
                    else:
                        newMaxEdgeLen = maxEdgeLen

                    if not refined:
                        continueRefining = False

                    j += 1
                if newMaxEdgeLen > maxFaceEdgeLen:
                    maxFaceEdgeLen = newMaxEdgeLen
            if maxFaceEdgeLen > maxObjectEdgeLen:
                maxObjectEdgeLen = maxFaceEdgeLen
            print("")

        print("Maximum edge length set       : %s" % minEdgeLength)
        print("Maximum edge length in object : %s" % maxObjectEdgeLen)
        print("Object refined with %s triangles added." % numTrianglesAdded)

    def SaveSTL(self, filepath):
        self.filepath = filepath
        f = open(self.filepath, "w")
        
        f.write("solid %s\n" % (os.path.splitext(os.path.split(filepath)[1])[0]))
        
        for faceNum in range(0, len(self.faces)):
            currentFace = self.faces[faceNum]
            for trinum in range(0, len(currentFace.triangles)):
                f.write("facet normal {} {} {}\n".format(currentFace.triangles[trinum].facenormal[0],
                                                         currentFace.triangles[trinum].facenormal[1],
                                                         currentFace.triangles[trinum].facenormal[2]))
                f.write("outer loop\n")
                for i in range(0, 3):
                    f.write("\tvertex {} {} {}\n".format(currentFace.triangles[trinum].vertices[i].point[0],
                                                         currentFace.triangles[trinum].vertices[i].point[1],
                                                         currentFace.triangles[trinum].vertices[i].point[2]))
                    pass
                f.write("endloop\n")
                f.write("endfacet\n")
                pass
            pass
        
        f.write("endsolid")
        f.close()
        pass
    pass

class DMFace:

    def __init__(self, **kwargs):
        self.dlid = 0
        self.surfID = 0
        self.parentObjID = None
        self.trimmed = None

        self.bBoxMin = np.zeros((3))
        self.bBoxMax = np.zeros((3))

        self.vertexFaces = []
        self.triangles = []

    def GetCommonFace(self, vi1, vi2, face):
        commonFaceIndex = -1
        for p in range(0, len(self.vertexFaces[vi1])):
            faceIndex1 = self.vertexFaces[vi1][p]
            if p != face:
                for q in range(0, len(self.vertexFaces[vi2])):
                    faceIndex2 = self.vertexFaces[vi2][q]
                    if faceIndex1 == faceIndex2:
                        commonFaceIndex = faceIndex1
        return commonFaceIndex


class DMTriangle:

    def __init__(self, **kwargs):
        self.vertices = [None, None, None]
        self.vertexIndex = None
        self.facenormal = None
        self.visibilityFactor = 0

        self.adjacentFaceIndex = np.zeros((3))

    def __getitem__(self, item):
        return self.vertices[item]

    def __setitem__(self, key, value):
        self.vertices[key] = value


class DMVertex:

    def __init__(self, **kwargs):
        self.point = None
        self.normal = np.zeros((3))


class Layer(object):
    """The Layer is a collection of LayerBodies representing 
    a single lamina

    Reference semantics: The Layer is a high-level object that is mutable 
    once created (e.g. fiber breakage)
"""
    
    Name = None
    #Layup = None  # angle in degrees (Don't need this in Geometry Kernel)
    BodyList= None # List of LayerBody objects
    # MoldList = None # List of LayerMold objects (pre-offset) DO WE REALLY NEED TO KEEP THIS? Better to concatentate the FaceListOrig's of the BodyList
    Type = None # "LAMINA" or "STIFFENER"
    # Direction = None # "NODIR", "ORIG", "OFFSET", "SIDE" but "SIDE" or "NODIR" would be inapplicable for a layer
    Thickness = None # Layer thickness
    
    # PairOffset: Layer object on this layer's offset side
    # PairOrig: Layer object on this layer's orig side

    RefMold = None # LayerMold that had been used for offset operation

    DMObj = None
    
    def __init__(self,**kwargs):
        self.BodyList=[]
        self.MoldList=[]
        self.Direction="NODIR"
        
        for argname in kwargs:
            if not hasattr(self,argname):
                raise ValueError("Invalid attribute")
            setattr(self,argname,kwargs[argname])
            pass
        pass

    def bodynames(self):

        return [ Body.Name for Body in self.BodyList ]

    def FindLayerBodyNameByPoint(self,Point,PointTolerance):
        PointVertex = BRepBuilderAPI.BRepBuilderAPI_MakeVertex(gp_Pnt(Point[0],Point[1],Point[2])).Vertex()

        for Body in self.BodyList:

            DistanceCalculator = BRepExtrema_DistShapeShape(Body.Shape, PointVertex)
            DistanceCalculator.Perform()

            if DistanceCalculator.NbSolution() > 0:
                CurrentDist = DistanceCalculator.Value()
                if CurrentDist < PointTolerance:
                    return Body.Name
                pass
            pass
        return None


    def OffsetMold(self):
        """Return a LayerMold based on the OFFSET side of this layer"""
        return LayerMold.FromFaceLists([ Body.FaceListOffset for Body in self.BodyList ])

    def OrigMold(self):
        """Return a LayerMold based on the ORIG side of this layer.
        (Note that the OFFSET direction of this LayerMold points away
        from this layer!) """
        return LayerMold.FromFaceLists([ Body.FaceListOrig for Body in self.BodyList ])

    def Split(self, crackWireFile, Tolerance):
        """Split the layer using the crackWire outline"""

        ReplacementLayerBodyList = []
        for LayerBody in self.BodyList:

            newlayerbodies = LayerBody.Split(crackWireFile,Tolerance)

            for newlayerbody in newlayerbodies:
                newlayerbody.Owner = self
                pass

            ReplacementLayerBodyList.extend(newlayerbodies)
            pass
        self.BodyList = ReplacementLayerBodyList
        pass



    pass


    @classmethod
    def CreateFromMold(cls,Name,Mold,Thickness,Direction,Tolerance,MeshSize=0.5):
        """Create a layer from a LayerMold."""

        if Direction=="OFFSET":
            OffsetDist=Thickness
            pass
        elif Direction=="ORIG":
            OffsetDist=-Thickness
            pass
        else:
            raise ValueError("Invalid Direction: %s" % (str(Direction)))

        # !!!*** When doing the offset operation
        # Should attempt to heal or merge split faces that share the same 
        # underlying geometry. How??? 
        
        mkOffset = BRepOffsetAPI.BRepOffsetAPI_MakeOffsetShape(Mold.Shape, OffsetDist, Tolerance,
                                                               BRepOffset.BRepOffset_Skin,
                                                               False, False,
                                                               GeomAbs_Arc)

        #step_writer2=STEPControl_Writer()
        #step_writer2.Transfer(Mold.Shape,STEPControl_ShellBasedSurfaceModel,True)
        #step_writer2.Write("/tmp/OffsetFace.step")

        assert (mkOffset.IsDone())

        OffsetShell = mkOffset.Shape()

        ## Convert result of offsetting operation (from mkOffset.Shape()) into NURBS
        #OffsetUnconverted = mkOffset.Shape()
        #NurbsConverter = BRepBuilderAPI.BRepBuilderAPI_NurbsConvert(OffsetUnconverted,False)
        #OffsetShell = NurbsConverter.Shape()

        

        # Build list OffsetFaces of faces in OffsetShell
        OffsetShellFacesExp = TopExp_Explorer(OffsetShell,TopAbs_FACE)
        OffsetFaces = []

        while OffsetShellFacesExp.More():
            OffsetFaces.append(topods_Face(OffsetShellFacesExp.Current()))
            OffsetShellFacesExp.Next()
            pass
        

        # Check for outline of original shape
        # FreeBounds = ShapeAnalysis_FreeBounds(Mold.Shape)
        # ClosedWires = FreeBounds.GetClosedWires()

        FreeCheck = ShapeAnalysis_FreeBoundsProperties(Mold.Shape)
        FreeCheck.Perform()

        print("Mold has %d free boundaries." % (FreeCheck.NbClosedFreeBounds()))
        assert (FreeCheck.NbClosedFreeBounds() >= 1)

        # sys.modules["__main__"].__dict__.update(globals())
        # sys.modules["__main__"].__dict__.update(locals())
        # raise ValueError("Break")

        # Build solid from original + offset
        bRepBuilder = BRep_Builder()
        perimeter = TopoDS_Compound()
        bRepBuilder.MakeCompound(perimeter)

        for nBoundary in range(1, FreeCheck.NbClosedFreeBounds() + 1):
            OrigBoundary = FreeCheck.ClosedFreeBound(nBoundary).GetObject()
            OrigWire = OrigBoundary.FreeBound()
            
            # Offset edges (alone?)
            offsetEdge = mkOffset.MakeOffset().OffsetEdgesFromShapes()
            
            # Offset Wire we are creating
            OffsetWire = TopoDS_Wire()
            bRepBuilder.MakeWire(OffsetWire)
            topExplorer = TopExp_Explorer(OrigWire, TopAbs_EDGE)
            
            # iterate over pieces of the original wire
            while topExplorer.More():
                
                if not offsetEdge.HasImage(topExplorer.Current()):
                    raise ValueError("Image of original wire segment is missing")
                
                currentImage = offsetEdge.Image(topExplorer.Current())
                edgeIterator = TopTools_ListIteratorOfListOfShape(currentImage)

                edgeCount = 0
                
                while edgeIterator.More():
                    # Only consider edges
                    if edgeIterator.Value().ShapeType() == TopAbs_EDGE:
                        edgeCount += 1
                        mappedEdge = edgeIterator.Value()
                        pass
                    
                    edgeIterator.Next()
                    pass
                
                if edgeCount != 1:
                    raise ValueError("Got more than one edge from offset operation")
                
                bRepBuilder.Add(OffsetWire, mappedEdge)
                
                topExplorer.Next()
                pass
            
            # Generate side faces
            sideGenerator = BRepOffsetAPI.BRepOffsetAPI_ThruSections()
            sideGenerator.AddWire(OrigWire)
            sideGenerator.AddWire(OffsetWire)
            sideGenerator.Build()
            
            if (not sideGenerator.IsDone()):
                raise ValueError("Side face generation failed\n")

            sideShape = sideGenerator.Shape()
            
            bRepBuilder.Add(perimeter, sideShape)
            pass

        # Now we have to sew all of these pieces together
        # (with a thread!)
        thread = BRepBuilderAPI.BRepBuilderAPI_Sewing()  # sewing tool
        thread.Add(Mold.Shape)
        thread.Add(perimeter)
        thread.Add(OffsetShell)
        thread.Perform()
        
        sewedShell = thread.SewedShape()
        
        # if it is a closed shell, turn it into a solid if possible
        if sewedShell.ShapeType() != TopAbs_SHELL or not sewedShell.Closed():
            raise ValueError("Resulting solid either not a shell or not closed")
        
        # SolidMaker = BRepBuilderAPI.BRepBuilderAPI_MakeSolid(TopoDS.Shell(ResultShape))
        solidMaker = BRepBuilderAPI.BRepBuilderAPI_MakeSolid()
        solidMaker.Add(topods_Shell(sewedShell))
        if not solidMaker.IsDone():
            raise ValueError("Solid maker failed")

        solidShape = solidMaker.Solid()
        
        if not BRepLib.breplib_OrientClosedSolid(solidShape):
            raise ValueError("Solid maker did not yield a closed solid")
        # We successfully got a closed solid

        DMObj = cls.CreateDMObject(Mold.Shape, MeshSize)
        
        #return [layerSolid, offsetSurface]
        NewLayer= cls(Name=Name,
                      Type="LAMINA",
                      Direction=Direction,
                      Thickness=Thickness,
                      RefMold=Mold,
                      DMObj=DMObj)


        NewLayerBody = LayerBody(Name="%s_LB1" % (Name),
                                 Owner=NewLayer,
                                 Shape=solidShape)

        # Extract faces and sort them into NewLayerBody.FaceListOrig, NewLayerBody.FaceListOffset, and NewLayerBody.FaceListSide
        
        # Iterate over all faces
        FaceExp=TopExp_Explorer(solidShape,TopAbs_FACE)
        while FaceExp.More():
            # Extract the Surface object (geometry, not topology) underlying this face
            tds_Face=topods_Face(FaceExp.Current())
            FaceSurf=BRep_Tool.Surface(tds_Face)

            
            # Search for this face in the mold
            MatchedInMold=False
            for MoldFace in Mold.FaceList: # Iterate over LayerBodyFaces in Mold
                if tds_Face.IsSame(MoldFace.Face):
                    if MatchedInMold:
                        raise ValueError("Same face matched twice in mold (!?)")  
                    MatchedInMold=True
                    pass
                pass


            # Search for this face in the offset surface
            MatchedInOffset=False
            for OffsetFace in OffsetFaces: # Iterate over LayerBodyFaces in Mold
                if tds_Face.IsSame(OffsetFace): 
                    if MatchedInOffset:
                        raise ValueError("Same face matched twice in offset (!?)")
                    MatchedInOffset=True
                    pass
                pass


            # Since these faces are extracted from the generated
            # closed object, that the normals generated in FromOCC()
            # should be outward-facing and thus we don't need to
            # provide an IsPointingInside lambda.
            
            if MatchedInMold and not(MatchedInOffset):
                # Create LayerBodyFace
                NewLayerBody.FaceListOrig.append(LayerBodyFace.FromOCC(tds_Face,"ORIG",Owner=NewLayerBody))
                pass
            elif MatchedInOffset and not(MatchedInMold):
                # Create LayerBodyFace
                NewLayerBody.FaceListOffset.append(LayerBodyFace.FromOCC(tds_Face,"OFFSET",Owner=NewLayerBody))
                pass
            elif MatchedInOffset and MatchedInMold:
                raise ValueError("Same surface matched in both offset and mold (!?)")
            else:
                # Must be a side face
                # Create LayerBodyFace
                NewLayerBody.FaceListSide.append(LayerBodyFace.FromOCC(tds_Face,"OFFSET",Owner=NewLayerBody))
                pass
            FaceExp.Next()
            pass
        
        
        NewLayer.BodyList.append(NewLayerBody)

        return NewLayer
    pass

    @classmethod
    def CreateDMObject(cls, solidShape, MeshSize):
        totalNumFaces = 0
        totalNumTriangles = 0
        modelPos = np.zeros((3))
        modelSize = 0.0

        tempObject = DMObject()
        tempObject.totalNumTriangles = 0
        tempObject.objID = 1

        meshMaxSize = Tessellate(tempObject, solidShape, MeshSize)

        faceID = 0

        # Iterate over all faces
        FaceExp=TopExp_Explorer(solidShape,TopAbs_FACE)

        while FaceExp.More():
            # Extract the Surface object (geometry, not topology) underlying this face
            tds_Face=topods_Face(FaceExp.Current())

            objFace = DMFace()
            objFace.dlid = tempObject.objID * 1000 + faceID + 1
            objFace.surfID = faceID
            objFace.parentObjID = tempObject.objID
            objFace.trimmed = False

            FB = Bnd_Box()

            brepbndlib_Add(tds_Face, FB)
            bBoxMinX, bBoxMinY, bBoxMinZ, bBoxMaxX, bBoxMaxY, bBoxMaxZ = FB.Get()
            faceOffset = 0.0
            objFace.bBoxMin[0] = bBoxMinX - faceOffset
            objFace.bBoxMin[1] = bBoxMinY - faceOffset
            objFace.bBoxMin[2] = bBoxMinZ - faceOffset
            objFace.bBoxMax[0] = bBoxMinX + faceOffset
            objFace.bBoxMax[1] = bBoxMinY + faceOffset
            objFace.bBoxMax[2] = bBoxMinZ + faceOffset

            L = TopLoc_Location()
            faceTriangulation = BRep_Tool().Triangulation(tds_Face, L).GetObject()
            numTriangles = faceTriangulation.NbTriangles()
            numVertices = faceTriangulation.NbNodes()

            tempVertices = []
            if numVertices > 0 and numTriangles > 0:
                tri = faceTriangulation.Triangles()

                for i in range(0, numVertices):
                    currentVertex = faceTriangulation.Nodes().Value(i+1)
                    tempVertex = DMVertex()
                    tempVertex.point = np.array([currentVertex.X(), currentVertex.Y(), currentVertex.Z()])
                    tempVertices.append(tempVertex)
                    objFace.vertexFaces.append([])

                for i in range(0, faceTriangulation.NbTriangles()):
                    trian = tri.Value(i+1)
                    index1, index2, index3 = trian.Get()

                    t = DMTriangle()
                    t.vertexIndex = np.array([index1, index2, index3])

                    t[0] = tempVertices[t.vertexIndex[0]-1]
                    t[1] = tempVertices[t.vertexIndex[1]-1]
                    t[2] = tempVertices[t.vertexIndex[2]-1]

                    side1 = t[1].point - t[0].point
                    side2 = t[2].point - t[0].point

                    t.facenormal = np.cross(side1, side2)
                    area = 0.5*np.linalg.norm(t.facenormal)
                    t.facenormal = VectorNormalize(t.facenormal)

                    tempVertices[t.vertexIndex[0]-1].normal += area * t.facenormal
                    tempVertices[t.vertexIndex[1]-1].normal += area * t.facenormal
                    tempVertices[t.vertexIndex[2]-1].normal += area * t.facenormal

                    t.visibilityFactor = 0

                    triangleNum = len(objFace.triangles)
                    objFace.vertexFaces[t.vertexIndex[0]-1].append(triangleNum)
                    objFace.vertexFaces[t.vertexIndex[1]-1].append(triangleNum)
                    objFace.vertexFaces[t.vertexIndex[2]-1].append(triangleNum)
                    objFace.triangles.append(t)
                print("Face %s triangulated with %s triangles." % (faceID, faceTriangulation.NbNodes()))
            else:
                print("Face %s triangulation failed!" % faceID)

            for i in range(0, len(tempVertices)):
                tempVertices[i].normal = VectorNormalize(tempVertices[i].normal)

            for i in range(0, len(objFace.triangles)):
                vertexIndex = objFace.triangles[i].vertexIndex
                objFace.triangles[i].triangleID = i
                objFace.triangles[i][0].normal = tempVertices[vertexIndex[0] - 1].normal
                objFace.triangles[i][1].normal = tempVertices[vertexIndex[1] - 1].normal
                objFace.triangles[i][2].normal = tempVertices[vertexIndex[2] - 1].normal

                objFace.triangles[i].adjacentFaceIndex[0] = objFace.GetCommonFace(vertexIndex[0] - 1,
                                                                                  vertexIndex[1] - 1, i)
                objFace.triangles[i].adjacentFaceIndex[1] = objFace.GetCommonFace(vertexIndex[1] - 1,
                                                                                  vertexIndex[2] - 1, i)
                objFace.triangles[i].adjacentFaceIndex[2] = objFace.GetCommonFace(vertexIndex[2] - 1,
                                                                                  vertexIndex[0] - 1, i)

                totalNumTriangles += 1
                tempObject.totalNumTriangles += 1

            tempObject.faces.append(objFace)
            del tempVertices
            faceID += 1
            totalNumFaces += 1

            FaceExp.Next()


        oldNumTriangles = tempObject.totalNumTriangles
        tempObject.RefineTessellation(meshMaxSize)
        totalNumTriangles += (tempObject.totalNumTriangles - oldNumTriangles)

        currentModelSize = np.linalg.norm(tempObject.bBoxMax - tempObject.bBoxMin) / 4
        if currentModelSize > modelSize:
            modelSize = currentModelSize
        currentModelPos = np.array([-(tempObject.bBoxMin[0] + tempObject.bBoxMax[0]),
                                    -(tempObject.bBoxMin[1] + tempObject.bBoxMax[1]),
                                    -(tempObject.bBoxMin[2] + tempObject.bBoxMax[2])])
        modelPos = (modelPos + currentModelPos) / 2.0

        return tempObject


class LayerBody(object):
    """ The LayerBody is a solid, defined as a boundary representation
    from a shell which in turn consists of many faces. It represents
    a portion of a layer or stiffener.

    ** MAYBE should make this derive from "Body" class 
    that just has a single list of all faces. 

    Object semantics: The LayerBody is mutable once created. But
    the faces and shape from which it is construct are not mutable
    and may need to be regenerated when the LayerBody is modified.  

"""
    Name = None
    FaceListOrig=None # Faces on the "ORIG" side (LayerBodyFace)
    FaceListOffset=None # Faces on the "OFFSET" side (LayerBodyFace)
    FaceListSide=None # Facees on the "SIDE" sides (LayerBodyFace)
    Owner = None # reference to Layer object of which this LayerBody MAY be a part
    #Mold = None  # LayerMold that was used to generate this body SHOULD WE REALLY HAVE THIS ATTRIBUTE? ALL OF THE RELEVANT INFO SHOULD BE IN FaceListOrig
    Shape = None # OpenCascade TopoDS Shape which is a Solid which consists of a single closed shell which consists of multiple faces.
    fiberorientation = None
    # Constant
    StepModelType = STEPControl_ManifoldSolidBrep # Used for OpenCascade STEP writer
    
    
    def __init__(self,**kwargs):
        self.FaceListOrig=[]
        self.FaceListOffset=[]
        self.FaceListSide=[]
        
        for argname in kwargs:
            if not hasattr(self,argname):
                raise ValueError("Invalid attribute")
            setattr(self,argname,kwargs[argname])
            pass
        pass

    def Rebuild_Shape(self):
        # Build/rebuild the .Shape attribute from the Face Lists: FaceListOrig, FaceListOffset, and FaceListSide

        # Now we have to sew all of these pieces together
        # (with a thread!)
        thread = BRepBuilderAPI.BRepBuilderAPI_Sewing()  # sewing tool

        for Face in self.FaceListOrig:            
            thread.Add(Face.Face)
            pass

        for Face in self.FaceListOffset:            
            thread.Add(Face.Face)
            pass

        for Face in self.FaceListSide:            
            thread.Add(Face.Face)
            pass

        thread.Perform()
        
        sewedShell = thread.SewedShape()
        
        # if it is a closed shell, turn it into a solid if possible
        if sewedShell.ShapeType() != TopAbs_SHELL or not sewedShell.Closed():
            raise ValueError("Resulting solid either not a shell or not closed")
        
        # SolidMaker = BRepBuilderAPI.BRepBuilderAPI_MakeSolid(TopoDS.Shell(ResultShape))
        solidMaker = BRepBuilderAPI.BRepBuilderAPI_MakeSolid()
        solidMaker.Add(topods_Shell(sewedShell))
        if not solidMaker.IsDone():
            raise ValueError("Solid maker failed")

        solidShape = solidMaker.Solid()
        
        if not BRepLib.breplib_OrientClosedSolid(solidShape):
            raise ValueError("Solid maker did not yield a closed solid")
        # We successfully got a closed solid
        self.Shape = solidShape
        pass


    def Split(self,crackWireFile,Tolerance):

        # Load the CSV file of the splitting line and create a face between the ORIG and OFFSET faces of the layer
        crack_wire_pt_list = []
        with open(crackWireFile) as csvfile:
            reader = csv.reader(csvfile, delimiter=',', quotechar='"')
            for row in reader:
                if len(row) != 3:
                    raise ValueError("Malformed row in CSV file %s: %s" % (crackWireFile, ",".join(row)))
                try:
                    x = float(row[0])
                    y = float(row[1])
                    z = float(row[2])
                    crack_wire_pt_list.append((x, y, z))
                    pass
                except ValueError:
                    pass
                pass
            pass

        if len(crack_wire_pt_list) == 0:
            raise ValueError("Could not parse any lines from CSV file %s" % (crackWireFile))

        # If there are n entries in the delam_outlist, one of which is doubled (start and end). There will be n-1 segments
        crack_wirepointsHArray = TColgp_HArray1OfPnt(1, len(crack_wire_pt_list))

        for pos in range(len(crack_wire_pt_list)):
            current_point = gp_Pnt(crack_wire_pt_list[pos][0], crack_wire_pt_list[pos][1], crack_wire_pt_list[pos][2])
            crack_wirepointsHArray.SetValue(pos + 1, current_point)
            pass

        # Interpolate the points to make a closed curve
        interpAPI = GeomAPI_Interpolate(crack_wirepointsHArray.GetHandle(), False, Tolerance)
        interpAPI.Perform()
        if interpAPI.IsDone():
            crack_curve = interpAPI.Curve()
        else:
            raise ValueError("Curve interpolation failed\n")

        # Convert a curve to edge and then to Shape
        crack_edge = BRepBuilderAPI_MakeEdge(crack_curve).Edge()
        WireBuilder = BRepBuilderAPI_MakeWire()
        WireBuilder.Add(crack_edge)
        CrackWireShape = WireBuilder.Shape()

        # step_writer2=STEPControl_Writer()
        # step_writer2.Transfer(CrackWireShape,STEPControl_GeometricCurveSet,True)
        # step_writer2.Write("../data/Wire.STEP")

        exp = TopExp_Explorer(CrackWireShape, TopAbs_EDGE)

        # Iterate over all edges
        edge_shapes = []
        while exp.More():
            edge_shapes.append(exp.Current())

            exp.Next()
            pass
        edge_edges = [topods_Edge(edge_shape) for edge_shape in edge_shapes]


        # Get the offset and orig faces of the layer body and project the edge to both
        # ***!!!! Really we should be iterating over FaceListOffset
        # and FaceListOrig and finding correspondences.
        # This is probably OK for most cases as long as
        # offsetFace and origFace share the same underlying surface
        # and we are just using these to create a tool
        # that will be used to cut the body.

        # There is an issue for curved faces when the tool does not completely extend beyond
        # the top and bottom of the layer. To deal with such faces, we will offset the orig
        # and offset faces in opposite directions to create a tool that is slightly bigger

        offsetFace = self.FaceListOffset[0].Face
        offsetOffsetFace = OffsetFaceMoreAlongDirection(offsetFace, 1, Tolerance)
        origFace = self.FaceListOrig[0].Face
        offsetOrigFace = OffsetFaceMoreAlongDirection(origFace, 1, Tolerance)

        ProjectionEdges_a = ProjectEdgesOntoFace(edge_edges, offsetOffsetFace)
        ProjectionEdges_b = ProjectEdgesOntoFace(edge_edges, offsetOrigFace)

        # Generate faces connecting original and offset projected edges.
        # We will use this as a tool to do the cut.

        # For the moment assume only one edge

        build = BRep_Builder()  # !!!*** Are build and Perimeter still necessary????
        Perimeter = TopoDS_Compound()
        build.MakeCompound(Perimeter)

        wire_a = TopoDS_Wire()
        build.MakeWire(wire_a)
        wire_b = TopoDS_Wire()
        build.MakeWire(wire_b)

        for edgecnt in range(len(edge_edges)):
            projectionedge_a = ProjectionEdges_a[edgecnt]
            projectionedge_b = ProjectionEdges_b[edgecnt]

            build.Add(wire_a, projectionedge_a)
            build.Add(wire_b, projectionedge_b)
            pass

        # Generate cutting tool face
        CutToolGenerator = BRepOffsetAPI.BRepOffsetAPI_ThruSections()
        CutToolGenerator.AddWire(wire_a)
        CutToolGenerator.AddWire(wire_b)
        CutToolGenerator.Build()

        if (not CutToolGenerator.IsDone()):
            raise ValueError("CutTool generation failed\n")

        CutToolShape = CutToolGenerator.Shape()

        build.Add(Perimeter, CutToolShape)

        GASplitter = GEOMAlgo_Splitter()
        GASplitter.AddArgument(self.Shape)
        GASplitter.AddTool(CutToolShape)
        GASplitter.Perform()

        # if (not GASplitter.IsDone()):
        #    raise ValueError("Splitting face failed\n")

        SplitBodies = GASplitter.Shape()

        #step_writer2=STEPControl_Writer()
        #step_writer2.Transfer(CutToolShape,STEPControl_ShellBasedSurfaceModel,True)
        #step_writer2.Transfer(self.Shape, STEPControl_ManifoldSolidBrep, True)
        #step_writer2.Transfer(SplitBodies,STEPControl_ManifoldSolidBrep,True)
        #step_writer2.Write("../data/allShapes.STEP")

        #sys.modules["__main__"].__dict__.update(globals())
        #sys.modules["__main__"].__dict__.update(locals())
        #raise ValueError("Break")

        # !!!*** Need to create two layerbodies, replacing the existing layerbody in the
        # layer structure. Need to generate and sort LayerBodyFaces into all of the right places.


        # Extract the bodies created after the split operation
        bodyIterator = TopoDS_Iterator(SplitBodies)
        bodyCount = 0
        SplitBodyList=[]
        while bodyIterator.More():
            # Only consider edges
            if bodyIterator.Value().ShapeType() == TopAbs_SOLID:
                bodyCount += 1
                SplitBodyList.append(bodyIterator.Value())
                pass
            bodyIterator.Next()
            pass
        # For each body extracted, find the faces and add them to the corresponding list
        NewLayerBodies=[]
        for (bodyNum, body) in enumerate(SplitBodyList):

            NewLayerBody = LayerBody(Name="%s_Split%d" % (self.Name,bodyNum+1),
                                     Shape=body)

            # Extract faces and sort them into NewLayerBody.FaceListOrig, NewLayerBody.FaceListOffset, and NewLayerBody.FaceListSide

            # Iterate over all faces
            FaceExp = TopExp_Explorer(body, TopAbs_FACE)
            while FaceExp.More():
                # Extract the Surface object (geometry, not topology) underlying this face
                tds_Face = topods_Face(FaceExp.Current())
                FaceSurf = BRep_Tool.Surface(tds_Face)

                # WARNING: This surface matching trick seems to be rather fragile
                # ... in particular we have seen these surfaces never match
                # when generated from a cutting operation. The symptom is you
                # get only "side faces" out. A better solution is to use
                # Face.IsSame() method, but that won't work here because the faces have been split
                
                # Search for this face in the mold
                MatchedInOrig = False
                for OrigFace in self.FaceListOrig:  # Iterate over LayerBodyFaces in pre-split LayerBody
                    OrigFaceSurf = BRep_Tool.Surface(OrigFace.Face)
                    if OrigFaceSurf == FaceSurf:  # Same underlying surface
                        if MatchedInOrig:
                            raise ValueError("Same surface matched twice in mold (!?)")
                        MatchedInOrig = True
                        pass
                    pass

                # Search for this face in the offset surface
                MatchedInOffset = False
                for OffsetFace in self.FaceListOffset:  # Iterate over LayerBodyFaces in pre-split LayerBody
                    OffsetFaceSurf = BRep_Tool.Surface(OffsetFace.Face)
                    if OffsetFaceSurf == FaceSurf:  # Same underlying surface
                        if MatchedInOffset:
                            raise ValueError("Same surface matched twice in offset (!?)")
                        MatchedInOffset = True
                        pass
                    pass

                # Since these faces are extracted from the generated
                # closed object, that the normals generated in FromOCC()
                # should be outward-facing and thus we don't need to
                # provide an IsPointingInside lambda.

                if MatchedInOrig and not (MatchedInOffset):
                    # Create LayerBodyFace
                    NewLayerBody.FaceListOrig.append(LayerBodyFace.FromOCC(tds_Face, "ORIG", Owner=NewLayerBody))
                    pass
                elif MatchedInOffset and not (MatchedInOrig):
                    # Create LayerBodyFace
                    NewLayerBody.FaceListOffset.append(LayerBodyFace.FromOCC(tds_Face, "OFFSET", Owner=NewLayerBody))
                    pass
                elif MatchedInOffset and MatchedInOrig:
                    raise ValueError("Same surface matched in both offset and orig (!?)")
                else:
                    # Must be a side face
                    # Create LayerBodyFace
                    NewLayerBody.FaceListSide.append(LayerBodyFace.FromOCC(tds_Face, "OFFSET", Owner=NewLayerBody))
                    pass
                FaceExp.Next()
                pass
            NewLayerBodies.append(NewLayerBody)
            pass
        return NewLayerBodies

    def GetOffsetEdge(self):
        # return an edge along the offset direction (between ORIG and OFFSET faces)
        # return (point, tangent)

        # If there are more than 1 side faces, then perform an intersection
        # if not, go through the face edges and choose the edge that is along the offset direction
        if (len(self.FaceListSide) > 1):
            # Choose the first side face and intersect it with any other adjacent face.
            # Find the intersection edge
            # Evaluate midpoint and tangent
            sideFace1 = self.FaceListSide[0].Face

            for side in range(1,len(self.FaceListSide)):
                sideFace2 =  self.FaceListSide[side].Face

                section = BRepAlgoAPI_Section(sideFace1, sideFace2)
                section.Build()

                # sectionCompoundEdges is a compound shape
                # Need to iterate and get all edges

                if section.IsDone():
                    sectionCompoundEdges = section.Shape()
                else:
                    raise ValueError("Could not compute Section!")

                exp = TopExp_Explorer(sectionCompoundEdges, TopAbs_EDGE)

                # Iterate over all edges
                # return a list of edges

                edge_list = []
                while exp.More():
                    current_edge = topods_Edge(exp.Current())
                    edge_list.append(current_edge)
                    exp.Next()
                    pass

                if (len(edge_list) == 0):
                    continue

                assert(len(edge_list) == 1)

                intersection_edge = edge_list[0]
                (intersection_curve, start, end) = BRep_Tool().Curve(intersection_edge)  # Handle_Geom_Curve
                #point = intersection_curve.GetObject().Value((start + end)*0.5) # gp_Pnt
                point = gp_Pnt()
                tangent = gp_Vec()
                intersection_curve.GetObject().D1((start + end)*0.5, point, tangent)
                break
        else:
            # Only one side face
            EdgeExplorer = TopExp_Explorer(self.FaceListSide[0].Face, TopAbs_EDGE)
            # Iterate over all faces, creating LayerBodyFaces
            while EdgeExplorer.More():
                edge = topods_Edge(EdgeExplorer.Current())
                VertexExplorer = TopExp_Explorer(edge, TopAbs_VERTEX)
                vertex1 = topods_Vertex(VertexExplorer.Current())
                VertexExplorer.Next()
                vertex2 = topods_Vertex(VertexExplorer.Current())

                if (not vertex1.IsSame(vertex2)):
                    (intersection_curve, start, end) = BRep_Tool().Curve(edge)  # Handle_Geom_Curve
                    # point = intersection_curve.GetObject().Value((start + end)*0.5) # gp_Pnt
                    point = gp_Pnt()
                    tangent = gp_Vec()
                    intersection_curve.GetObject().D1((start + end) * 0.5, point, tangent)
                    break

                EdgeExplorer.Next()
            pass

        return (np.array((point.X(),point.Y(),point.Z()),dtype='d'),np.array((tangent.X(),tangent.Y(),tangent.Z()),dtype='d'))

    # !!!*** Need API to provide points for a wire segment to do the splitting,
    # Need to connect the wire segment to the domain boundary
    # Need to bond the region between segment and boundary




    pass


# Defined an implicit interface where
# LayerBodies and LayerMolds
# all have .Name members (string) NOTE: Removed from LayerMolds because we don't expect
#                                       to need to save LayerMolds. Probably we want a
#                                       different class for shells, which would have
#                                       these characteristics. 
# and have .Shape members (TopoDS Shape)
# and StepModelType members with the appropriate data type for STEPControl_Writer()


class LayerMold(object):
    """ The LayerMold represents the collection of faces 
    from which we will generate an offset surface and then side faces, 
    and then construct a single layer. 
    
    The collection of faces is represented as a TopoDS_Shape which is a
    TopoDS_Shell

    Value semantics: The LayerMold is immutable once created. If 
    there is a need to change it, create a new LayerMold being sure 
    to deep copy anything being changed

"""

    # Owner = None   # reference to owning layer... DON'T THINK THIS SHOULD BE NECESSARY ANYMORE. MANY/ALL LAYERMOLDS ARE NOT OWNED BY A LAYER 
    Direction = None # Direction, i.e. "ORIG" or "OFFSET". if "ORIG" then growth direction should be flipped from the face normals
                     # when offsetting in the "OFFSET" direction,
                     # if "OFFSET" then growth direction should be flipped if offsetting in the "ORIG" direction
    Shape = None # OpenCascade TopoDS Shape which is a Shell (?) which consists of multiple Faces
    FaceList = None # List of LayerBodyFace
    
    # Constant:
    StepModelType = STEPControl_ShellBasedSurfaceModel # Used for OpenCascade STEP writer
    
    
    
    def __init__(self,**kwargs):
        self.FaceList=[]
        self.Direction="NODIR"

        for argname in kwargs:
            if not hasattr(self,argname):
                raise ValueError("Invalid attribute")
            setattr(self,argname,kwargs[argname])
            pass
        pass

    @classmethod
    def FromShell(cls,Shell,OrigDirPoint=None,OrigDirNormal=None,OrigDirTolerance=None):
        """Create a LayerMold from a TopAbs_Shell. If OrigDirPoint and normal are 
        None, then it assumes the Shell is already oriented with the desired normal. 
        Otherwise will flip the normal as described for FromFile()."""
        Direction = "ORIG"

        if OrigDirPoint is not None:
            OrigDirPointVertex=BRepBuilderAPI.BRepBuilderAPI_MakeVertex(gp_Pnt(OrigDirPoint[0],OrigDirPoint[1],OrigDirPoint[2])).Vertex()
            pass
        
        FaceExp=TopExp_Explorer(Shell,TopAbs_FACE)
        # Iterate over all faces, creating LayerBodyFaces
        FaceList=[]
        ClosestDist=np.Inf
        ClosestNormal=None
        while FaceExp.More():
            
            FaceList.append(LayerBodyFace.FromOCC(topods_Face(FaceExp.Current()),"NODIR"))

            if OrigDirPoint is not None:
                # We will have to flip Direction if nearest point on the object to
                # OrigDirPoint has outward normal component away from the
                # OrigDirNormal direction
                
                # This is an example of how to evaluate point-to-object
                # distance and find (U,V) coordinates.

                # Find the closest point by this method: 
                # https://www.opencascade.com/content/closest-point-step-object
                DistCalc=BRepExtrema_DistShapeShape(FaceExp.Current(),OrigDirPointVertex)
                DistCalc.Perform()
                ThisDist=DistCalc.Value()
                
                if DistCalc.NbSolution() > 0 and ThisDist < ClosestDist:
                    ClosestDist=ThisDist

                    DistSupport = DistCalc.SupportOnShape1(1).ShapeType()
                    if DistSupport != TopAbs_FACE:
                        # This test face found nearest point to be
                        # an edge or a vertex, not inside the face
                        # itself
                        continue
                        #raise ValueError("Nearest point on shell to selected point %s for determining ORIGINAL direction is a "
                        #                 "vertex or and edge, not an interior point on the face. Select a different OrigDirPoint." % (str(OrigDirPoint)))
                    
                    (ClosestU,ClosestV)=DistCalc.ParOnFaceS1(1) # Evaluate (u,v) coordinates on this face of closest point
                    
                    # Here is an example of how to extract
                    # a normal vector given (U,V) coordinates
                    ThisFaceSurface=BRep_Tool().Surface(topods_Face(FaceExp.Current()))
                    ClosestUVProps = GeomLProp_SLProps(ThisFaceSurface,ClosestU,ClosestV,1,OrigDirTolerance)
                    ClosestNormalVec=gp_Vec(ClosestUVProps.Normal())
                    ClosestNormal=np.array((ClosestNormalVec.X(),ClosestNormalVec.Y(),ClosestNormalVec.Z()),dtype='d')
                    if FaceExp.Current().Orientation==TopAbs_REVERSED:
                        # Face is reversed from underlying surface -> we need to flip the normal
                        ClosestNormal=-ClosestNormal
                        pass
                    pass
                pass

            FaceExp.Next()
            
            pass

        if OrigDirPoint is not None:
            if ClosestNormal is None:
                raise ValueError("No closest point on loaded shape to vertex at %s. Try passing OrigDirPoint and OrigDirNormal parameters to specify the point where the ORIGINAL direction is defined." % (str(OrigDirPoint)))
            
            
            NormalDir = ClosestNormal
            if np.dot(NormalDir,OrigDirNormal) < 0.0:
                # Normal is backwards from the direction our caller requested,
                # flip it, by switching Direction from "ORIG" to "OFFSET"
                Direction="OFFSET"
                pass
            pass
        
        # Instantiate the class, calling constructor
        return cls(Direction=Direction,
                   Shape=Shell,
                   FaceList=FaceList)

    @classmethod
    def CutMoldFromShell(cls, shellfilename, toolfilename, OrigDirPoint=np.array((0.0, 0.0, 0.0)), OrigDirNormal=np.array((0.0, 0.0, -1.0)),
                 OrigDirTolerance=1e-6):
        """Create a LayerMold from a STEP, IGES, or BREP file and cut the mold using the tool
        The LayerMold containing a surface. OrigDirPoint and OrigDirNormal
        are used to define the "ORIG" direction. The closest point
        on the surface to OrigDirPoint is found. The ORIG direction
        is the side of the surface with positive component in the
        OrigDirNormal direction.
        """

        ShellShapes = loaders.load_byfilename(shellfilename)
        ToolShapes = loaders.load_byfilename(toolfilename)

        ShellExp = TopExp_Explorer(ShellShapes, TopAbs_SHELL)
        # Iterate over all shells
        ShellList = []
        while ShellExp.More():
            CurrentShellShape = ShellExp.Current()
            if CurrentShellShape.ShapeType() != TopAbs_SHELL:
                raise ValueError("File %s type is %s, not TopAbs_SHELL" % (shellfilename, str(CurrentShellShape.ShapeType)))
            ShellList.append(CurrentShellShape)
            ShellExp.Next()

        if (len(ShellList) > 1):
            raise ValueError("File %s contains more than one TopAbs_SHELL" % (shellfilename))
        ShellShape = ShellList[0]

        ToolExp = TopExp_Explorer(ToolShapes, TopAbs_SHELL)
        # Iterate over all shells
        ToolList = []
        while ToolExp.More():
            CurrentToolShape = ToolExp.Current()
            if CurrentToolShape.ShapeType() != TopAbs_SHELL:
                raise ValueError("File %s type is %s, not TopAbs_SHELL" % (shellfilename, str(CurrentToolShape.ShapeType)))
            ToolList.append(CurrentToolShape)
            ToolExp.Next()

        if (len(ToolList) > 1):
            raise ValueError("File %s contains more than one TopAbs_SHELL" % (shellfilename))

        ToolShape = ToolList[0]

        GASplitter = GEOMAlgo_Splitter()
        GASplitter.AddArgument(ShellShape)
        GASplitter.AddTool(ToolShape)
        GASplitter.Perform()

        SplitShells = GASplitter.Shape()

        SplitFaceExp = TopExp_Explorer(SplitShells, TopAbs_FACE)
        # Iterate over all shells
        SplitFaceList = []
        while SplitFaceExp.More():
            CurrentShellShape = SplitFaceExp.Current()
            SplitFaceList.append(CurrentShellShape)
            SplitFaceExp.Next()

        GProps1 = GProp_GProps()
        brepgprop_SurfaceProperties(SplitFaceList[0], GProps1)
        SurfArea1 = GProps1.Mass()

        GProps2 = GProp_GProps()
        brepgprop_SurfaceProperties(SplitFaceList[1], GProps2)
        SurfArea2 = GProps2.Mass()

        if SurfArea1 > SurfArea2:
            ShellModel = SplitFaceList[0]
            MoldShape = SplitFaceList[1]
            pass
        else:
            ShellModel = SplitFaceList[1]
            MoldShape = SplitFaceList[0]
            pass

        Name = os.path.splitext(os.path.split(shellfilename)[1])[0]  # Based on filename with no extension
        # Owner = None

        # Create a compound face
        # compoundBuilder = BRep_Builder()
        # MoldCompound = TopoDS_Compound()
        # compoundBuilder.MakeCompound(MoldCompound)
        # compoundBuilder.Add(MoldCompound, MoldShape);

        # Create a shell from the face
        shellBuilder = BRep_Builder()
        MoldShell = TopoDS_Shell()
        shellBuilder.MakeShell(MoldShell)
        shellBuilder.Add(MoldShell, MoldShape)

        # step_writer=STEPControl_Writer()
        # step_writer.Transfer(MoldShell,STEPControl_ShellBasedSurfaceModel,True)
        # #step_writer.Transfer(ShellModel,STEPControl_ShellBasedSurfaceModel,True)
        # step_writer.Write("../data/SplitShell.STEP")
        #
        # sys.modules["__main__"].__dict__.update(globals())
        # sys.modules["__main__"].__dict__.update(locals())
        # raise ValueError("Break")

        return (cls.FromShell(MoldShell, OrigDirPoint, OrigDirNormal, OrigDirTolerance),
                ShellModel)

    @classmethod
    def CutMoldFromSolid(cls, solidfilename, toolfilename, OrigDirPoint=np.array((0.0, 0.0, 0.0)), OrigDirNormal=np.array((0.0, 0.0, -1.0)),
                         OrigDirPointTolerance=1e-5,OrigDirNormalTolerance=1e-6):
        """Create a LayerMold from a STEP, IGES, or BREP file and cut the mold using the tool
        The LayerMold containing a surface. OrigDirPoint and OrigDirNormal
        are used to define the "ORIG" direction. The closest point
        on the surface to OrigDirPoint is found. The ORIG direction
        is the side of the surface with positive component in the
        OrigDirNormal direction.
        """

        SolidShaoes = loaders.load_byfilename(solidfilename)
        ToolShapes = loaders.load_byfilename(toolfilename)

        SolidExp = TopExp_Explorer(SolidShaoes, TopAbs_SOLID)
        # Iterate over all solids
        SolidList = []
        while SolidExp.More():
            CurrentSolidShape = SolidExp.Current()
            if CurrentSolidShape.ShapeType() != TopAbs_SOLID:
                raise ValueError("File %s type is %s, not TopAbs_SOLID" % (solidfilename, str(CurrentSolidShape.ShapeType)))
            SolidList.append(CurrentSolidShape)
            SolidExp.Next()

        if (len(SolidList) > 1):
            raise ValueError("File %s contains more than one TopAbs_SOLID" % (solidfilename))
        OrigSolidShape = SolidList[0]

        ToolExp = TopExp_Explorer(ToolShapes, TopAbs_SHELL)
        # Iterate over all shells
        ToolList = []
        while ToolExp.More():
            CurrentToolShape = ToolExp.Current()
            if CurrentToolShape.ShapeType() != TopAbs_SHELL:
                raise ValueError("File %s type is %s, not TopAbs_SHELL" % (shellfilename, str(CurrentToolShape.ShapeType)))
            ToolList.append(CurrentToolShape)
            ToolExp.Next()

        if (len(ToolList) > 1):
            raise ValueError("File %s contains more than one TopAbs_SHELL" % (shellfilename))

        ToolShape = ToolList[0]

        GASplitter = GEOMAlgo_Splitter()
        GASplitter.AddArgument(OrigSolidShape)
        GASplitter.AddTool(ToolShape)
        GASplitter.Perform()

        SplitSolids = GASplitter.Shape()

        SplitSolidExp = TopExp_Explorer(SplitSolids, TopAbs_SOLID)
        # Iterate over all solids
        SplitSolidList = []
        while SplitSolidExp.More():
            CurrentSolidShape = SplitSolidExp.Current()
            SplitSolidList.append(CurrentSolidShape)
            SplitSolidExp.Next()
            pass

        GProps1 = GProp_GProps()
        brepgprop_VolumeProperties(SplitSolidList[0], GProps1)
        SurfArea1 = GProps1.Mass()

        GProps2 = GProp_GProps()
        brepgprop_VolumeProperties(SplitSolidList[1], GProps2)
        SurfArea2 = GProps2.Mass()

        if SurfArea1 > SurfArea2:
            HollowedSolidModel = SplitSolidList[0]
            LayerSolidModel = SplitSolidList[1]
            pass
        else:
            HollowedSolidModel = SplitSolidList[1]
            LayerSolidModel = SplitSolidList[0]
            pass

        # step_writer=STEPControl_Writer()
        # # step_writer.Transfer(SolidModel,STEPControl_ManifoldSolidBrep,True)
        # step_writer.Transfer(LayerSolidModel,STEPControl_ManifoldSolidBrep,True)
        # step_writer.Write("../data/SplitSolid.STEP")


        # Iterate through all the faces and identify the ORIG face based on the input
        # point and normal. Extract the face and send as input to create the layer
        # structure

        # ... So really we should iterate through the side faces of the resulting cut pieces
        # by identifying which faces exist in both pieces. Then that set of side faces
        # topologically splits the ORIG and the OFFSET sides... At which point we can
        # identify the ORIG side by one of its faces matching OrigDirPoint and OrigDirNormal.
        #
        # For the time being we assume that the ORIG and OFFSET sides are singles faces and we
        # can identify the ORIG directly by using OrigDirPoint and OrigDirNormal. 

        OrigFace = SelectFaceByPointNormal(LayerSolidModel,OrigDirPoint,OrigDirNormal,OrigDirPointTolerance,OrigDirNormalTolerance)
        if OrigFace is None:
            raise ValueError("Point %s not found any surface of LayerSolidModel cut from %s by %s" % (str(OrigDirPoint),solidfilename, toolfilename))


        # Create a shell from the face
        shellBuilder = BRep_Builder()
        OrigShell = TopoDS_Shell()
        shellBuilder.MakeShell(OrigShell)
        shellBuilder.Add(OrigShell, OrigFace)

        # step_writer=STEPControl_Writer()
        # step_writer.Transfer(MoldShell,STEPControl_ShellBasedSurfaceModel,True)
        # #step_writer.Transfer(ShellModel,STEPControl_ShellBasedSurfaceModel,True)
        # step_writer.Write("../data/SplitShell.STEP")
        #
        # sys.modules["__main__"].__dict__.update(globals())
        # sys.modules["__main__"].__dict__.update(locals())
        # raise ValueError("Break")

        return (cls.FromShell(OrigShell, OrigDirPoint, OrigDirNormal, OrigDirPointTolerance),
                solid.Solid.FromOCC(os.path.split(solidfilename)[1],OrigSolidShape,PointTolerance=OrigDirPointTolerance,NormalTolerance=OrigDirNormalTolerance))
    
    @classmethod
    def FromFile(cls,filename,OrigDirPoint=np.array((0.0,0.0,0.0)),OrigDirNormal=np.array((0.0,0.0,1.0)),OrigDirTolerance=1e-6):
        """Create a LayerMold from a STEP, IGES, or BREP file 
        containing a surface. OrigDirPoint and OrigDirNormal 
        are used to define the "ORIG" direction. The closest point
        on the surface to OrigDirPoint is found. The ORIG direction 
        is the side of the surface with positive component in the
        OrigDirNormal direction. 
        """


        Shape = loaders.load_byfilename(filename)

        if Shape.ShapeType() != TopAbs_SHELL:
            raise ValueError("File %s type is %s, not TopAbs_SHELL" % (filename,str(Shape.ShapeType)))
        
        Name = os.path.splitext(os.path.split(filename)[1])[0] # Based on filename with no extension
        #Owner = None


        return cls.FromShell(topods_Shell(Shape),OrigDirPoint,OrigDirNormal,OrigDirTolerance)
    
    @classmethod
    def FromFaceLists(cls,FaceLists):
        """Build a LayerMold from a collection of FaceLists, containing faces that can be attached
        and that are presumed to already be oriented correctly and consistently for the 
        "ORIG" direction of the LayerMold"""

        # Per https://www.opencascade.com/content/how-create-shell since our faces
        # should already share edges and vertices, we should be able to use BRep_Builder.
        # ... if this fails then could maybe try BRepOffsetAPI_Sewing

        Builder=BRep_Builder();
        Shell=TopoDS_Shell();
        Builder.MakeShell(Shell);
        for FaceList in FaceLists:
            for Face in FaceList:  # Face is a LayerBodyFace
                Builder.Add(Shell, Face.Face)
                pass
            pass

        #!!!!****
        # Use BRepFeat_Gluer https://www.opencascade.com/doc/occt-6.9.1/refman/html/class_b_rep_feat___gluer.html
        # iterating over pairs of faces in the face list, finding shared wire edges,
        # gluing those two components, then repeat with those two merged,
        # hopfully accommodating any components that won't merge.
        
        # OCC.BRepFeat.BRepFeat_Gluer
        
        ## Remove unnecessary internal edges
        #Fuser = BRepLib_FuseEdges(Shell);
        #Fuser.Perform();

        #BRepFeat.BRepFeat_Gluer()
        
        return cls.FromShell(Shell)
        #return cls.FromShell(Fuser.Shape())
    
    pass


class LayerBodyFace(object): # Formerly LayerSurface
    """The LayerBodyFace represents a TopoDS_Face, one of (potentially) 
    many representing the faces on a LayerBody or LayerMold

    Multiple LayerBodyFace objects may exist representing the same
    actual face, or subregions of such a face. There is no guarantee
    of uniqueness EXCEPT THAT within a particular layer faces may 
    only exist once, and within multiple layers being assembled, 
    faces may only exist twice (Additional LayerBodies may exist 
    but they shouldn't be part of layers). 

    Value semantics: The LayerBodyFace is immutable once created. If 
    there is a need to change it, create a new LayerBodyFace being sure 
    to deep copy anything being changed
    
    """
    Face = None  # TopoDS_Face (not unique; other TopoDS_Faces may exist that also represent this face) 
    Point = None # Numpy Array representing coordinates of a point on the face
    Normal = None # Numpy array representing unit normal pointing away from the
                  # layerbody at the location of Point. Direction is arbitrary if this face is not
                  # part of a LayerBody
    ParPoint = None  # Numpy Array representing parametric coordinates of the above point on the face
                     # This is initialized to (0.05,0.05) and passed to OCCFindPointNormal, where
                     # it is updated if the reference parameter is outside the trimmed region
    # There had been an "InitialSurface" boolean... I think this info
    # should be covered by Direction=="ORIG"
    Direction = None # "ORIG", "OFFSET", "SIDE", or "NODIR", representing the side of the Owning LayerBody which corresponds to this Face

    Owner = None # If this LayerBodyFace may be part of a LayerBody, this is the LayerBody of which it might be a part.
    # ****!!!! NOTE: as of 2/8/19, not always updated when we do delamination splits!!!***
    
    BCType = None # Formerly DelaminationType: None, "NODELAM" "NOMODEL", "COHESIVE", "CONTACT" or, "TIE"
    # MatchingFace = None # Formerly SurfPair, This would be the matching LayerBodyFace in the adjacent (or non-adjacent)
    # layer or stiffener, assigned by adjacent_layer_boundary_conditions() ... To maintain immutability of these objects suggest that
    # instead of this attribute, have adjacent_layer_boundary_conditions() create a dictionary
    # that can be used to look up matching faces. 

    # StiffenerGenerated and StiffnerPaired members removed (should no longer be necessary)
    # CreatedFrom = None   # (probably shouldn't be necessary anymore) but would be the origin LayerBodyFace


    def __init__(self,**kwargs):
        self.Direction="NODIR"

        for argname in kwargs:
            if not hasattr(self,argname):
                raise ValueError("Invalid attribute")
            setattr(self,argname,kwargs[argname])
            pass
        pass

    def __eq__(self,other):
        # Test equality of faces. This checks whether they the faces
        # share the same underlying geometry (surface object)
        #
        # Ideally, it would check for geometric equivalence...
        # but OpenCascade has no function for this. 
        #
        # Since the layer bodies, etc should be built from the same
        # underlying geometry we should in general be OK.
        # We do have to be careful, for example, that if we do an
        # imprint operation on a surface that we don't do a second
        # imprint operation on the mating surface... because the
        # 2nd imprint operation would give us another copy of the
        # geometry that would then not show pointer equivalence
        
        #ThisSurface = BRep_Tool().Surface(self.Face)
        #OtherSurface = BRep_Tool().Surface(other.Face)
        #print("ThisSurface=%s" % (str(ThisSurface)))
        #print("OtherSurface=%s" % (str(OtherSurface)))
        
        #print("Equality operator returns %s" % (str(ThisSurface.DumpToString() == OtherSurface.DumpToString())))
        
        #return ThisSurface.DumpToString() == OtherSurface.DumpToString() and self.Face.IsEqual(other.Face)


        # Update: Use OpenCascade IsSame operator
        # That tests for same underlying TShape (?) with same location but not necessarily orientation ... normal may be flipped.

        # Note that this also works for testing equality
        # with ImmutableSolidFace objects. 
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
    # Returning a constant here means that dictionaries indexed by LayerBodyFace
    # will perform like lists (but at least you will be able to do direct lookups)
    def __hash__(self):
        return 0
    
    @classmethod
    def FromOCC(cls,Face,Direction,Owner=None,BCType=None,IsPointingInside=None,OrigPointTolerance=1e-5,OrigNormalTolerance=1e-6):
        """ Create a LayerBodyFace from a TopoDS_Face. 
        IsPointingInside should be None (if this LayerBodyFace is 
        not actually part of a LayerBody, or if the OCC-Evaluated normal
        is guaranteed to be outward-facing) or a function/method/lambda 
        that, given a Normal, returns whether  that normal is pointing inside
        the LayerBody."""

        (Point,Normal,ParPoint)=FindOCCPointNormal(Face,PointTolerance=OrigPointTolerance,NormalTolerance=OrigNormalTolerance)

        # For a LayerBodyFace inside a LayerBody, the normal
        # should be pointing outward. 
        if IsPointingInside is not None and IsPointingInside(Normal):
            Normal=-Normal
            pass
        
        return cls(Face=Face,
                   Point=Point,
                   Normal=Normal,
                   ParPoint=ParPoint,
                   Direction=Direction,
                   Owner=Owner,
                   BCType=BCType)
        

    pass

    
if __name__=="__main__":

    pointTolerance = 1e-6
    Mold = LayerMold.FromFile(os.path.join("..","data","FlatMold3.STEP"))
    Layer1=Layer.CreateFromMold("Layer 1",Mold,2.0,"OFFSET",pointTolerance)
    Layer2=Layer.CreateFromMold("Layer 2",Layer1.OffsetMold(),2.0,"OFFSET",pointTolerance)
    Layer3=Layer.CreateFromMold("Layer 3",Layer2.OffsetMold(),2.0,"OFFSET",pointTolerance)

    Layer2.Split(os.path.join("..","data","SplitLine.csv"), pointTolerance)

    step_writer=STEPControl_Writer()
    
    step_writer.Transfer(Layer1.BodyList[0].Shape,STEPControl_ManifoldSolidBrep,True)
    step_writer.Transfer(Layer2.BodyList[0].Shape,STEPControl_ManifoldSolidBrep,True)
    step_writer.Transfer(Layer2.BodyList[1].Shape,STEPControl_ManifoldSolidBrep,True)
    step_writer.Transfer(Layer3.BodyList[0].Shape,STEPControl_ManifoldSolidBrep,True)
    step_writer.Write("../data/Layers.step")

    pass
