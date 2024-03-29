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
import copy
import time
import os.path
import csv
import numpy as np
import OCC

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
    from OCC.Core.TopoDS import TopoDS_Compound
    from OCC.Core.TopoDS import TopoDS_Wire
    from OCC.Core.TopoDS import TopoDS_Vertex
    from OCC.Core.TopoDS import TopoDS_Shell
    from OCC.Core.TopoDS import topods_Shell
    from OCC.Core.TopoDS import topods_Face
    from OCC.Core.TopoDS import topods_Edge
    from OCC.Core.TopoDS import topods_Wire
    from OCC.Core.GC import GC_MakeSegment
    from OCC.Core.BRep import BRep_Builder
    from OCC.Core.BRep import BRep_Tool
    from OCC.Core.BRepExtrema import BRepExtrema_DistShapeShape
    from OCC.Core import BRepLib
    from OCC.Core import BRepOffsetAPI
    from OCC.Core import BRepOffset
    from OCC.Core.Geom import Geom_OffsetCurve
    from OCC.Core.Geom import Geom_Curve
    from OCC.Core.Geom import Geom_Surface
    from OCC.Core import BRepBuilderAPI
    from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge
    from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeWire
    from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
    #from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_FaceError
    #from OCC.Core.BRepClass import BRepClass_FacePassiveClassifier
    from OCC.Core.BRepClass import BRepClass_FaceExplorer
    from OCC.Core.BRepClass import BRepClass_FClassifier
    from OCC.Core.ShapeAnalysis import ShapeAnalysis_FreeBoundsProperties
    from OCC.Core.BRepTools import breptools_Read
    from OCC.Core.BRepTools import breptools_Write
    from OCC.Core.TopExp import TopExp_Explorer
    from OCC.Core.TopAbs import TopAbs_ON
    from OCC.Core.TopAbs import TopAbs_IN
    from OCC.Core.TopAbs import TopAbs_OUT
    from OCC.Core.TopAbs import TopAbs_FACE
    from OCC.Core.TopAbs import TopAbs_VERTEX
    from OCC.Core.TopAbs import TopAbs_EDGE
    from OCC.Core.TopAbs import TopAbs_SHELL
    from OCC.Core.TopAbs import TopAbs_FORWARD
    from OCC.Core.TopAbs import TopAbs_REVERSED
    from OCC.Core.GeomAbs import GeomAbs_Arc
    from OCC.Core.GeomAbs import GeomAbs_Intersection
    from OCC.Core.Geom import Geom_Line
    from OCC.Core.TopTools import TopTools_ListIteratorOfListOfShape
    from OCC.Core.TopTools import TopTools_ListOfShape
    from OCC.Core.GeomLProp import GeomLProp_SLProps
    from OCC.Core.GProp import GProp_GProps
    from OCC.Core.BRepGProp import brepgprop_SurfaceProperties
    from OCC.Core.gp import gp_Pnt2d
    from OCC.Core.gp import gp_Vec
    from OCC.Core.gp import gp_Dir
    from OCC.Core.gp import gp_Pnt
    from OCC.Core.gp import gp_Pln

    use_BOPAlgo=False
    try:
        import OCC.Core.GEOMAlgo
        pass
    except ImportError:
        use_BOPAlgo=True
        pass

    if use_BOPAlgo: 
        from OCC.Core.BOPAlgo import BOPAlgo_Splitter as GEOMAlgo_Splitter
        pass
    else:
        from OCC.Core.GEOMAlgo import GEOMAlgo_Splitter
        pass
    from OCC.Core.Geom2d import Geom2d_Curve
    from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Fuse
    from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Section
    from OCC.Core import GeomProjLib
    from OCC.Core.TColgp import TColgp_Array1OfPnt
    from OCC.Core.TColgp import TColgp_HArray1OfPnt
    from OCC.Core.GeomAPI import (GeomAPI_Interpolate, GeomAPI_PointsToBSpline)
    from OCC.Core.ShapeFix import ShapeFix_Edge
    from OCC.Core.ShapeAnalysis import ShapeAnalysis_WireOrder

    from OCC.Core.STEPControl import STEPControl_Reader
    from OCC.Core.STEPControl import STEPControl_Writer
    from OCC.Core.STEPControl import STEPControl_ManifoldSolidBrep
    from OCC.Core.STEPControl import STEPControl_Writer,STEPControl_ShellBasedSurfaceModel,STEPControl_GeometricCurveSet

    pass
else: # !use_OCC_core
    from OCC.TopoDS import topods
    from OCC.TopoDS import TopoDS_Face
    from OCC.TopoDS import TopoDS_Shape
    from OCC.TopoDS import TopoDS_Compound
    from OCC.TopoDS import TopoDS_Wire
    from OCC.TopoDS import TopoDS_Vertex
    from OCC.TopoDS import TopoDS_Shell
    from OCC.TopoDS import topods_Shell
    from OCC.TopoDS import topods_Face
    from OCC.TopoDS import topods_Edge
    from OCC.TopoDS import topods_Wire
    from OCC.GC import GC_MakeSegment
    from OCC.BRep import BRep_Builder
    from OCC.BRep import BRep_Tool
    from OCC.BRepExtrema import BRepExtrema_DistShapeShape
    from OCC import BRepLib
    from OCC import BRepOffsetAPI
    from OCC import BRepOffset
    from OCC.Geom import Geom_OffsetCurve
    from OCC.Geom import Geom_Curve
    from OCC.Geom import Geom_Surface
    from OCC import BRepBuilderAPI
    from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeEdge
    from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeWire
    from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeFace
    #from OCC.BRepBuilderAPI import BRepBuilderAPI_FaceError
    #from OCC.BRepClass import BRepClass_FacePassiveClassifier
    from OCC.BRepClass import BRepClass_FaceExplorer
    from OCC.BRepClass import BRepClass_FClassifier
    from OCC.ShapeAnalysis import ShapeAnalysis_FreeBoundsProperties
    from OCC.BRepTools import breptools_Read
    from OCC.BRepTools import breptools_Write
    from OCC.TopExp import TopExp_Explorer
    from OCC.TopAbs import TopAbs_ON
    from OCC.TopAbs import TopAbs_IN
    from OCC.TopAbs import TopAbs_OUT
    from OCC.TopAbs import TopAbs_FACE
    from OCC.TopAbs import TopAbs_VERTEX
    from OCC.TopAbs import TopAbs_EDGE
    from OCC.TopAbs import TopAbs_SHELL
    from OCC.TopAbs import TopAbs_FORWARD
    from OCC.TopAbs import TopAbs_REVERSED
    from OCC.GeomAbs import GeomAbs_Arc
    from OCC.GeomAbs import GeomAbs_Intersection
    from OCC.Geom import Geom_Line
    from OCC.TopTools import TopTools_ListIteratorOfListOfShape
    from OCC.TopTools import TopTools_ListOfShape
    from OCC.GeomLProp import GeomLProp_SLProps
    from OCC.GProp import GProp_GProps
    from OCC.BRepGProp import brepgprop_SurfaceProperties
    from OCC.gp import gp_Pnt2d
    from OCC.gp import gp_Vec
    from OCC.gp import gp_Dir
    from OCC.gp import gp_Pnt
    from OCC.gp import gp_Pln
    from OCC.GEOMAlgo import GEOMAlgo_Splitter
    from OCC.Geom2d import Geom2d_Curve
    from OCC.BRepAlgoAPI import BRepAlgoAPI_Fuse
    from OCC.BRepAlgoAPI import BRepAlgoAPI_Section
    from OCC import GeomProjLib
    from OCC.TColgp import TColgp_Array1OfPnt
    from OCC.TColgp import TColgp_HArray1OfPnt
    from OCC.GeomAPI import (GeomAPI_Interpolate, GeomAPI_PointsToBSpline)
    from OCC.ShapeFix import ShapeFix_Edge
    from OCC.ShapeAnalysis import ShapeAnalysis_WireOrder

    from OCC.STEPControl import STEPControl_Reader
    from OCC.STEPControl import STEPControl_Writer
    from OCC.STEPControl import STEPControl_ManifoldSolidBrep
    from OCC.STEPControl import STEPControl_Writer,STEPControl_ShellBasedSurfaceModel,STEPControl_GeometricCurveSet
    pass


from .layer import LayerBody,LayerBodyFace
from .tools import ProjectEdgesOntoFace,FindOCCPointNormal,SelectFaceByPointNormal,OCCPointInFace

from . import loaders
from . import layer


def FaceFaceIntersect(face1, face2):

    # Adarsh: Try reimplementing with GeomAlgo_Splitter ***!!!
    section = BRepAlgoAPI_Section(face1, face2)
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
        raise ValueError("Intersection failed! Check gap width for no model zone.")


    # build = BRep_Builder()
    # Perimeter = TopoDS_Compound()
    # build.MakeCompound(Perimeter)
    # wire = TopoDS_Wire()
    # build.MakeWire(wire)
    # while exp.More():
    #     current_edge = exp.Current()
    #     build.Add(wire, current_edge)
    #     exp.Next()
    #     pass

    # WireBuilder = BRepBuilderAPI_MakeWire()
    # WireBuilder.Add(wire)
    # WireShape = WireBuilder.Shape()
    # assert (WireShape.Closed())

    # step_writer2 = STEPControl_Writer()
    # step_writer2.Transfer(face1,STEPControl_ShellBasedSurfaceModel,True)
    # step_writer2.Transfer(face2,STEPControl_ShellBasedSurfaceModel,True)
    # step_writer2.Transfer(sectionCompoundEdges,STEPControl_GeometricCurveSet,True)
    # step_writer2.Write("../data/allShapes.STEP")
    #
    # sys.modules["__main__"].__dict__.update(globals())
    # sys.modules["__main__"].__dict__.update(locals())
    # raise ValueError("Break")

    return edge_list


def CreateReferenceFace(edges, face, scale, tolerance):
    # Evaluate points on projected curve and evaluate the parametric points
    totalNumPoints = 0;
    if (len(edges) == 0):
        raise ValueError("No input edges provided")

    curveParPts = []
    for edgecnt in range(len(edges)):
        edge = edges[edgecnt]
        numEdgePoints = 250
        totalNumPoints = totalNumPoints + numEdgePoints

        # Make sure this edge has a pcurve (2d projected curve) on this face
        edgefixer = ShapeFix_Edge()
        edgefixer.FixAddPCurve(edge, face, False)

        (curve, parStart, parEnd) = BRep_Tool().CurveOnSurface(edge, face)
        #print(parStart, parEnd)

        if not hasattr(curve,"Value"):
            curve = curve.GetObject() # backward compatibility
            pass
        
        for indU in range(numEdgePoints):
            u = indU * (1.0 / numEdgePoints) * (parEnd - parStart) + parStart
            curvePoint = curve.Value(u)
            u = curvePoint.X()
            v = curvePoint.Y()
            #print(u, v)
            curveParPts.append([u*scale, v*scale, 0])
            pass

    curvePointTemp = [curveParPts[0][0], curveParPts[0][1], curveParPts[0][2]]
    curveParPts.append(curvePointTemp)

    # If there are n entries in the delam_outlist, one of which is doubled (start and end). There will be n-1 segments
    parPointsHArray = TColgp_HArray1OfPnt(1, len(curveParPts))

    for pos in range(len(curveParPts)):
        current_point = gp_Pnt(curveParPts[pos][0], curveParPts[pos][1], curveParPts[pos][2])
        parPointsHArray.SetValue(pos + 1, current_point)
        pass

    if hasattr(parPointsHArray,"GetHandle"): # Backward compatibility
        parPointsHArray=parPointsHArray.GetHandle()
        pass
    
    # Interpolate the points to make a closed curve
    interpAPI = GeomAPI_Interpolate(parPointsHArray, False, tolerance)
    interpAPI.Perform()
    if interpAPI.IsDone():
        delam_par_curve = interpAPI.Curve()
    else:
        raise ValueError("Parametric curve interpolation failed")

    # Convert a curve to edge and then to Shape
    delam_par_edge = BRepBuilderAPI_MakeEdge(delam_par_curve).Edge()
    WireBuilder = BRepBuilderAPI_MakeWire()
    WireBuilder.Add(delam_par_edge)
    ParWireShape = WireBuilder.Shape()
    assert (ParWireShape.Closed())

    P = gp_Pln()
    FaceBuilder = BRepBuilderAPI_MakeFace(P, topods_Wire(ParWireShape), True)
    referenceFace = FaceBuilder.Face()

    print ("Created reference face")

    #error = FaceBuilder.Error()
    #if (error != BRepBuilderAPI_FaceError.BRepBuilderAPI_FaceDone):
    #    print("Reference face generation failed!")
    #    pass

    # step_writer2 = STEPControl_Writer()
    # step_writer2.Transfer(referenceFace,STEPControl_ShellBasedSurfaceModel,True)
    # step_writer2.Transfer(ParWireShape, STEPControl_GeometricCurveSet, True)
    # step_writer2.Write("../data/allShapes.STEP")

    # sys.modules["__main__"].__dict__.update(globals())
    # sys.modules["__main__"].__dict__.update(locals())
    # raise ValueError("Break")

    return referenceFace


def CreateDelaminationWire(delam_outline, tolerance):
    # Create a delamination wire from the delam_outline file

    delam_outlist = []
    with open(delam_outline) as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='"')
        for row in reader:
            if len(row) != 3:
                raise ValueError("Malformed row in CSV file %s: %s" % (delam_outline, ",".join(row)))
            try:
                x = float(row[0])
                y = float(row[1])
                z = float(row[2])
                delam_outlist.append((x, y, z))
                pass
            except ValueError:
                pass
            pass
        pass

    if len(delam_outlist) == 0:
        raise ValueError("Could not parse any lines from CSV file %s" % (delam_outline))

    if delam_outlist[0] != delam_outlist[-1]:
        raise ValueError(
            "Delamination outline from %s does not form a closed wire (first and last vertices do not match)" % (
                delam_outline))

    # If there are n entries in the delam_outlist, one of which is doubled (start and end). There will be n-1 segments
    # We need to ignore the last point if we want to make a periodic spline
    delam_outpointsHArray = TColgp_HArray1OfPnt(1, len(delam_outlist)-1)

    for pos in range(len(delam_outlist)-1):
        current_point = gp_Pnt(delam_outlist[pos][0], delam_outlist[pos][1], delam_outlist[pos][2])
        delam_outpointsHArray.SetValue(pos + 1, current_point)
        pass

    if hasattr(delam_outpointsHArray,"GetHandle"): # Backward compatibility
        delm_outpointsHArray=delam_outpointsHArray.GetHandle()
        pass

    # Interpolate the points to make a closed curve
    interpAPI = GeomAPI_Interpolate(delam_outpointsHArray, True, tolerance)
    interpAPI.Perform()
    if interpAPI.IsDone():
        delam_curve = interpAPI.Curve()
    else:
        raise ValueError("Curve interpolation failed\n")

    if not hasattr(delam_curve,"FirstParameter"):
        # Backwards Compatibility with old pythonocc
        delam_curve = delam_curve.GetObject()
        pass
    
    # Analyze the delam_curve for minimum curvature
    numPts = 500
    uStart = delam_curve.FirstParameter()
    uEnd = delam_curve.LastParameter()
    minCurvature = 100
    for i in range(numPts):
        u = i/(numPts*(uEnd-uStart)) + uStart
        p = gp_Pnt(0,0,0)
        v1 = gp_Vec(0,0,0)
        v2 = gp_Vec(0,0,0)
        delam_curve.D2(u,p,v1,v2)
        currentCurvature = v2.Magnitude()

        if (currentCurvature < minCurvature):
            minCurvature = currentCurvature
            pass

        pass
    print("Curve parameter range : %f to %f"%(uStart, uEnd))
    print("Minimum curvature of the delamination outline is : %f"%(minCurvature))

    # Convert the curve to edge and then to Shape
    delam_out_edge = BRepBuilderAPI_MakeEdge(delam_curve).Edge()
    WireBuilder = BRepBuilderAPI_MakeWire()
    WireBuilder.Add(delam_out_edge)
    WireShape = WireBuilder.Shape()
    assert (WireShape.Closed())

    return WireShape


def ReplaceFacesWithImprintedSubfacesInLayerBodyFaceList(ImprintedFaces,LayerBodyFaceList,PointTolerance):
    """ImprintedFaces are faces from a "Fuse" operation, that are equivalent to or subfaces of
    elements in the LayerBodyFaceList. This function identifies all matching faces in the LayerBodyFaceList
    corresponding to all ImprintedFaces, and performs  replacement of the LayerBodyFaceList
    elements, with one or more elements from ImprintedFaces replacing each matched entry in LayerBodyFaceList.
    Returns the updated LayerBodyFaceList """


    AddToLayerBodyFaceList=[]
    RemoveFromLayerBodyFaceListIndices=set([])
    for faceShape in ImprintedFaces:
        layerBodyFace = LayerBodyFace.FromOCC(topods_Face(faceShape),"OFFSET")

        #print("Imprint: Face of %s: Got point %s" % (LayerBodyFaceList[0].Owner.Name,str(layerBodyFace.Point)))

        
        # Loop through all faces in this list
        for layer1BodyFaceIndex in range(len(LayerBodyFaceList)):
            layer1BodyFace=LayerBodyFaceList[layer1BodyFaceIndex]
            
            # Check if the reference point from fused face matches any other face
            pointClassification = OCCPointInFace(layerBodyFace.Point, layer1BodyFace.Face,PointTolerance)
            #print(pointClassification, TopAbs_IN, TopAbs_OUT, TopAbs_ON)
            #print("Face of %s: pointClassification=%d" % (LayerBodyFaceList[layer1BodyFaceIndex].Owner.Name,pointClassification))
            if (pointClassification == TopAbs_IN):
                
                #if (layerBodyFace.Point[0] >-5.94 and layerBodyFace.Point[0] <-5.93):
                #    print("Got troublesome point")
                #    if "Split2" in LayerBodyFaceList[layer1BodyFaceIndex].Owner.Name:
                #        print("writing debug.step")
                #        step_writer=STEPControl_Writer()
                #        step_writer.Transfer(layerBodyFace.Face,STEPControl_ShellBasedSurfaceModel,True)
                #        step_writer.Transfer(layer1BodyFace.Face,STEPControl_ShellBasedSurfaceModel,True)
                #        step_writer.Write('/tmp/debug.step')
                #
                #        pass
                
                #print("Found a matched face in %s "%layer1BodyFace.Owner.Name)
                layerBodyFace.Direction = layer1BodyFace.Direction
                layerBodyFace.Owner = layer1BodyFace.Owner
                layerBodyFace.BCType = layer1BodyFace.BCType
                AddToLayerBodyFaceList.append(layerBodyFace)
                RemoveFromLayerBodyFaceListIndices.add(layer1BodyFaceIndex)
                pass
            
            pass
        
        pass

    ReturnFaceList=copy.copy(LayerBodyFaceList)
    # Go through RemoveFromLayerBodyFaceListIndices from largest to smallest, removing them
    for RemoveIndex in sorted(RemoveFromLayerBodyFaceListIndices,reverse=True):
        del ReturnFaceList[RemoveIndex]
        pass

    # Add in all the faces to be added
    ReturnFaceList.extend(AddToLayerBodyFaceList)
    return ReturnFaceList

class OCCModelBuilder(object):
    """The OCCModelBuilder contains global parameters 
    for constructing laminates using the OpenCascade geometry kernel plus various methods.
    
    When bonding layers: 
     * call OCCMB.applydelaminations(layer1,layer2,delamlist) for all layer pairs in the bonding step
       to perform the layer imprint and face region assignment operations. 
     * call FAL=OCCMB.adjacentlayers(layer1,layer2,"COHESIVE") to read out the face adjacency list. 
     * Once the FAL is read out, neither layer can have its layerbodies split (fiber breakage) because
       that would require changing the layerbody names stored in the FAL. 
     * Thus any layer splitting must occur prior to the first adjacentlayers() call for that layer. 
     * Layer (body) splitting would also require an imprint operation on any adjoining layers. 
    """

    # Class members
    PointTolerance = None # Point positioning tolerance in mm
    NormalTolerance = None # Tolerance used for normal vectors
    Debug = None  # Flag for enabling troubleshooting output
    NextUnique = None # Next number to be generated by the Unique() method
    GapWidth = None # size of NOMODEL gap
    STLMeshSize = None
    
    # Constructor
    def __init__(self,**kwargs):
        # Default values:
        self.PointTolerance=1e-5
        self.NormalTolerance=1e-6
        self.Debug=False
        self.NextUnique=0
        self.GapWidth=0.5 # default gap of 0.5 mm
        self.STLMeshSize=3.0
        
        for argname in kwargs:
            if not hasattr(self,argname):
                raise ValueError("Invalid attribute")
            setattr(self,argname,kwargs[argname])
            pass
        pass

    def ProjectEdgesOntoFace(self,edge_edges,face):
        edge_curves = [ BRep_Tool.Curve(edge) for edge in edge_edges ]

        surface = BRep_Tool.Surface(face)
        
        
        # Note: edge_curves[i][1] and edge_curves[i][2] appear to be start and end u coordinates for curve

        Projections = [ GeomProjLib.geomprojlib_Project(edge_curve[0],surface) for edge_curve in edge_curves ]
        # If this fails with a RunTimeError, it probably means that your delamination outline does not project down onto the layer

        # Right here we should be trimming our projection to line up with layerbodyface1 and the unprojected edge (element of edge_edges)
        # But it's probably OK not to, because we are using the projection to make a tool that will be used to cut the face
        # and the extension of the tool beyond the face boundary shouldn't cause any problems, at least so long as thath
        # geometry doesn't get too weird
        
        ProjectionEdges = [ BRepBuilderAPI.BRepBuilderAPI_MakeEdge(Projection).Edge() for Projection in Projections ]

        # If we did trimmming, we would need to construct wire from the edge(s) that actually projected to something within the face,
        # with any gaps filled by appropriately trimmed edges from the face.

        #ProjectedWireBuilder = BRepBuilderAPI.BRepBuilderAPI_MakeWire()
        #
        #for ProjectionEdge in ProjectionEdges:
        #    ProjectedWireBuilder.add(ProjectionEdge)
        #    pass
        #
        # ProjectedWire = ProjectedWireBuilder.Wire()

        # Need to take ProjectionEdges, which are located on layerbodysurface1
        # and perform a offset from the surface
        return ProjectionEdges
    
    def OffsetFaceInBothDirections(self,face):
        OffsetDist = 100000.0*self.PointTolerance
        #OffsetDist = 1.0  # For debug
        # (do we need to convert face into a shell?) 
        mkOffset1 = BRepOffsetAPI.BRepOffsetAPI_MakeOffsetShape(face, OffsetDist, self.PointTolerance,
                                                                BRepOffset.BRepOffset_Skin,
                                                                False, False,
                                                                GeomAbs_Arc)
        assert (mkOffset1.IsDone())
        
        OffsetShell1 = mkOffset1.Shape()
        
        OffsetShell1FacesExp = TopExp_Explorer(OffsetShell1,TopAbs_FACE)
        OffsetFaces1 = []

        while OffsetShell1FacesExp.More():
            OffsetFaces1.append(topods_Face(OffsetShell1FacesExp.Current()))
            OffsetShell1FacesExp.Next()
            pass
        
        assert(len(OffsetFaces1)==1) # Offset of a single face should give a single face


        mkOffset2 = BRepOffsetAPI.BRepOffsetAPI_MakeOffsetShape(face, -OffsetDist, self.PointTolerance,
                                                                BRepOffset.BRepOffset_Skin,
                                                                False, False,
                                                                GeomAbs_Arc)
        assert (mkOffset2.IsDone())
        
        OffsetShell2 = mkOffset2.Shape()
        
        OffsetShell2FacesExp = TopExp_Explorer(OffsetShell2,TopAbs_FACE)
        OffsetFaces2 = []

        while OffsetShell2FacesExp.More():
            OffsetFaces2.append(topods_Face(OffsetShell2FacesExp.Current()))
            OffsetShell2FacesExp.Next()
            pass
        
        assert(len(OffsetFaces2)==1) # Offset of a single face should give a single face

        return (OffsetFaces1[0],OffsetFaces2[0])
    

    
    # Class methods
    def Unique(self):
        """Return a unique number (increasing sequence starting at 0)"""
        
        UniqueVal=self.NextUnique
        self.NextUnique += 1
        return UniqueVal
    

    def eval_face_pairs(self,layerbody1,layerbody2):
        """Evaluate the face pairs common between layerbody1 and layerbody2.
        Return a dictionary indexed by the common faces from layerbody1 that 
        identifies the corresponding faces from layerbody2"""

        # Right now we use equivalent surface geometry to determine equality of faces
        # this is NOT adequate.
        #
        # What we need is:
        #   * Equality comes from same surface, with same number of edges, where edges
        #     come from same wires, with same vertices.
        # Moreover:
        #   * If the geometry is equivalent, but the faces are not, then we may need
        #     to do an imprint.
        #   * One possibility is to attempt to imprint all faces of layerbody1 onto
        #     all faces of layerbody2 with the same underlying geometry, then
        #     vice-versa, then match up those common faces. 
        #   * Better option: See case #12 from 
        #     https://www.opencascade.com/doc/occt-7.0.0/overview/html/occt_user_guides__boolean_operations.html#occt_algorithms_9_4_12
        #     (fuse operations cause imprinting) so if we simply try to fuse the faces 
        #     with common geometry from both sides. The outcome should include all faces
        #     imprinted        

        # we created imprint_layers
        # that will imprint faces from lb1 onto lb2 and vice versa,
        # and then rebuilds the layerbodies of both lb1 and lb2 
        #
        # It current does this based on the assumption that the underlying geometries
        # of lb1 and lb2 come from the same source (same underlying surface object)
        # It might be better to do this by first comparing geometry, then
        # doing a rigorous equality comparison via edges and vertices (as described above).
        # In any case, it 
        # takes all faces from layerbody1 that share and underlying surface with layerbody2
        # and all of those faces from layerbody2, Then performing a fuse operation, yielding 
        # a large number of pieces. Then it reconstructs layerbodies 1 and 2 from the correct
        # subpieces. A correct subpiece is one that, given a non-boundary point on the subpiece,
        # the point lies inside a pre-existing face of that layerbody
        #
        # Then imprint_layers updates the face lists accordingly and sews the correct subpieces
        # back together into the layerbody.
        #

        # eval_face_pairs that assumes
        # and SHOULD check that faces are properly imprinted against each
        # other 
        # *!!!*** TODO: Check that all faces are properly imprinted


        
        FaceListTotal1 = layerbody1.FaceListOrig + layerbody1.FaceListOffset + layerbody1.FaceListSide
        FaceListTotal2 = layerbody2.FaceListOrig + layerbody2.FaceListOffset + layerbody2.FaceListSide
        
        CommonFaces = {}
        for face1 in FaceListTotal1:
            if face1 in FaceListTotal2:
                face2=FaceListTotal2[FaceListTotal2.index(face1)]

                # Note that face1 and face2 pass equality (== operator) because
                # equality operator identifies equivalent geometry. But the face1 and
                # face2 objects are actually different!
                
                CommonFaces[face1]=face2
                pass
            pass

        return CommonFaces

    def _imprint_delamination_split_nomodel(self,split_faces, BCTypes,NoModelToolShape,NoModelRefParamFace,parScale):

        assert(len(split_faces)==1)  # Assume we are starting with a single face
        
        
        # Use NoModelToolShape to do the split, operate on split_faces[0]
        NoModelSplitter=GEOMAlgo_Splitter()
        NoModelSplitter.AddArgument(topods_Face(split_faces[0]))
        NoModelSplitter.AddTool(NoModelToolShape)
        NoModelSplitter.Perform()
        
        # Empty out split_faces and BCTypes list. We will refill them from the pieces
        split_faces = []
        BCTypes=[]
        
        NoModelSplitCompound = NoModelSplitter.Shape()
        
        # Iterate over the pieces that have been split. Some of these will
        # be for CONTACT b.c.'s, and some NOMODEL.
        NoModel_split_exp=TopExp_Explorer(NoModelSplitCompound,TopAbs_FACE)
        
        # We tell the difference by using a previously created reference face in parametric coordinates
        # that includes solely the CONTACT zone. We determine a parametric coordinates point for
        # each of these faces, and check it against the previously determined reference face
        while NoModel_split_exp.More():
            nomodel_split_face_shape = topods_Face(NoModel_split_exp.Current())
            
            (NoModel_Split_Point,NoModel_Split_Normal,NoModel_Split_ParPoint) = FindOCCPointNormal(nomodel_split_face_shape,self.PointTolerance,self.NormalTolerance)
            split_faces.append(nomodel_split_face_shape)
            
            # Debugging only
            # nmsp_vertex = BRepBuilderAPI.BRepBuilderAPI_MakeVertex(gp_Pnt(NoModel_Split_Point[0],NoModel_Split_Point[1],NoModel_Split_Point[2])).Vertex()
            # end_point = np.array((NoModel_Split_ParPoint[0]*parScale, NoModel_Split_ParPoint[1]*parScale, 0)) + np.array((0,0,1))*0.1
            # nmsp_line = GC_MakeSegment(gp_Pnt(NoModel_Split_ParPoint[0]*parScale,NoModel_Split_ParPoint[1]*parScale,0.0),gp_Pnt(end_point[0],end_point[1],end_point[2])).Value()
            # nmsp_edge = BRepBuilderAPI_MakeEdge(nmsp_line)
            # nmsp_edge.Build()
            # print("Saving breps!")
            #
            # # Create new compound and make it valid.
            # compoundBuilder = BRep_Builder()
            # compoundShape = TopoDS_Compound()
            # compoundBuilder.MakeCompound(compoundShape)
            #
            # # The arguments order is: where to add, what to add.
            # compoundBuilder.Add(compoundShape, NoModelRefParamFace)
            # compoundBuilder.Add(compoundShape, RefParamFace)
            # compoundBuilder.Add(compoundShape, nmsp_edge.Shape())
            #
            # breptools_Write(compoundShape,"../data/lineface%d.brep" % (len(BCTypes)))
            # #breptools_Write(nmsp_vertex,"/data/point.brep")
            # #breptools_Write(NoModelRefParamFace,"../data/face%d.brep" % (len(BCTypes)))
            #
            # end debugging
            
            if (OCCPointInFace((NoModel_Split_ParPoint[0]*parScale,NoModel_Split_ParPoint[1]*parScale,0.0),NoModelRefParamFace,self.PointTolerance) == TopAbs_IN):
                # This particular nomodel_split_face is a contact zone.
                BCTypes.append("CONTACT")
                pass
            else:
                BCTypes.append("NONE")
                pass
            
            NoModel_split_exp.Next()
            pass
        return (split_faces,BCTypes)
    

    def _imprint_delaminations_assign_bcs(self,layerbodyface,SplitFace,ToolShapes,parScale):
        """ SplitFace is an OpenCascade face that has been split into pieces by imprinting. 
        Need to identify and apply the appropriate boundary conditions"""
        
        #step_writer2=STEPControl_Writer()
        #step_writer2.Transfer(SideShape,STEPControl_ShellBasedSurfaceModel,True)
        ##step_writer2.Transfer(layerbody.Shape, STEPControl_ManifoldSolidBrep, True)
        ##step_writer2.Transfer(layerbody2.Shape, STEPControl_ManifoldSolidBrep, True)
        #step_writer2.Transfer(SplitFace,STEPControl_ShellBasedSurfaceModel,True)
        #step_writer2.Write("../data/allShapes.STEP")

        split_face_exp=TopExp_Explorer(SplitFace,TopAbs_FACE)
        # Iterate over all faces
        numsplitfaces = 0
        split_face_shapes=[]
        split_layerbodyfaces=[]

        #step_writer = STEPControl_Writer()
        # TODO: sdh: refactor from here down, add NOMODEL switch. ***!!!

        # Go through the split faces
        while split_face_exp.More():
            split_face_shape=split_face_exp.Current()
            split_face_shapes.append(split_face_shape)

            
            #step_writer.Transfer(split_face_shape, STEPControl_ShellBasedSurfaceModel, True)

            # Since this face may be subdivided by a NOMODEL zone split,
            # need a place to keep track of the pieces, this
            # split_faces list.
            # Also keep a parallel list of BCTypes
            split_faces = [ topods_Face(split_face_shape) ]
            BCTypes= [ "TIE" ] # default
            
            (Point,Normal,ParPoint) = FindOCCPointNormal(split_faces[0],self.PointTolerance,self.NormalTolerance)

            # Match the split face with all of the ToolShapes to see if this particular split_face is part of
            # any of the tools (delaminations). If it is, then it is inside a delaminated region and needs to
            # be split into its CONTACT and NOMODEL zones, and those pieces must be identified and
            # assigned BCTypes of CONTACT or NOMODEL

            for (ToolShape, NoModelToolShape, RefParamFace, NoModelRefParamFace) in ToolShapes:
                
                # RefParamFace is in a 2D world of the (u,v) parameter space of the underlying surface,
                # mapped to the (x,y) plane. 
                if (OCCPointInFace((ParPoint[0]*parScale,ParPoint[1]*parScale,0.0),RefParamFace,self.PointTolerance) == TopAbs_IN):
                    # Matched! This particular split_face is a delamination zone.
                    # Need to do another split for the no model zone

                    if self.GapWidth == 0: # GapWidth of 0 disables NOMODEL zone
                        BCTypes[0]="CONTACT"
                        pass
                    else: 
                        (split_faces, BCTypes) = self._imprint_delamination_split_nomodel(split_faces, BCTypes,NoModelToolShape,NoModelRefParamFace,parScale)
                        pass
                    # Since we've confired that this subface is a delaminated region and created a NOMODEL zone (if applicable),
                    # we don't need to look further at whether this subface is a delaminated region

                    # ***!!! Possible bug: If this subface is part of two delaminations, then wouldn't those
                    # give different NoModelToolShapes? ... I think this will only work correctly for
                    # non-overlapping delaminations
                    break
                
                pass


            # go through list of subfaces and corresponding bctypes,
            # and create and append a LayerBodyFace to the split_layerbodyfaces list for each
            
            for facecnt in range(len(split_faces)):
                split_face=split_faces[facecnt]
                BCType=BCTypes[facecnt]

                (Point,Normal,ParPoint) = FindOCCPointNormal(split_face,self.PointTolerance,self.NormalTolerance)

                split_layerbodyfaces.append(LayerBodyFace(Face=split_face,
                                                          Point=Point,
                                                          Normal=Normal,
                                                          ParPoint=ParPoint,
                                                          Direction=layerbodyface.Direction,
                                                          Owner=layerbodyface.Owner,
                                                          BCType=BCType)) 
                numsplitfaces = numsplitfaces +1
                pass
            

            split_face_exp.Next()
            pass


        print("Number of split faces %d"%(numsplitfaces))

        if numsplitfaces==0:
            raise ValueError("No split faces")

        return split_layerbodyfaces
        # !!!***  Need to set BCTType on each generate LayerBodyFace !!!***

        # Create a reference face using the ProjectionEdges and layerbodyface.Face
        # Project the outline onto the face to figure out Boundary Conditions
        #projectionEdges = self.ProjectEdgesOntoFace(edge_edges, layerbodyface.Face)

        #Make wire from edges
        #WireBuilder = BRep_Builder()  # !!!*** Are build and Perimeter still necessary????
        #WirePerimeter = TopoDS_Compound()
        #WireBuilder.MakeCompound(WirePerimeter)

        #delamWire = TopoDS_Wire()
        #WireBuilder.MakeWire(delamWire)

        #for edgecnt in range(len(edge_edges)):
        #    projectionEdge = projectionEdges[edgecnt]
        #    WireBuilder.Add(delamWire, projectionEdge)
        #    pass

        #delamSurface = BRep_Tool.Surface(layerbodyface.Face)
        #FaceBuilder = BRepBuilderAPI_MakeFace(delamSurface, self.PointTolerance)
        #delamFace = FaceBuilder.Face()
        #FaceBuilder.Add(delamWire)

        #error = FaceBuilder.Error()
        #print("Error : %d"%(error))
        #if (error != BRepBuilderAPI_Error.BRepBuilderAPI_FaceDone):
        #    print("Face generation failed!")

        #step_writer2=STEPControl_Writer()
        #step_writer2.Transfer(delamFace,STEPControl_ShellBasedSurfaceModel,True)
        #step_writer2.Transfer(delamWire,STEPControl_GeometricCurveSet,True)
        #step_writer2.Write("../data/allShapes.STEP")

        # (Could also do similar process on side faces, but how could we ever get a delamination on the side faces???)
    
    def _imprint_delaminations_update_layerbody(self,layerbody,layerbodyface,split_layerbodyfaces):
        """ Update the layerbody, replacing the single face "layerbodyface" 
        with the collection of split_layerbodyfaces. layerbodyface might be on 
        either ORIG or OFFSET sides, so we try both possibilities"""
        
        if layerbodyface in layerbody.FaceListOrig:
            # Remove original face
            print("Removing original face from orig side")
            del layerbody.FaceListOrig[layerbody.FaceListOrig.index(layerbodyface)]
            # Add new faces
            for split_layerbodyface in split_layerbodyfaces:
                layerbody.FaceListOrig.append(split_layerbodyface)
                pass
            pass
        
        if layerbodyface in layerbody.FaceListOffset:
            # Remove offset face
            print("Removing original face from offset side")
            del layerbody.FaceListOffset[layerbody.FaceListOffset.index(layerbodyface)]
            # Add new faces
            for split_layerbodyface in split_layerbodyfaces:
                layerbody.FaceListOffset.append(split_layerbodyface)
                pass
            pass



        layerbody.Rebuild_Shape()
        pass
    
    def _imprint_delaminations_build_NOMODEL_tool(self,ToolShape,layerbodyface,parScale):
        """ Build the tool for cutting out the NOMODEL zone. ToolShape was the surface that we 
        intersected with the face to cut out the original delamination. The NOMODEL tool is made by 
        offsetting ToolShape by gapwidth. We try offsetting in both directions and select the one which
        gives us lower surface area, i.e. the inside offset. We then intersect the NOMODEL tool with 
        layerbodyface.face to get the intersection curve on the face, then obtain the intersection 
        curve in (u,v) coordinates. We use the (u,v) version to construct the nomodel reference param face
        (RefNoModelParamFace)... which is used in identifying CONTACT vs NOMODEL boundary condition.""" 

        # !!!*** (maybe we should just identify NOMODEL zone by moving epsilon from original cutting curve?)
        
        
        # Try offsetting both ways and see which gives the smaller surface area 
        mkOffset1 = BRepOffsetAPI.BRepOffsetAPI_MakeOffsetShape(ToolShape, -self.GapWidth, self.PointTolerance,
                                                                BRepOffset.BRepOffset_Skin,
                                                                True, True,
                                                                GeomAbs_Arc)
        assert (mkOffset1.IsDone())
        
        mkOffset2 = BRepOffsetAPI.BRepOffsetAPI_MakeOffsetShape(ToolShape, self.GapWidth, self.PointTolerance,
                                                                BRepOffset.BRepOffset_Skin,
                                                                True, True,
                                                                GeomAbs_Arc)
        
        assert (mkOffset2.IsDone())
        
        GProps1=GProp_GProps()
        brepgprop_SurfaceProperties(mkOffset1.Shape(),GProps1)
        SurfArea1=GProps1.Mass()
        
        GProps2=GProp_GProps()
        brepgprop_SurfaceProperties(mkOffset2.Shape(),GProps2)
        SurfArea2=GProps2.Mass()
        
        if SurfArea1 > SurfArea2:
            NoModelToolShells = mkOffset2.Shape()
            pass
        else:
            NoModelToolShells = mkOffset1.Shape()
            pass
        
        
        # The NoModelToolShape is of type Shell.
        # Need to extract the face from it and make sure there is only one face
        
        faceExplorer = TopExp_Explorer(NoModelToolShells, TopAbs_FACE)
        
        # Iterate over all edges
        # return a list of edges
        
        face_list = []
        while faceExplorer.More():
            current_face = topods_Face(faceExplorer.Current())
            face_list.append(current_face)
            faceExplorer.Next()
            pass
        
        if (len(face_list) != 1):
            raise ValueError("Face offset for NoModelToolShape created multiple faces!")
        
        NoModelToolShape = face_list[0]
        
        # step_writer2=STEPControl_Writer()
        # step_writer2.Transfer(layerbodyface.Face,STEPControl_ShellBasedSurfaceModel,True)
        # step_writer2.Transfer(NoModelToolShape,STEPControl_ShellBasedSurfaceModel,True)
        # #
        # step_writer2.Transfer(mkOffset2.Shape(),STEPControl_ShellBasedSurfaceModel,True)
        # #
        # step_writer2.Transfer(RefParamFace,STEPControl_ShellBasedSurfaceModel,True)
        # #
        # step_writer2.Transfer(RefNoModelParamFace,STEPControl_ShellBasedSurfaceModel,True)
        # step_writer2.Write("../data/OffsetTest.STEP")
        #
        # sys.modules["__main__"].__dict__.update(globals())
        # sys.modules["__main__"].__dict__.update(locals())
        # raise ValueError("Break")
        
        # Intersect the NoModelToolShape with the face to create the NoModelRefParamFace
        
        NoModelWireEdges = FaceFaceIntersect(NoModelToolShape, layerbodyface.Face)
        
        # Create reference parametric face for the no model zone using the NoModelWireShape
        RefNoModelParamFace = CreateReferenceFace(NoModelWireEdges, layerbodyface.Face, parScale, self.PointTolerance)
        
        return (NoModelToolShape,RefNoModelParamFace)

    
    def imprint_delaminations(self,layerbody,layerbodyface,delam_outlines):
        """Given a first layerbody and corresponding face, and a list of delamination outlines 
        (loop of 3D coordinates, hopefully projected onto the faces): 
         A. Loop over each delam outline
           1. Identify the regions of layerbodyface that are inside delam_outlines
           2. Create an offset curve a distance of self.GapWidth inside delam_outlines
           3. Imprint both delam_outlinse and the offset curves
           4. Identify the regions of layerbodyface that are between delam_outlines and the offset
              curves to have "NOMODEL" BCType unless they already had "CONTACT" BCType 
           5. Identify the regions of facebody1/2 that are inside the offset curves to have "CONTACT"
              BCType
          B. Update layerbody with newly constructed
              imprinted faces, marked as determined above. 
           1. The updated layerbody may have layerbodyface 
              replaced/subdivided, but other faces in layerbody remain unchanged. 
"""
        
        # NOTE: May need additional parameters (adjacent surfaces or faces?) to do the
        # delam_outline offset curve?

        # NOTE: When regenerating layerbodies, do NOT give them new names unless they are being
        # split (which they aren't from this function)

        parScale = 100

        if len(delam_outlines)==0:
            return # nothing to do if there are no delaminations

        # Offset Face in both directions to create new faces for projection
        (bounding_face_a,bounding_face_b) = self.OffsetFaceInBothDirections(layerbodyface.Face)
        
        ToolShapes=[]
        
        for delam_outline in delam_outlines:
            # ***!!! Still need to implement in-plane offsets

            # Loading WireShape directly from a STEP file
            #WireShape = loaders.load_byfilename(os.path.join("..","data","Delam1.STEP"))

            # delam_outline is the file name
            WireShape = CreateDelaminationWire(delam_outline, self.PointTolerance)

            # Offset wire to create no model shape

            # Offset the tool shape to create the inner tool shape for the no model zone
            #NoModelDist = -1
            #mkOffset = BRepOffsetAPI.BRepOffsetAPI_MakeOffset(topods_Wire(WireShape), GeomAbs_Arc, False)
            #mkOffset.Perform(NoModelDist)

            #assert (mkOffset.IsDone())
            #WireInShape = mkOffset.Shape()

            #step_writer2=STEPControl_Writer()
            #step_writer2.Transfer(WireShape,STEPControl_GeometricCurveSet,True)
            #step_writer2.Transfer(WireInShape,STEPControl_GeometricCurveSet,True)
            #step_writer2.Write("../data/allShapes.STEP")

            #sys.modules["__main__"].__dict__.update(globals())
            #sys.modules["__main__"].__dict__.update(locals())
            #raise ValueError("Break")


            exp=TopExp_Explorer(WireShape,TopAbs_EDGE)
            
            # Iterate over all edges
            edge_shapes=[]
            while exp.More():
               edge_shapes.append(exp.Current())

               exp.Next()
               pass
            edge_edges = [ topods_Edge(edge_shape) for edge_shape in edge_shapes ]

            
            ProjectionEdges_a = self.ProjectEdgesOntoFace(edge_edges,bounding_face_a)
            ProjectionEdges_b = self.ProjectEdgesOntoFace(edge_edges,bounding_face_b)

            # Project the curve on to the face to figure out the parametric curve
            ProjectionEdges = self.ProjectEdgesOntoFace(edge_edges, layerbodyface.Face)

            # Use the parametric curve to create a reference face that is inside the curve, i.e.
            # inside the delaminated zone, so that we can compare it with the faces that will be generated
            # by the splitting tool 
            RefParamFace = CreateReferenceFace(ProjectionEdges, layerbodyface.Face, parScale, self.PointTolerance)

            # Generate faces connecting original and projected edges.
            # We will use this as a tool to do the cut. 
            # For the moment assume only one edge
        
            build=BRep_Builder()
            Perimeter=TopoDS_Compound()
            build.MakeCompound(Perimeter)
            
            wire_a = TopoDS_Wire()
            build.MakeWire(wire_a)
            wire_b = TopoDS_Wire()
            build.MakeWire(wire_b)
        
            for edgecnt in range(len(edge_edges)):
                projectionedge_a = ProjectionEdges_a[edgecnt]
                projectionedge_b = ProjectionEdges_b[edgecnt]
                
                
                build.Add(wire_a,projectionedge_a)
                build.Add(wire_b,projectionedge_b)
                pass
            
            # Generate side faces
            ToolGenerator = BRepOffsetAPI.BRepOffsetAPI_ThruSections()
            ToolGenerator.AddWire(wire_a)
            ToolGenerator.AddWire(wire_b)
            ToolGenerator.Build()
            
            if (not ToolGenerator.IsDone()):
                raise ValueError("Tool side face generation failed\n")
        
            ToolShape = ToolGenerator.Shape()

            build.Add(Perimeter,ToolShape)

            # Offset the tool shape to create the inner tool shape for the no model zone

            if self.GapWidth > 0: # only create NOMODEL tool if we have a nonzero GapWidth

                (NoModelToolShape,RefNoModelParamFace) = self._imprint_delaminations_build_NOMODEL_tool(ToolShape,layerbodyface,parScale)

                pass
            else:
                # self.GapWidth == 0 -> We don't need the NOMODEL zone tools 
                NoModelToolShape = None
                RefNoModelParamFace = None
                pass
            
            # Create a Tuple to store the ToolShape and the NoModelToolShape, and RefParamFace -- used to identify the inside region
            ToolShapes.append((ToolShape, NoModelToolShape, RefParamFace, RefNoModelParamFace))
                
            pass
            
        # Split the face using the delamination zones
        GASplitter=GEOMAlgo_Splitter()
        GASplitter.AddArgument(topods_Face(layerbodyface.Face))
        for ToolShape in ToolShapes:
            GASplitter.AddTool(ToolShape[0])
            pass

        # This commented-out code is used to time and record the parameters
        # to the calls to the splitter, which seem to take a fraction
        # of a second for example 1 where the underlying layer is a
        # single piece, but about a minute each in example 2 where the
        # delamination crosses the body (fiber breakage) boundary (!)
        
        #global splittercall
        #if not "splittercall" in globals():
        #    splittercall=0
        #    pass
        #print("splitter call #%d" % (splittercall))
        #
        #breptools_Write(layerbodyface.Face,'/tmp/splitcall%d_face.brep' % (splittercall))
        #for ToolShapeNum in range(len(ToolShapes)):
        #    breptools_Write(ToolShapes[ToolShapeNum][0],'/tmp/splitcall%d_tool%d.brep' % (splittercall,ToolShapeNum))
        #    pass
        #starttime=time.clock()
        
        GASplitter.Perform()
        #elapsed=time.clock()-starttime
        #print("splitter call #%d took %f seconds" % (splittercall,elapsed))
        
        #splittercall += 1
        
        #if (not GASplitter.IsDone()):
        #    raise ValueError("Splitting face failed\n")

        SplitFace = GASplitter.Shape()
        # Hopefully this did not damage layerbodyface


        # Now that the face has been imprinted by the various delaminations
        # need to assign the boundary conditions and imprint nomodel zones, etc. 

        split_layerbodyfaces = self._imprint_delaminations_assign_bcs(layerbodyface,SplitFace,ToolShapes,parScale)
        
        
        # split_layerbodyfaces should have two or more faces (TopoDS_Shape of type Face

        # Modify layerbody in-place with the original layerbodyfaces except for this one,

        self._imprint_delaminations_update_layerbody(layerbody,layerbodyface,split_layerbodyfaces)
        

        #step_writer = STEPControl_Writer()
        #step_writer.Transfer(delamSolidShape, STEPControl_ManifoldSolidBrep, True)
        #step_writer.Write("../Data/Layers.step")

        # Also need to repeat the process for the other face

        # (so actual content of this function should be abstracted into
        # a new function we can call twice). 

        return 

    def apply_delaminations(self,layer1,layer2,delaminationlist):
        """Iterate over the layerbodies of layer1 and layer2. For each pair of layerbodies, 
        evaluate the common faces with eval_face_pairs(). For each common face pair and 
        each delamination, update the corresponding layerbodies with process_delamination. 
        
        When done, rerun eval_face_pairs and use the resulting dictionary to construct and 
        return a face adjacency list (FAL)"""

        self.imprint_layers(layer1,layer2) # make sure that layers are imprinted against each other

        for lb1cnt in range(len(layer1.BodyList)):
            lb1=layer1.BodyList[lb1cnt]

            for lb2cnt in range(len(layer2.BodyList)):
                lb2=layer2.BodyList[lb2cnt]

                # !!!*** modify eval_face_pairs into match_face_pairs
                # that will imprint faces from lb1 onto lb2 and vice versa,
                # and then return a dictionary indexed by the face in lb1
                # of the corresponding imprinted face in lb2
                CommonFaces=self.eval_face_pairs(lb1,lb2)

                for CommonFace in CommonFaces:
                    #print ("Imprinting delaminations onto CommonFace %s" % (str(CommonFace)))
                    self.imprint_delaminations(lb1,CommonFace,delaminationlist)
                    #print ("Imprinting delaminations onto CommonFaces[CommonFace] %s" % (str(CommonFaces[CommonFace])))
                    self.imprint_delaminations(lb2,CommonFaces[CommonFace],delaminationlist)
                    #raise ValueError("test")
                    pass

                pass
            pass
        pass

    def imprint_layers(self,layer1,layer2):
        """ Imprint the layerbodies of layer1 and layer2 onto each other
        First find all faces that share a common surface
        Perform Boolean Fuse operation on the faces
        Sort the faces based on underlying layerbody
        Reconstruct the layerbody and replace them in the layer
        """
        # Reconstruct layerbodies 1 and 2 from the correct
        # subpieces. A correct subpiece is one that, given a non-boundary point on the subpiece,
        # the point lies inside a pre-existing face of that layerbody
        #
        # Then the correct subpieces can be sewn back together into the layerbody,
        # and the face lists updates accordingly.

        # ASSUMPTION: layer1 offset faces are fused to layer 2 orig faces

        #layer1FaceList = TopTools_ListOfShape()
        #layer2FaceList = TopTools_ListOfShape()

        layer1FaceList=[]
        layer2FaceList=[]
        for LB in layer1.BodyList:
            for LayerBodyFaceList in [ LB.FaceListOrig, LB.FaceListOffset, LB.FaceListSide]: # Do we really need side?
                layer1FaceList.extend(LayerBodyFaceList)
                #for LayerBodyFace in LayerBodyFaceList:
                #    layer1FaceList.Append(LayerBodyFace.Face)
                #    pass
                pass
            pass

        for LB in layer2.BodyList:
            for LayerBodyFaceList in [ LB.FaceListOrig, LB.FaceListOffset, LB.FaceListSide]: # Do we really need side? 
                layer2FaceList.extend(LayerBodyFaceList)
                #for LayerBodyFace in LayerBodyFaceList:
                #    layer2FaceList.Append(LayerBodyFace.Face)
                #    pass
                pass            
            pass

        # Find faces with matched underlying surfaces...
        # (This can potentially be eliminated but it will
        # make the operation very slow)

        FusedFaces=[]
        
        SurfaceDict1={}
        for face in layer1FaceList:
            # OpenCascade surface provides a DumpToString/DumpJsonToString method
            # that gives a (presumably unique) identification of the
            # underlying surface. The string is hashable, so usable
            # as a dictionary key
            if hasattr(Geom_Surface,"DumptoString"):
                key = BRep_Tool().Surface(face.Face).DumpToString()
                pass
            else:
                key = BRep_Tool().Surface(face.Face).DumpJsonToString()
                pass
            if key not in SurfaceDict1:
                SurfaceDict1[key]=[]
                pass
            
            SurfaceDict1[key].append(face.Face)
            pass
        
        SurfaceDict2={}
        for face in layer2FaceList:
            # OpenCascade surface provides a DumpToString/DumpJsonToString method
            # that gives a (presumably unique) identification of the
            # underlying surface. The string is hashable, so usable
            # as a dictionary key
            if hasattr(Geom_Surface,"DumptoString"):
                key = BRep_Tool().Surface(face.Face).DumpToString()
                pass
            else:
                key = BRep_Tool().Surface(face.Face).DumpJsonToString()
                pass
            
            if key not in SurfaceDict2:
                SurfaceDict2[key]=[]
                pass
            
            SurfaceDict2[key].append(face.Face)
            pass

        for surfacekey in SurfaceDict1:
            if surfacekey in SurfaceDict2:
                # Found surface with faces both in layer 1 and 2
                layer1surfacefaces = TopTools_ListOfShape()
                for face in SurfaceDict1[surfacekey]:
                    layer1surfacefaces.Append(face)
                    pass
                
                layer2surfacefaces = TopTools_ListOfShape()
                for face in SurfaceDict2[surfacekey]:
                    layer2surfacefaces.Append(face)
                    pass
 
                
                # Fuse operation on all faces.
                fuser = BRepAlgoAPI_Fuse()
                fuser.SetArguments(layer1surfacefaces)
                fuser.SetTools(layer2surfacefaces)
                fuser.Build()
                fusedShape = fuser.Shape()
                # Separate out the faces from the fusedShape and add
                # to FusedFaces

                topExplorer = TopExp_Explorer(fusedShape, TopAbs_FACE)

                # iterate over pieces of the original wire
                while topExplorer.More():
                    FusedFaces.append(topExplorer.Current())
                    topExplorer.Next()
                    pass

                pass
            pass

        # !!!*** Need to loop through fused faces, and use them to
        # replace faces in layerbodies of layer1 and layer2
        # based on finding identification for the new faces and
        # testing whether that identification works on old faces.
        # ... Then replace the old faces with identified new faces.

        for LB in layer1.BodyList:
            LB.FaceListOrig=ReplaceFacesWithImprintedSubfacesInLayerBodyFaceList(FusedFaces,LB.FaceListOrig,self.PointTolerance)
            LB.FaceListOffset=ReplaceFacesWithImprintedSubfacesInLayerBodyFaceList(FusedFaces,LB.FaceListOffset,self.PointTolerance)
            LB.FaceListSide=ReplaceFacesWithImprintedSubfacesInLayerBodyFaceList(FusedFaces,LB.FaceListSide,self.PointTolerance)

            
            # Rebuild layerbody from faces
            LB.Rebuild_Shape()
            

            pass

        for LB in layer2.BodyList:
            LB.FaceListOrig=ReplaceFacesWithImprintedSubfacesInLayerBodyFaceList(FusedFaces,LB.FaceListOrig,self.PointTolerance)
            LB.FaceListOffset=ReplaceFacesWithImprintedSubfacesInLayerBodyFaceList(FusedFaces,LB.FaceListOffset,self.PointTolerance)
            LB.FaceListSide=ReplaceFacesWithImprintedSubfacesInLayerBodyFaceList(FusedFaces,LB.FaceListSide,self.PointTolerance)

            
            # Rebuild layerbody from faces
            LB.Rebuild_Shape()
            

            pass


        
        pass





        #step_writer=STEPControl_Writer()
        #step_writer.Transfer(layer1OffsetFaceList[0].Face,STEPControl_ShellBasedSurfaceModel,True)
        #step_writer.Transfer(layer2OrigFaceList[0].Face,STEPControl_ShellBasedSurfaceModel,True)
        #step_writer.Transfer(fusedShape,STEPControl_ShellBasedSurfaceModel,True)
        #for shape in FusedFaces:
        #    step_writer.Transfer(shape,STEPControl_ShellBasedSurfaceModel,True)
        #    pass
        
        #step_writer.Write("../data/fusedFace.STEP")
        
        #sys.modules["__main__"].__dict__.update(globals())
        #sys.modules["__main__"].__dict__.update(locals())
        #raise ValueError("Break")





        pass

    def adjacent_layer_boundary_conditions(self,layer1,layer2,bc_map=None):
        """ Once adjacent_layer_boundary_conditions() is called, the LayerBody's in EITHER layer can't be split
        anymore -- because then they might need new names,
        and the return values contain the layer body names that will be used
        to apply the boundary conditions. Returns the face adjacency list"""

        FAL = [] # Face Adjacency List

        self.imprint_layers(layer1,layer2)  # imprint_layers operates on the full layers because the bodies might be
        # subdivided differently on both sides
        for lb1 in layer1.BodyList:
            for lb2 in layer2.BodyList:
                #print("lb1 = %s; lb2=%s" % (lb1.Name,lb2.Name))
                FacePairs=self.eval_face_pairs(lb1,lb2)
                for postimprint_CommonFace in FacePairs:
                    BCType=postimprint_CommonFace.BCType

                    assert(lb1.Name is not None)
                    assert(lb2.Name is not None)
                    
                    if BCType is None: # BCType not otherwise set... insert default
                        BCType="TIE"
                        pass

                    if bc_map is not None:
                        # apply user-supplied BC mapping
                        BCType=bc_map[BCType]
                        pass
                    
                    print ("Appending to FAL: %s" % ({ "name1": lb1.Name,
                                                       "name2": lb2.Name,
                                                       "bcType": BCType,
                                                       "point1": postimprint_CommonFace.Point,
                                                       "normal1": postimprint_CommonFace.Normal,}))
                    
                    FAL.append( { "name1": lb1.Name,
                                  "name2": lb2.Name,
                                  "bcType": BCType,
                                  "point1": postimprint_CommonFace.Point,
                                  "normal1": postimprint_CommonFace.Normal,
                                  #"point2": FacePairs[postimprint_CommonFace].Point,
                                  #"normal2": FacePairs[postimprint_CommonFace].Normal,
                    })
                    pass
                
                pass
            
            pass
        
        return FAL

    def save_layer_surface_stl(self,stlfilename,layer1,layer2):
        """ Save an .STL of the boundary between layer1 and layer2...
        Must be called after finalization of BOTH layer1 and layer2."""
        self.imprint_layers(layer1,layer2)  # imprint_layers operates on the full layers because the bodies might be
        # subdivided differently on both sides

        shellbuilder = BRep_Builder();
        shell=TopoDS_Shell();
        shellbuilder.MakeShell(shell);

        # Should be able to connect these faces
        # just with a BRep_Builder because they should
        # already be topologically joined (otherwise
        # we may need to sew them together
        # with BRepBuilderAPI_Sewing())

        for lb1 in layer1.BodyList:
            for lb2 in layer2.BodyList:
                #print("lb1 = %s; lb2=%s" % (lb1.Name,lb2.Name))
                FacePairs=self.eval_face_pairs(lb1,lb2)

                               
                for postimprint_CommonFace in FacePairs:
                    # postimprint_CommonFace is the LayerBodyFace from lb1... use it to build a shell 
                    shellbuilder.Add(shell, postimprint_CommonFace.Face);
                    pass
                pass
            pass

        MeshDMObject = layer.Layer.CreateDMObject(shell,self.STLMeshSize)
        MeshDMObject.SaveSTL(stlfilename)
        pass
    
        
        
    def save(self,cad_file_name,to_be_saved):
        step_writer=STEPControl_Writer()
        BodyNameList=[]
        BodyNameSet=set([])

        #for layerobj in to_be_saved:
        #    for layerbodyobj in layerobj.BodyList:
                
        for layerbody_or_solid_obj in to_be_saved:
            step_writer.Transfer(layerbody_or_solid_obj.Shape,STEPControl_ManifoldSolidBrep,True)
            
            assert(not(layerbody_or_solid_obj.Name in BodyNameSet)) # verify nome layerbodyname reuse!
            
            BodyNameList.append(layerbody_or_solid_obj.Name)
            
            BodyNameSet.add(layerbody_or_solid_obj.Name)
            
            pass
        
        step_writer.Write(cad_file_name)
        
        
        return BodyNameList
    
    pass


if __name__=="__main__":
    MB=OCCModelBuilder(PointTolerance=1e-5,NormalTolerance=1e-6)
    
    pointTolerance = 1e-6
    Mold = layer.LayerMold.FromFile(os.path.join("..","data","CurvedMold1.STEP"))
    #Mold = layer.LayerMold.FromFile(os.path.join("..","data","FlatMold3.STEP"))
    Layer1=layer.Layer.CreateFromMold("Layer 1",Mold,2.0,"OFFSET",pointTolerance)
    #Layer1=layer.Layer.CreateFromMold("Layer 1",Mold,2.0,"ORIG",1e-6)
    Layer2=layer.Layer.CreateFromMold("Layer 2",Layer1.OffsetMold(),2.0,"OFFSET",pointTolerance)
    #Layer2=layer.Layer.CreateFromMold("Layer 2",Mold,2.0,"OFFSET",1e-6)
    Layer3=layer.Layer.CreateFromMold("Layer 3",Layer2.OffsetMold(),2.0,"OFFSET",pointTolerance)

    #Layer2.Split(os.path.join("..","data","SplitLine2.csv"), pointTolerance)


    #delaminationlist = [ ]
    delaminationlist = [ os.path.join("..","data","nasa-delam12-1.csv"), os.path.join("..","data","nasa-delam12-2.csv") ]

    #MB.imprint_layers(Layer1, Layer2)
    #MB.imprint_layers(Layer2, Layer3)
    #defaultBCType = 2
    #FAL = MB.adjacent_layer_boundary_conditions(Layer1,Layer2,defaultBCType)
    #MB.apply_delaminations(Layer1,Layer2,delaminationlist)


    step_writer=STEPControl_Writer()
    
    step_writer.Transfer(Layer1.BodyList[0].Shape,STEPControl_ManifoldSolidBrep,True)
    step_writer.Transfer(Layer2.BodyList[0].Shape,STEPControl_ManifoldSolidBrep,True)
    #step_writer.Transfer(Layer2.BodyList[1].Shape,STEPControl_ManifoldSolidBrep,True)
    step_writer.Transfer(Layer3.BodyList[0].Shape,STEPControl_ManifoldSolidBrep,True)
    step_writer.Write(os.path.join("..","data","Layers.step"))

    pass
