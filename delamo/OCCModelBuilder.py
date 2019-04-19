import sys
import copy
import os.path
import csv
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
from OCC.BRep import BRep_Builder
from OCC.BRep import BRep_Tool
from OCC.BRepExtrema import BRepExtrema_DistShapeShape
from OCC import BRepLib
from OCC import BRepOffsetAPI
from OCC import BRepOffset
from OCC import BRepBuilderAPI
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeEdge
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeWire
#from OCC.BRepClass import BRepClass_FacePassiveClassifier
from OCC.BRepClass import BRepClass_FaceExplorer
from OCC.BRepClass import BRepClass_FClassifier
from OCC.ShapeAnalysis import ShapeAnalysis_FreeBoundsProperties
from OCC.BRepTools import breptools_Read
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
from OCC.Geom import Geom_Line
from OCC.TopTools import TopTools_ListIteratorOfListOfShape
from OCC.TopTools import TopTools_ListOfShape
from OCC.GeomLProp import GeomLProp_SLProps
from OCC.gp import gp_Pnt2d
from OCC.gp import gp_Vec
from OCC.gp import gp_Dir
from OCC.gp import gp_Pnt
from OCC.GEOMAlgo import GEOMAlgo_Splitter
from OCC.BRepAlgoAPI import BRepAlgoAPI_Fuse
from OCC import GeomProjLib
from OCC.TColgp import TColgp_Array1OfPnt
from OCC.TColgp import TColgp_HArray1OfPnt
from OCC.GeomAPI import (GeomAPI_Interpolate, GeomAPI_PointsToBSpline)

from OCC.STEPControl import STEPControl_Reader
from OCC.STEPControl import STEPControl_Writer
from OCC.STEPControl import STEPControl_ManifoldSolidBrep
from OCC.STEPControl import STEPControl_Writer,STEPControl_ShellBasedSurfaceModel,STEPControl_GeometricCurveSet

from layer import LayerBody,LayerBodyFace
from layer import OCCPointInFace

import loaders
import layer


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
        # Loop through all faces in this list
        for layer1BodyFaceIndex in range(len(LayerBodyFaceList)):
            layer1BodyFace=LayerBodyFaceList[layer1BodyFaceIndex]
            
            # Check if the reference point from fused face matches any other face
            pointClassification = OCCPointInFace(layerBodyFace.Point, layer1BodyFace.Face,PointTolerance)
            #print(pointClassification, TopAbs_IN, TopAbs_OUT, TopAbs_ON)
            if (pointClassification == TopAbs_IN):
                #print("Found a matched face in %s "%layer1BodyFace.Owner.Name)
                layerBodyFace.Direction = layer1BodyFace.Direction
                layerBodyFace.Owner = layer1BodyFace.Owner
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
    
    # Constructor
    def __init__(self,**kwargs):
        # Default values:
        self.PointTolerance=1e-5
        self.NormalTolerance=1e-6
        self.Debug=False
        self.NextUnique=0
        
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

    def imprint_delaminations(self,layerbody,layerbodyface,delam_outlines):
        """Given a first layerbody and corresponding face, and a second 
        layerbody and corresponding face, and a list of delamination outlines 
        (loop of 3D coordinates, hopefully projected onto the faces): 
         A. Loop over each delam outline
           1. Identify the regions of facebody1/2 that are inside delam_outline
           2. Create an offset curve a distance of self.GapWidth inside delam_outline
           3. Imprint both delam_outline and the offset curve
           4. Identify the regions of facebody1/2 that are between delam_outline and the offset
              curve to have "NOMODEL" BCType unless they already had "CONTACT" BCType 
           5. Identify the regions of facebody1/2 that are inside the offset curve to have "CONTACT"
              BCType
          B. Generate new replacements for layerbody1 and layerbody2 with newly constructed
              imprinted faces, marked as determined above. 
           1. The replaced layerbody1 and layerbody2 may have facebody1 and facebody2 
              replaced/subdivided, but other faces in layerbodies1/2 should remain unchanged. 
          C. Return (replacement_layerbody1, replacement_layerbody2)"""

        # NOTE: May need additional parameters (adjacent surfaces or faces?) to do the
        # delam_outline offset curve?

        # NOTE: When regenerating layerbodies, do NOT give them new names unless they are being
        # split (which they aren't from this function)


        SideShapes=[]
        
        for delam_outline in delam_outlines:
            # ***!!! Temporarily load curve rather than constructing wire from delamination outline

            # ***!!! Still need to implement in-plane offsets
            # delam_outline is the file name
            delam_outlist = []
            with open(delam_outline) as csvfile:
                reader=csv.reader(csvfile,delimiter=',',quotechar='"')
                for row in reader:
                    if len(row) != 3:
                        raise ValueError("Malformed row in CSV file %s: %s" % (delam_outline,",".join(row)))
                    try:
                        x=float(row[0])
                        y=float(row[1])
                        z=float(row[2])
                        delam_outlist.append((x,y,z))
                        pass
                    except ValueError:
                        pass
                    pass
                pass

            if len(delam_outlist) == 0:
                raise ValueError("Could not parse any lines from CSV file %s" % (delam_outline))

            if delam_outlist[0] != delam_outlist[-1]:
                raise ValueError("Delamination outline from %s does not form a closed wire (first and last vertices do not match)" % (delam_outline))
            
            # If there are n entries in the delam_outlist, one of which is doubled (start and end). There will be n-1 segments
            delam_outpointsHArray = TColgp_HArray1OfPnt(1, len(delam_outlist))

            for pos in range(len(delam_outlist)):
                current_point = gp_Pnt(delam_outlist[pos][0],delam_outlist[pos][1],delam_outlist[pos][2])
                delam_outpointsHArray.SetValue(pos+1, current_point)
                pass

            # Interpolate the points to make a closed curve
            interpAPI = GeomAPI_Interpolate(delam_outpointsHArray.GetHandle(), False, self.PointTolerance)
            interpAPI.Perform()
            if interpAPI.IsDone():
                 delam_out_curve = interpAPI.Curve()
            else:
                raise ValueError("Curve interpolation failed\n")

            # Convert a curve to edge and then to Shape
            delam_out_edge = BRepBuilderAPI_MakeEdge(delam_out_curve).Edge()
            WireBuilder = BRepBuilderAPI_MakeWire()
            WireBuilder.Add(delam_out_edge)
            WireShape = WireBuilder.Shape()
            assert(WireShape.Closed())

            #step_writer2=STEPControl_Writer()
            #step_writer2.Transfer(WireShape,STEPControl_GeometricCurveSet,True)
            #step_writer2.Write("../data/Wire.STEP")
            #
            #sys.modules["__main__"].__dict__.update(globals())
            #sys.modules["__main__"].__dict__.update(locals())
            #raise ValueError("Break")

            # Loading WireShape directly from a STEP file
            #WireShape = loaders.load_byfilename(os.path.join("..","data","Delam1.STEP"))

            exp=TopExp_Explorer(WireShape,TopAbs_EDGE)
            
            # Iterate over all edges
            edge_shapes=[]
            while exp.More():
               edge_shapes.append(exp.Current())

               exp.Next()
               pass
            edge_edges = [ topods_Edge(edge_shape) for edge_shape in edge_shapes ]

            # Offset Face in both directions to create new faces for projection
            (bounding_face_a,bounding_face_b) = self.OffsetFaceInBothDirections(layerbodyface.Face)
            
            ProjectionEdges_a = self.ProjectEdgesOntoFace(edge_edges,bounding_face_a)
            ProjectionEdges_b = self.ProjectEdgesOntoFace(edge_edges,bounding_face_b)

            # Generate faces connecting original and projected edges.
            # We will use this as a tool to do the cut. 
            
            # For the moment assume only one edge
        
            build=BRep_Builder()  # !!!*** Are build and Perimeter still necessary????
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
            SideGenerator = BRepOffsetAPI.BRepOffsetAPI_ThruSections()
            SideGenerator.AddWire(wire_a)
            SideGenerator.AddWire(wire_b)
            SideGenerator.Build()
            
            if (not SideGenerator.IsDone()):
                raise ValueError("Side face generation failed\n")
        
            SideShape = SideGenerator.Shape()
            
            build.Add(Perimeter,SideShape)
            SideShapes.append(SideShape)
            pass
        
            
        GASplitter=GEOMAlgo_Splitter()
        GASplitter.AddArgument(topods_Face(layerbodyface.Face))
        for SideShape in SideShapes:
            GASplitter.AddTool(SideShape)
            pass
        
        GASplitter.Perform()

        #if (not GASplitter.IsDone()):
        #    raise ValueError("Splitting face failed\n")

        SplitFace= GASplitter.Shape()
        # Hopefully this did not damage layerbodyface
        
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
        step_writer = STEPControl_Writer()

        while split_face_exp.More():
            split_face_shape=split_face_exp.Current()
            split_face_shapes.append(split_face_shape)
            step_writer.Transfer(split_face_shape, STEPControl_ShellBasedSurfaceModel, True)

            split_face = topods_Face(split_face_shape)
            (Point,Normal,ParPoint) = layer.FindOCCPointNormal(split_face,self.PointTolerance,self.NormalTolerance)
            
            split_layerbodyfaces.append(LayerBodyFace(Face=split_face,
                                                      Point=Point,
                                                      Normal=Normal,
                                                      ParPoint=ParPoint,
                                                      Direction=layerbodyface.Direction,
                                                      Owner=layerbodyface.Owner,
                                                      BCType="TIE")) # !!!*** BCType needs to be set correctly ***!!!
            
            numsplitfaces = numsplitfaces +1

            split_face_exp.Next()
            pass
        step_writer.Write("/tmp/split_faces.step")

        print("Number of split faces %d"%(numsplitfaces))

        # split_face_shapes now should have two or more faces (TopoDS_Shape of type Face
        # Need to create a new layerbody with the original layerbodyfaces except for this one, and two new layerbodyfaces


        # Modify layerbody in-place

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
        
        

        # !!!***  Need to set BCTType on each generate LayerBodyFace !!!***
        
        # (Could also do similar process on side faces, but how could we ever get a delamination on the side faces???)

        #sys.modules["__main__"].__dict__.update(globals())
        #sys.modules["__main__"].__dict__.update(locals())
        #raise ValueError("Break")

        layerbody.Rebuild_Shape()
        

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
            # OpenCascade surface provides a DumpToString method
            # that gives a (presumably unique) identification of the
            # underlying surface. The string is hashable, so usable
            # as a dictionary key
            key = BRep_Tool().Surface(face.Face).DumpToString()

            if key not in SurfaceDict1:
                SurfaceDict1[key]=[]
                pass
            
            SurfaceDict1[key].append(face.Face)
            pass
        
        SurfaceDict2={}
        for face in layer2FaceList:
            # OpenCascade surface provides a DumpToString method
            # that gives a (presumably unique) identification of the
            # underlying surface. The string is hashable, so usable
            # as a dictionary key
            key = BRep_Tool().Surface(face.Face).DumpToString()

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
                FacePairs=self.eval_face_pairs(lb1,lb2)
                for postimprint_CommonFace in FacePairs:
                    BCType=postimprint_CommonFace.BCType

                    assert(lb1.Name is not None)
                    assert(lb2.Name is not None)
                    
                    if BCType is None: # BCType not otherwise set... insert default
                        BCType="TIE"
                        pass
                    elif bc_map is not None:
                        # apply user-supplied BC mapping
                        BCType=bc_map[BCType]
                        pass
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


    def save(self,cad_file_name,to_be_saved):
        step_writer=STEPControl_Writer()
        BodyNameList=[]
        BodyNameSet=set([])

        #for layerobj in to_be_saved:
        #    for layerbodyobj in layerobj.BodyList:
                
        for layerbodyobj in to_be_saved:
            step_writer.Transfer(layerbodyobj.Shape,STEPControl_ManifoldSolidBrep,True)
            
            assert(not(layerbodyobj.Name in BodyNameSet)) # verify nome layerbodyname reuse!
            
            BodyNameList.append(layerbodyobj.Name)
            
            BodyNameSet.add(layerbodyobj.Name)
            
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

    Layer2.Split(os.path.join("..","data","SplitLine2.csv"), pointTolerance)


    #delaminationlist = [ ]
    #delaminationlist = [ os.path.join("..","data","nasa-delam12-1.csv"), os.path.join("..","data","nasa-delam12-2.csv") ]

    MB.imprint_layers(Layer1, Layer2)
    MB.imprint_layers(Layer2, Layer3)
    #defaultBCType = 2
    #FAL = MB.adjacent_layer_boundary_conditions(Layer1,Layer2,defaultBCType)
    #MB.apply_delaminations(Layer1,Layer2,delaminationlist)


    step_writer=STEPControl_Writer()
    
    step_writer.Transfer(Layer1.BodyList[0].Shape,STEPControl_ManifoldSolidBrep,True)
    step_writer.Transfer(Layer2.BodyList[0].Shape,STEPControl_ManifoldSolidBrep,True)
    step_writer.Transfer(Layer2.BodyList[1].Shape,STEPControl_ManifoldSolidBrep,True)
    step_writer.Transfer(Layer3.BodyList[0].Shape,STEPControl_ManifoldSolidBrep,True)
    step_writer.Write(os.path.join("..","data","Layers.step"))

    pass
