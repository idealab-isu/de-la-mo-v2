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

from OCC.STEPControl import STEPControl_Reader
from OCC.STEPControl import STEPControl_Writer
from OCC.STEPControl import STEPControl_ManifoldSolidBrep
from OCC.STEPControl import STEPControl_Writer,STEPControl_ShellBasedSurfaceModel,STEPControl_GeometricCurveSet

from layer import LayerBody,LayerBodyFace

import loaders
import layer




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
        OffsetDist = 100.0*self.PointTolerance
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

        FaceListTotal1 = layerbody1.FaceListOrig + layerbody1.FaceListOffset + layerbody1.FaceListSide
        FaceListTotal2 = layerbody2.FaceListOrig + layerbody2.FaceListOffset + layerbody2.FaceListSide
        
        CommonFaces = {}
        for face1 in FaceListTotal1:
            if face1 in FaceListTotal2:
                face2=FaceListTotal2[FaceListTotal2.index(face1)]

                # Note that face1 and face2 pass equality (== operator) because
                # equality operator identifices equivalent geometry. But the face1 and
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

            # step_writer2=STEPControl_Writer()
            # step_writer2.Transfer(WireShape,STEPControl_GeometricCurveSet,True)
            # step_writer2.Write("../data/Wire.STEP")
            #
            # sys.modules["__main__"].__dict__.update(globals())
            # sys.modules["__main__"].__dict__.update(locals())
            # raise ValueError("Break")

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
        
        # step_writer2=STEPControl_Writer()
        # step_writer2.Transfer(SideShape,STEPControl_ShellBasedSurfaceModel,True)
        # #step_writer2.Transfer(layerbody.Shape, STEPControl_ManifoldSolidBrep, True)
        # #step_writer2.Transfer(layerbody2.Shape, STEPControl_ManifoldSolidBrep, True)
        # #step_writer2.Transfer(SplitFace,STEPControl_ShellBasedSurfaceModel,True)
        # step_writer2.Write("../data/allShapes.STEP")

        split_face_exp=TopExp_Explorer(SplitFace,TopAbs_FACE)
        # Iterate over all faces
        numsplitfaces = 0
        split_face_shapes=[]
        split_layerbodyfaces=[]
        while split_face_exp.More():
            split_face_shape=split_face_exp.Current()
            split_face_shapes.append(split_face_shape)

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

        print("Number of split faces %d"%(numsplitfaces))

        # split_face_shapes now should have two or more faces (TopoDS_Shape of type Face
        # Need to create a new layerbody with the original layerbodyfaces except for this one, and two new layerbodyfaces


        # Create a new layerbody with the split face
        DelamLayerBody = copy.copy(layerbody)  # non-deep copy
        # in the copy we are initializing, assign the face lists with
        # the face we are replacing removed, and with the new faces
        # we have created, added.

        DelamLayerBody.FaceListOrig = layerbody.FaceListOrig
        if layerbodyface in DelamLayerBody.FaceListOrig:
            # Remove original face
            del DelamLayerBody.FaceListOrig[DelamLayerBody.FaceListOrig.index(layerbodyface)]
            # Add new faces
            for split_layerbodyface in split_layerbodyfaces:
                DelamLayerBody.FaceListOffset.append(split_layerbodyface)
                pass
            pass

        DelamLayerBody.FaceListOffset = layerbody.FaceListOffset
        if layerbodyface in DelamLayerBody.FaceListOffset:
            # Remove offset face
            del DelamLayerBody.FaceListOffset[DelamLayerBody.FaceListOffset.index(layerbodyface)]
            # Add new faces
            for split_layerbodyface in split_layerbodyfaces:
                DelamLayerBody.FaceListOffset.append(split_layerbodyface)
                pass
            pass

        

        # !!!***  Need to set BCTType on each generate LayerBodyFace !!!***
        
        # (Could also do similar process on side faces, but how could we ever get a delamination on the side faces???)

        #sys.modules["__main__"].__dict__.update(globals())
        #sys.modules["__main__"].__dict__.update(locals())
        #raise ValueError("Break")

        DelamLayerBody._Initializing_Layerbody_Construct_Shape()
        

        #step_writer = STEPControl_Writer()
        #step_writer.Transfer(delamSolidShape, STEPControl_ManifoldSolidBrep, True)
        #step_writer.Write("../Data/Layers.step")

        # Also need to repeat the process for the other face

        # (so actual content of this function should be abstracted into
        # a new function we can call twice). 

        return DelamLayerBody

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
                
                CommonFaces=self.eval_face_pairs(lb1,lb2)

                for CommonFace in CommonFaces:
                    replacement_lb1=self.imprint_delaminations(lb1,CommonFace,delaminationlist)
                    replacement_lb2=self.imprint_delaminations(lb2,CommonFaces[CommonFace],delaminationlist)
                    layer1.BodyList[lb1cnt]=replacement_lb1
                    layer2.BodyList[lb2cnt]=replacement_lb2
                    pass

                pass
            pass
        pass
    
    def adjacent_layers(self,layer1,layer2,defaultbc,bc_map=None):
        """ Once adjacent_layers() is called, the LayerBody's in EITHER layer can't be split
        anymore -- because then they might need new names,
        and the return values contain the layer body names that will be used
        to apply the boundary conditions"""
        
        FAL = [] # Face Adjacency List
        
        for lb1 in layer1.BodyList:
            for lb2 in layer2.BodyList:
                FacePairs=self.eval_face_pairs(lb1,lb2)
                for postimprint_CommonFace in FacePairs:
                    BCType=postimprint_CommonFace.BCType

                    assert(lb1.Name is not None)
                    assert(lb2.Name is not None)
                    
                    if BCType is None: # BCType not otherwise set... insert default
                        BCType=defaultbc
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
    
    Mold = layer.LayerMold.FromFile(os.path.join("..","data","CurvedMold1.STEP"))
    #Mold = layer.LayerMold.FromFile(os.path.join("..","data","FlatMold3.STEP"))
    Layer1=layer.Layer.CreateFromMold("Layer 1",Mold,2.0,"OFFSET",1e-6)
    #Layer1=layer.Layer.CreateFromMold("Layer 1",Mold,2.0,"ORIG",1e-6)
    Layer2=layer.Layer.CreateFromMold("Layer 2",Layer1.OffsetMold(),2.0,"OFFSET",1e-6)
    #Layer2=layer.Layer.CreateFromMold("Layer 2",Mold,2.0,"OFFSET",1e-6)

    #delaminationlist = [ ]
    delaminationlist = [ os.path.join("..","data","nasa-delam12-1.csv"), os.path.join("..","data","nasa-delam12-2.csv") ]

    #defaultBCType = 2
    #FAL = MB.adjacent_layers(Layer1,Layer2,defaultBCType)
    MB.apply_delaminations(Layer1,Layer2,delaminationlist)


    step_writer=STEPControl_Writer()
    
    step_writer.Transfer(Layer1.BodyList[0].Shape,STEPControl_ManifoldSolidBrep,True)
    step_writer.Transfer(Layer2.BodyList[0].Shape,STEPControl_ManifoldSolidBrep,True)
    step_writer.Write(os.path.join("..","data","Layers.step"))

    pass
