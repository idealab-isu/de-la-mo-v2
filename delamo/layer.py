import sys
import os
import os.path

import numpy as np



from OCC.TopoDS import topods

from OCC.TopoDS import TopoDS_Face
from OCC.TopoDS import TopoDS_Shape
from OCC.TopoDS import TopoDS_Compound
from OCC.TopoDS import TopoDS_Wire
from OCC.TopoDS import TopoDS_Vertex
from OCC.TopoDS import TopoDS_Shell
from OCC.TopoDS import topods_Shell
from OCC.TopoDS import topods_Face
from OCC.BRep import BRep_Builder
from OCC.BRep import BRep_Tool
from OCC.BRepExtrema import BRepExtrema_DistShapeShape
from OCC import BRepLib
from OCC import BRepOffsetAPI
from OCC import BRepOffset
from OCC import BRepBuilderAPI
from OCC.BRepClass import BRepClass_FacePassiveClassifier
from OCC.ShapeAnalysis import ShapeAnalysis_FreeBoundsProperties
from OCC.BRepTools import breptools_Read
from OCC.TopExp import TopExp_Explorer
from OCC.TopAbs import TopAbs_FACE
from OCC.TopAbs import TopAbs_VERTEX
from OCC.TopAbs import TopAbs_EDGE
from OCC.TopAbs import TopAbs_SHELL
from OCC.TopAbs import TopAbs_FORWARD
from OCC.TopAbs import TopAbs_REVERSED
from OCC.GeomAbs import GeomAbs_Arc
from OCC.TopTools import TopTools_ListIteratorOfListOfShape
from OCC.GeomLProp import GeomLProp_SLProps
from OCC.gp import gp_Pnt2d
from OCC.gp import gp_Vec
from OCC.gp import gp_Dir
from OCC.gp import gp_Pnt

from OCC.STEPControl import STEPControl_Reader
from OCC.STEPControl import STEPControl_Writer,STEPControl_ShellBasedSurfaceModel
from OCC.STEPControl import STEPControl_Writer,STEPControl_ManifoldSolidBrep
from OCC.STEPControl import STEPControl_Writer,STEPControl_GeometricCurveSet
from OCC.IGESControl import IGESControl_Reader
from OCC.IFSelect import IFSelect_RetDone, IFSelect_ItemsByEntity

import loaders

def FindOCCPointNormal(Face, refParPoint, OrigPointTolerance):
    """ Given a face, find and return a (point, normal, parPoint)
    that uniquely identifies the face. The point must be 
    significantly farther than tolerance from any edge"""

    # Evaluate reference parametric point
    faceSurface = BRep_Tool().Surface(Face)
    refPointProps = GeomLProp_SLProps(faceSurface, refParPoint[0], refParPoint[1], 1, OrigPointTolerance)
    faceNormal = gp_Vec(refPointProps.Normal())
    facePoint = refPointProps.Value()

    # Check if original reference point is inside tha Face
    C=BRepClass_FacePassiveClassifier()
    C.Perform(Face, facePoint, OrigPointTolerance)
    if (C.State()==TopAbs_ON):
        print('Point in Face')


    sys.modules["__main__"].__dict__.update(globals())
    sys.modules["__main__"].__dict__.update(locals())
    raise ValueError("Break")

    # Find the closest point by this method:
    # https://www.opencascade.com/content/closest-point-step-object
    # DistCalc = BRepExtrema_DistShapeShape(Face, OrigDirPointVertex)
    # DistCalc.Perform()
    # ThisDist = DistCalc.Value()
    #
    # if DistCalc.NbSolution() > 0 and ThisDist < ClosestDist:
    #     ClosestDist = ThisDist
    #     (ClosestU, ClosestV) = DistCalc.ParOnFaceS1(1)  # Evaluate (u,v) coordinates on this face of closest point
    #
    #     # Here is an example of how to extract
    #     # a normal vector given (U,V) coordinates
    #     ThisFaceSurface = BRep_Tool().Surface(topods_Face(FaceExp.Current()))
    #     ClosestUVProps = GeomLProp_SLProps(ThisFaceSurface, ClosestU, ClosestV, 1, OrigDirTolerance)
    #     ClosestNormalVec = gp_Vec(ClosestUVProps.Normal())
    #     ClosestNormal = np.array((ClosestNormalVec.X(), ClosestNormalVec.Y(), ClosestNormalVec.Z()), dtype='d')
    #     if FaceExp.Current().Orientation == TopAbs_REVERSED:
    #         # Face is reversed from underlying surface -> we need to flip the normal
    #         ClosestNormal = -ClosestNormal
    #         pass
    #     pass
    # pass

    return ("FindOCCPointNormal() not implemented","FindOCCPointNormal() not implemented")



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

    def OffsetMold(self):
        """Return a LayerMold based on the OFFSET side of this layer"""
        return LayerMold.FromFaceLists([ Body.FaceListOffset for Body in self.BodyList ])

    def OrigMold(self):
        """Return a LayerMold based on the ORIG side of this layer.
        (Note that the OFFSET direction of this LayerMold points away
        from this layer!) """
        return LayerMold.FromFaceLists([ Body.FaceListOrig for Body in self.BodyList ])

    
    @classmethod
    def CreateFromMold(cls,Name,Mold,Thickness,Direction,Tolerance):
        """Create a layer from a LayerMold."""

        if Direction=="OFFSET":
            OffsetDist=Thickness
            pass
        elif Direction=="ORIG":
            OffsetDist=-Thickness
            pass
        else:
            raise ValueError("Invalid Direction: %s" % (str(Direction)))

        
        
        mkOffset = BRepOffsetAPI.BRepOffsetAPI_MakeOffsetShape(Mold.Shape, OffsetDist, Tolerance,
                                                               BRepOffset.BRepOffset_Skin,
                                                               False, False,
                                                               GeomAbs_Arc)
        assert (mkOffset.IsDone())
        OffsetShell = mkOffset.Shape()

        # Build list OffsetFaces of faces in OffsetShell
        OffsetShellFacesExp = TopExp_Explorer(OffsetShell,TopAbs_FACE)
        OffsetFaces = []

        while OffsetShellFacesExp.More():
            OffsetFaces.append(topods_Face(OffsetShellFacesExp.Current()))
            OffsetShellFacesExp.Next()
            pass
        

        # Check for outline of original shape
        FreeCheck = ShapeAnalysis_FreeBoundsProperties(Mold.Shape)
        FreeCheck.Perform()
        assert (FreeCheck.NbClosedFreeBounds() >= 1)
        
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

        #sys.modules["__main__"].__dict__.update(globals())
        #sys.modules["__main__"].__dict__.update(locals())
        #raise ValueError("Break")
        solidShape = solidMaker.Solid()
        
        if not BRepLib.breplib_OrientClosedSolid(solidShape):
            raise ValueError("Solid maker did not yield a closed solid")
        # We successfully got a closed solid

        
        #return [layerSolid, offsetSurface]
        NewLayer= cls(Name=Name,
                      Type="LAMINA",
                      Direction=Direction,
                      Thickness=Thickness,
                      RefMold=Mold)


        NewLayerBody = LayerBody(Owner=NewLayer,
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
                MoldFaceSurf=BRep_Tool.Surface(MoldFace.Face)
                if MoldFaceSurf==FaceSurf: # Same underlying surface
                    if MatchedInMold:
                        raise ValueError("Same surface matched twice in mold (!?)")
                    MatchedInMold=True
                    pass
                pass


            # Search for this face in the offset surface
            MatchedInOffset=False
            for OffsetFace in OffsetFaces: # Iterate over LayerBodyFaces in Mold
                OffsetFaceSurf=BRep_Tool.Surface(OffsetFace)
                if OffsetFaceSurf==FaceSurf: # Same underlying surface
                    if MatchedInOffset:
                        raise ValueError("Same surface matched twice in offset (!?)")
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


class LayerBody(object):
    """ The LayerBody is a solid, defined as a boundary representation
    from a shell which in turn consists of many faces. It represents
    a portion of a layer or stiffener.

    Value semantics: The LayerBody is immutable once created. If 
    there is a need to change it, create a new LayerBody, being
    sure to deep copy anything being changed, and remove the 
    old LayerBod(ies) from the owning Layer. The FaceLists will 
    in general need to be re-created with their owners pointing
    to the new LayerBody
"""
    Name = None
    FaceListOrig=None # Faces on the "ORIG" side (LayerBodyFace)
    FaceListOffset=None # Faces on the "OFFSET" side (LayerBodyFace)
    FaceListSide=None # Facees on the "SIDE" sides (LayerBodyFace)
    Owner = None # reference to Layer object of which this LayerBody MAY be a part
    #Mold = None  # LayerMold that was used to generate this body SHOULD WE REALLY HAVE THIS ATTRIBUTE? ALL OF THE RELEVANT INFO SHOULD BE IN FaceListOrig
    Shape = None # OpenCascade TopoDS Shape which is a Solid which consists of a single closed shell which consists of multiple faces.
    
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
                raise ValueError("No closest point on loaded shape to vertex at %s" % (str(OrigDirPoint)))

            
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
        
        return cls.FromShell(Shell)

    
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
                  # layerbody. Direction is arbitrary if this face is not
                  # part of a LayerBody
    ParPoint = None  # Numpy Array representing parametric coordinates of a point on the face
                     # This is initialized to (0.05,0.05) and passed to OCCFindPointNormal, where
                     # it is updated if the reference parameter is outside the trimmed region
    # There had been an "InitialSurface" boolean... I think this info
    # should be covered by Direction=="ORIG"
    Direction = None # "ORIG", "OFFSET", "SIDE", or "NODIR", representing the side of the Owning LayerBody which corresponds to this Face

    Owner = None # If this LayerBodyFace may be part of a LayerBody, this is the LayerBody of which it might be a part. 
    
    BCType = None # Formerly DelaminationType: None, "NODELAM" "NOMODEL", "COHESIVE", "CONTACT" or, "TIE"
    # MatchingFace = None # Formerly SurfPair, This would be the matching LayerBodyFace in the adjacent (or non-adjacent)
    # layer or stiffener, assigned by adjacent_layers() ... To maintain immutability of these objects suggest that
    # instead of this attribute, have adjacent_layers() create a dictionary
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
    @classmethod
    def FromOCC(cls,Face,Direction,Owner=None,BCType=None,IsPointingInside=None):
        """ Create a LayerBodyFace from a TopoDS_Face. 
        IsPointingInside should be None (if this LayerBodyFace is 
        not actually part of a LayerBody, or if the OCC-Evaluated normal
        is guaranteed to be outward-facing) or a function/method/lambda 
        that, given a Normal, returns whether  that normal is pointing inside
        the LayerBody."""

        RefParPoint = np.array([0.01,0.01])
        (Point,Normal,ParPoint)=FindOCCPointNormal(Face,RefParPoint,OrigPointTolerance=1e-5)

        # For a LayerBodyFace inside a LayerBody, the normal
        # should be pointing outward. 
        if IsPointingInside is not None and IsPointingInside(Normal):
            Normal=-Normal
            pass
        
        return cls(Face=Face,
                   Point=Point,
                   Normal=Normal,
                   Direction=Direction,
                   Owner=Owner,
                   BCType=BCType)
        

    pass

    
if __name__=="__main__":
    Mold = LayerMold.FromFile("../Data/CurvedMold1.STEP")
    Layer1=Layer.CreateFromMold("Layer 1",Mold,2.0,"OFFSET",1e-6)
    Layer2=Layer.CreateFromMold("Layer 2",Layer1.OffsetMold(),2.0,"OFFSET",1e-6)

    step_writer=STEPControl_Writer()
    
    step_writer.Transfer(Layer1.BodyList[0].Shape,STEPControl_ManifoldSolidBrep,True)
    step_writer.Transfer(Layer2.BodyList[0].Shape,STEPControl_ManifoldSolidBrep,True)
    step_writer.Write("../Data/Layers.step")

    pass
