import sys
import os.path
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
from OCC.TopTools import TopTools_ListIteratorOfListOfShape
from OCC.GeomLProp import GeomLProp_SLProps
from OCC.gp import gp_Pnt2d
from OCC.gp import gp_Vec
from OCC.gp import gp_Dir
from OCC.gp import gp_Pnt

from OCC.STEPControl import STEPControl_Reader
from OCC.STEPControl import STEPControl_Writer
from OCC.STEPControl import STEPControl_ManifoldSolidBrep

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

    def process_delamination(self,layerbody1,facebody1,layerbody2,facebody2,delam_outline):
        """Given a first layerbody and corresponding face, and a second 
        layerbody and corresponding face, and a delamination outline 
        (loop of 3D coordinates, hopefully projected onto the faces): 
           1. Identify the regions of facebody1/2 that are inside delam_outline
           2. Create an offset curve a distance of self.GapWidth inside delam_outline
           3. Imprint both delam_outline and the offset curve
           4. Identify the regions of facebody1/2 that are between delam_outline and the offset
              curve to have "NOMODEL" BCType unless they already had "CONTACT" BCType 
           5. Identify the regions of facebody1/2 that are inside the offset curve to have "CONTACT"
              BCType
           6. Generate new replacements for layerbody1 and layerbody2 with newly constructed
              imprinted faces, marked as determined above. 
           7. The replaced layerbody1 and layerbody2 may have facebody1 and facebody2 
              replaced/subdivided, but other faces in layerbodies1/2 should remain unchanged. 
           8. Return (replacement_layerbody1, replacement_layerbody2)"""

        # NOTE: May need additional parameters (adjacent surfaces or faces?) to do the
        # delam_outline offset curve?

        # NOTE: When regenerating layerbodies, do NOT give them new names unless they are being
        # split (which they aren't from this function)






        
        raise NotImplementedError()

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
                    for delamination in delaminationlist:
                        (replacement_lb1,replacement_lb2)=self.process_delamination(lb1,CommonFace,lb2,CommonFaces[CommonFace],delamination)
                        layer1.BodyList[lb1cnt]=replacement_lb1
                        layer2.BodyList[lb2cnt]=replacement_lb2
                        pass
                    pass

                pass
            pass
        pass
    
    def adjacent_layers(self,layer1,layer2,defaultbc):
        """ Once adjacent_layers() is called, the LayerBody's in EITHER layer can't be split
        anymore -- because then they might need new names,
        and the return values contain the layer body names that will be used
        to apply the boundary conditions"""
        
        FAL = [] # Face Adjacency List
        
        for lb1 in layer1.BodyList:
            for lb2 in layer1.BodyList:
                FacePairs=self.eval_face_pairs(lb1,lb2)
                for postimprint_CommonFace in FacePairs:
                    BCType=postimprint_CommonFace.BCType

                    if BCType is None: # BCType not otherwise set... insert default
                        BCType=defaultbc
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
    
    pass


if __name__=="__main__":
    MB=OCCModelBuilder(PointTolerance=1e-5,NormalTolerance=1e-6)
    
    Mold = layer.LayerMold.FromFile(os.path.join("..","data","CurvedMold1.STEP"))
    Layer1=layer.Layer.CreateFromMold("Layer 1",Mold,2.0,"OFFSET",1e-6)
    Layer2=layer.Layer.CreateFromMold("Layer 2",Layer1.OffsetMold(),2.0,"OFFSET",1e-6)

    delaminationlist = [ ]
    #delaminationlist = [ os.path.join("..","data","nasa-delam-12-1.csv") ]

    defaultBCType = 2
    FAL = MB.adjacent_layers(Layer1,Layer2,defaultBCType)
    #MB.apply_delaminations(Layer1,Layer2,delaminationlist)


    step_writer=STEPControl_Writer()
    
    step_writer.Transfer(Layer1.BodyList[0].Shape,STEPControl_ManifoldSolidBrep,True)
    step_writer.Transfer(Layer2.BodyList[0].Shape,STEPControl_ManifoldSolidBrep,True)
    step_writer.Write(os.path.join("..","Data","Layers.step"))

    pass
