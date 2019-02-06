import delamo.CADwrap
from delamo.api import DelamoModeler 
from delamo.api import Layer
from delamo.api import bond_layers

# NOTE: simple.laminate() needs some updating...

def laminate(initscript,controlpoints,layup,fixedpoint=None,forcesurfacepoint=None,delamination=None,delam_index=0):
    # initscript="abqparams.py
    #controlpoints="ControlPointsPlanar1.txt"
    DM=DelamoModeler.Initialize(facepointtolerancefactor=3.0,
                                       normaltolerance=100e-6)
    
    # LaminaSection represents the material properties, etc. of a single 
    # CFRP lamina
    LaminaSection=DM.assemblyinstrs.preexisting_variable("LaminaSection")

    # CohesiveInteraction represents the cohesive failure model between lamina
    CohesiveInteraction=DM.assemblyinstrs.preexisting_variable("CohesiveInteraction")
    
    # ContactInteraction represents the contact model used between 
    # delaminated layers
    ContactInteraction=DM.assemblyinstrs.preexisting_variable("ContactInteraction")

    # ApplyForceStep is a finite element modeling step, in particular the 
    # step where we apply forces. 
    
    ApplyForceStep=DM.assemblyinstrs.preexisting_variable("ApplyForceStep")

    # BendingJob is the job to run once the model is fully constructed. 
    # It is referenced through DM.runinstrs to that it doesn't get 
    # called until all of the setup instructions are complete. 
    BendingJob=DM.runinstrs.preexisting_variable("BendingJob")
    
    
    # Meshing parameters (used below)
    # -------------------------------
    MeshElemTypes=DM.assemblyinstrs.preexisting_variable("MeshElemTypes")
    ElemShape=DM.assemblyinstrs.preexisting_variable("ElemShape")
    ElemTetShape=DM.assemblyinstrs.preexisting_variable("ElemTetShape")
    ElemTechnique=DM.assemblyinstrs.preexisting_variable("ElemTechnique")
    meshsize=DM.assemblyinstrs.preexisting_variable("meshsize",unwrappedclass=float)
    finemeshsize=DM.assemblyinstrs.preexisting_variable("finemeshsize",unwrappedclass=float)
    
    # Note that it matters whether you define the preexisting variable
    # it on initinstrs, assemblyinstrs, contactinstrs etc. ... 
    # The instructions  used for the definition determine which context any
    # method calls would be run in.
    
    
    
    # This script then generates both a CAD file and a Python script.
    # The Python script can be run from Abaqus. It includes the 
    # initialization script referenced above, and also opens 
    # the CAD file and builds the model. 
    
    # CAD file to write
    cad_file_name = "CBlock.sat"
    
    # Python script file to write 
    script_to_generate="genabqscript.py"
    
    # Set delamination on or off
    is_delam = False
    
    # set stiffener on or off
    is_stiffener = True
    
    # Set layer thickness
    thickness = 0.199
    
    # Create the NURBS mold surface 
    # over which the lamina will be laid
    mold = delamo.CADwrap.NURBS()
    mold.degree_u(3)
    mold.degree_v(3)
    knot_vector_u = [0, 0, 0, 0, 1, 2, 3, 3, 3, 3]
    knot_vector_v = [0, 0, 0, 0, 1, 2, 3, 3, 3, 3]
    mold.knotvector_u(knot_vector_u)
    mold.knotvector_v(knot_vector_v)
    mold.read_ctrlpts(controlpoints)
    default_weight = 1.0
    weights = [default_weight for i in range(0, mold.ctrlpts_len())]
    mold.weights(weights)
    
    prev_layer=None
    nextlayernum=1
    for orientation in layup: 
        if prev_layer is None:
            layer=Layer.CreateFromMold(DM,mold,thickness,"Layer %d" % (nextlayernum),LaminaSection,orientation)
            pass
        else: 
            layer = Layer.CreateFromLayer(DM,prev_layer.gk_layer,delamo.CADwrap.OFFSET_DIRECTION,thickness,"Layer %d" % (nextlayernum), LaminaSection,orientation)
            
            pass
        layer.fe_layer_meshing.MeshSimple(MeshElemTypes,meshsize,ElemShape,ElemTechnique)
        
        if nextlayernum==1:
            layer.FixedFace(DM,ApplyForceStep,(layer.GetFaceFromPoint(DM,fixedpoint),))
            pass

        if prev_layer is not None:
            bond_layers(DM,prev_layer,layer,delamo.CADwrap.Delamination_COHESIVE,CohesiveInteraction,ContactInteraction)
            pass
            


        nextlayernum+=1
        prev_layer=layer

        pass

    # Add pressure applied to final layer 
    PressureMagnitude=-5e-2 # Units of MPa
    layer.Pressure(DM,ApplyForceStep,(layer.GetFaceFromPoint(DM,forcesurfacepoint),),PressureMagnitude,name="AppliedPressure")
    DM.Finalize(script_to_generate,cad_file_path_from_script)
