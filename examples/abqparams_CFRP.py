# This script is run in the Abaqus context, but wrapped proxies for variables defined here are available in the modeling context

# Create a new model database  (clear things out!)
abq.Mdb()

# Change ABAQUS working directory to the location of the script we are running
(scriptpath,scriptname)=os.path.split(getCurrentScriptPath())
if scriptpath != '':
    os.chdir(scriptpath)
    pass

# disable display (per scripting book, page 124)
abq.session.viewports['Viewport: 1'].setValues(displayedObject=None)

# Create model, give it a nicer name
abq.mdb.models.changeKey(fromName='Model-1',toName='FEModel')
FEModel=abq.mdb.models['FEModel']  

# Define LaminateAssembly -- shorthand for FEModel.rootAssembly 
LaminateAssembly=FEModel.rootAssembly



# Approximate material definitions for CFRP
# ** REMEMBER TO ALWAYS INCLUDE A DECIMAL POINT TO ENSURE FLOATING POINT
#    REPRESENTATION IN PYTHON 2.x **
CFRPDensity=1.5650   # N*s^2/mm
# Order of CFRPLaminaProps corresponds to abqC.LAMINA -- i.e. the ordering
# in the material editor for elastic 'Type=Lamina' is selected
#CFRPLaminaProps=( 1.415e11,  # E1, Pascals
#                  8.5e9, # E2, Pascals
#                  0.33, # nu12
#                  5.02e9, # G12, Pa
#                  5.02e9, # G13, Pa
#                  2.35e9) # G23, Pa 
CFRPEngineeringProps= ( 1.415e11/1.e6,  # E1, MPa
                        8.5e9/1.e6, # E2, MPa
                        8.5e9/1.e6, # E3, MPa
                        0.33, # nu12,
                        0.33, # nu13,
                        0.33, # nu23,  #!!! Pulled from thin air!!! Must recheck
                        5.02e9/1.e6, # G12, MPa
                        5.02e9/1.e6, # G13, MPa
                        2.35e9/1.e6) # G23,MPa  


#From Denizhan Yavas 4/7/16
#
#Yavas and Coker,  Composite Structures or Composites A or B, to be published
cohesiveslope=1.0e5 #  N/mm^3
tIC = 40.0 # MPa 
tIIC = 53.0 # MPa
GIC = .3753 # N/mm
GIIC = 1.4671 # N/mm
BK_Eta = 2.25 # mode mixing fitting factor for BK criterion

# From Song, Davila, and Rose, Guidelines and Parameter Selection for the Simulation of Progressive Delamination, 2008 Abaqus Users' Conference
E1=CFRPEngineeringProps[0]
E2=CFRPEngineeringProps[1]
Ne=3
M=1.0

l_cz_I = M*E2*GIC/(tIC**2)
l_cz_II = M*E2*GIIC/(tIIC**2)

l_e = min(l_cz_I/Ne,l_cz_II/Ne)

Ta = np.sqrt(E2*GIC/(Ne*l_e))
Sa = np.sqrt(E2*GIIC/(Ne*l_e))

Tbar = min(Ta,tIC)
Sbar = min(Sa,tIIC)

l_cz_I_expanded = M*E2*GIC/(Tbar**2)
l_cz_II_expanded = M*E2*GIIC/(Sbar**2)

l_e_expanded = min(l_cz_I_expanded/Ne,l_cz_II_expanded/Ne)

 
# Boundary condition params

DispBC_ZDispl=85.0 #80.0 # mm
meshsize=float(1.3) # mm
finemeshsize=float(0.25) # mm

# Create cohesive model
CohesiveInteraction=FEModel.ContactProperty("CohesiveInteraction")
CohesiveInteraction.CohesiveBehavior(
    defaultPenalties=abqC.OFF,
    table=((cohesiveslope,cohesiveslope,cohesiveslope),))
CohesiveInteraction.Damage(initTable=((tIC,tIIC,tIIC),),
                           useEvolution=abqC.ON,
                           evolutionType=abqC.ENERGY,
                           evolTable=((GIC,GIIC,GIIC),),
                           useMixedMode=abqC.ON,
                           mixedModeType=abqC.BK,
                           exponent=BK_Eta)

# Create contact model
ContactInteraction=FEModel.ContactProperty("ContactInteraction")
ContactInteraction.TangentialBehavior(formulation=abqC.FRICTIONLESS)
ContactInteraction.NormalBehavior(pressureOverclosure=abqC.HARD,
                                  allowSeparation=abqC.ON,
                                  constraintEnforcementMethod=abqC.DEFAULT)


# Define material properties for a single lamina
CFRPLaminaMat=FEModel.Material(name='CFRPLaminaMat') 
CFRPLaminaMat.Density(table=((CFRPDensity,),))
#CFRPLaminaMat.Elastic(type=abqC.LAMINA,table=(CFRPLaminaProps,))
CFRPLaminaMat.Elastic(type=abqC.ENGINEERING_CONSTANTS,table=(CFRPEngineeringProps,))

# Define section
LaminaSection=FEModel.HomogeneousSolidSection(name='LaminaSection',material=CFRPLaminaMat.name,thickness=None)


# Define force application step
ApplyForceStep=FEModel.StaticStep(name="ApplyForceStep",
                                        previous="Initial",
                                        nlgeom=abqC.ON,
                                        timeIncrementationMethod=AUTOMATIC,
                                        timePeriod=1.0,
                                        initialInc=0.05,
                                        minInc=1e-5,
                                        maxInc=0.05,
                                        description="Apply concentrated force")

# Alternative for Riks solver (maybe better convergence?)
# Standard solver sems to have much better convergence. -Bryan

## Define force application step
#ApplyForceStep=FEModel.StaticRiksStep(name="ApplyForceStep",
#                                        previous="Initial",
#                                        nlgeom=abqC.ON,   # may not have parameter.... seems to be limited to nonlinear
#                                        timeIncrementationMethod=AUTOMATIC,
#                                        #timePeriod=1.0,
#                                        #initialArcInc=0.05,
#                                        #minArcInc=1e-5,
#                                        #maxArcInc=0.1,
#                                        description="Apply concentrated force")




# Define field output
#FEModel.fieldOutputRequests.values()[0].setValues(variables=('S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'RF', 'CF', 'P', 'CSTRESS', 'CDISP')) # Should we define a new fieldOutput object or editing the existing default generated object

# Create a new output request with pressure as the requested output
# Could redefine the existing field output request to include pressure

#existing_request=FEModel.fieldOutputRequests.keys()[0]
#FEModel.fieldOutputRequests[existing_request].setValues(variables=('S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'RF', 'CF', 'P', 'CSTRESS', 'CDISP'))
FEModel.FieldOutputRequest(name='Output_Request', createStepName=ApplyForceStep.name, variables=('S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'RF', 'CF', 'P', 'CSTRESS', 'CDISP','DAMAGET','CSDMG')) 

# Instead of requesting output at every time step, request at every tenth of a time step
# With a time step of 1, this should give 10 frames of output
# Seems to affect convergence
for key in FEModel.fieldOutputRequests.keys():
    #FEModel.fieldOutputRequests[key].setValues(timeInterval=0.1)
    FEModel.fieldOutputRequests[key].setValues(frequency=1)
    pass


# Define mesh element type
# C3D10 for tetrahedral elements and C3D8R for hexahedral elements
# If ElemTechnique is Hex, only hex elements will be used
# There are mesh techniques where we can use a combination of elements, thus we specify
# element codes for different element types

MeshElemTypes=(#mesh.ElemType(elemCode=abqC.C3D10,
               #            elemLibrary=abqC.STANDARD,
               #            kinematicSplit=abqC.AVERAGE_STRAIN,
               #            secondOrderAccuracy=abqC.OFF, 
               #            hourglassControl=abqC.DEFAULT,
               #            distortionControl=abqC.DEFAULT),
               #mesh.ElemType(elemCode=abqC.C3D8,
               #            elemLibrary=abqC.STANDARD,
               #            kinematicSplit=abqC.AVERAGE_STRAIN,
               #            secondOrderAccuracy=abqC.OFF, 
               #            hourglassControl=abqC.DEFAULT,
               #              distortionControl=abqC.DEFAULT),
               mesh.ElemType(elemCode=abqC.C3D20, # quadratic mesh element
                           elemLibrary=abqC.STANDARD,
                           kinematicSplit=abqC.AVERAGE_STRAIN,
                           secondOrderAccuracy=abqC.OFF, 
                           hourglassControl=abqC.DEFAULT,
                             distortionControl=abqC.DEFAULT),
               mesh.ElemType(elemCode=abqC.C3D4)), # Tet element


# Possible values are: TET, HEX, HEX_DOMINATED 
# Default is HEX
# Hex should be more accurate, while Tet is more stable
# We are using 2nd order Tet to help with accuracy

#ElemShape=abqC.HEX_DOMINATED
#ElemTetShape=abqC.TET

# Possible values are: FREE, STRUCTURED, SWEEP 
# If you want Abaqus to define the meshing technique internally
# leave this value as None
#ElemTechnique=abqC.STRUCTURED
#ElemTechnique=abqC.FREE
FreeElemTechnique=abqC.FREE
ElemTechnique=abqC.SYSTEM_ASSIGN


LaminateAssembly.DatumCsysByDefault(abqC.CARTESIAN) # ?



# Create job
BendingJob=abq.mdb.Job(name="BendingJob",
                       model=FEModel.name,
                       description="Bending analysis job",
                       type=abqC.ANALYSIS,
                       atTime=None,
                       waitMinutes=0,
                       waitHours=0,
                       queue=None,
                       memory=90,
                       memoryUnits=abqC.PERCENTAGE,
                       getMemoryFromAnalysis=True,
                       explicitPrecision=abqC.SINGLE, 
                       nodalOutputPrecision=abqC.SINGLE,
                       echoPrint=abqC.OFF, modelPrint=abqC.OFF, 
                       contactPrint=abqC.OFF, historyPrint=abqC.OFF,
                       userSubroutine='', scratch='', 
                       resultsFormat=abqC.ODB,
                       multiprocessingMode=abqC.THREADS,
                       numCpus=4, numDomains=4, numGPUs=0)

