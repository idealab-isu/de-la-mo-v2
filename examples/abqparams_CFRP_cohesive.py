# Define material properties for cohesive layer
cohesive_damage_initiation_strain_normal=1e-5
cohesive_damage_initiation_strain_shear=1e-5
cohesive_fracture_energy=.3753   # Same as GIC (???)

CohesiveE=8.5e9/1.e6  # E, MPa
CohesiveG=3.2e9/1.e6 # G, MPa

CFRPCohesiveMat=FEModel.Material(name='CFRPCohesiveMat') 
CFRPCohesiveMat.Density(table=((CFRPDensity,),))
#CFRPCohesiveMat.Elastic(type=abqC.LAMINA,table=(CFRPLaminaProps,))
#CFRPCohesiveMat.Elastic(type=abqC.ENGINEERING_CONSTANTS,table=(CFRPEngineeringProps,))
CFRPCohesiveMat.Elastic(type=abqC.TRACTION,table=((CohesiveE, CohesiveG, CohesiveG),))
CFRPCohesiveMat.MaxeDamageInitiation(table=((
    cohesive_damage_initiation_strain_normal,
    cohesive_damage_initiation_strain_shear,
    cohesive_damage_initiation_strain_shear), ))
CFRPCohesiveMat.maxeDamageInitiation.DamageEvolution(
    type=abqC.ENERGY, table=((cohesive_fracture_energy, ), ))

CohesiveSection=FEModel.CohesiveSection(name='CohesiveSection', 
    material=CFRPCohesiveMat.name, response=abqC.TRACTION_SEPARATION, 
    outOfPlaneThickness=None)


# Update output request fields -- Add "STATUS" field, makes cohesive
# elements disappear once they fail
FEModel.fieldOutputRequests['Output_Request'].setValues(
    variables=('S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'RF', 'CF', 'P', 
               'CSTRESS', 'CDISP', 'DAMAGET','CSDMG','STATUS'))
