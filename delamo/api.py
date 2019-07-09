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

"""
Class hierarchy
---------------
 * DelamoModeler stores the various behind-the-scenes information required for model building
 * Part represents a CAD and finite element part
 * LayerPart is a particular kind of part
 * Assembly is a collection of parts and/or subassemblies
 * Layer is a particular type of Assembly, containing one or more LayerParts, representing a single layer
 * LaminaContact defines the contact/adhesion/cohesion between two lamina layeers
 * shell_solid_coupling stores the parameters needed for shell-to-solid coupling
 * CoordSys is an abstract class representing a coordinate system
 * SimpleCoordSys is an implementation of CoordSys that is a fixed Cartesian frame
"""

import collections
import sys
import os
import os.path
import ast
import numpy as np
import traceback

from delamo.codegen import codestore
from delamo.codegen import namedbinding_wrapper

from delamo.abqscript import write_abq_script
from delamo.abqscript import _capture_assignments_as_variables
from delamo.OCCModelBuilder import OCCModelBuilder
from delamo.layer import Layer as OCCLayer
import delamo

from delamo.autofiber.generator import AutoFiber


#import delamo.CADwrap

# Abaqus constants that need to be globally
# accessible can be wrapped by the namedbinding_wrapper
# (note that this won't work for anything that calls
# or accepts parameters from the primary codegen
# wrapper)

# if you want to access abqC from a user-defined module,
# use "from delamo.api import abqC"
abqC_wrapper=namedbinding_wrapper("abqC")
abqC=abqC_wrapper.wrapped






class DelamoModeler(object):
    """ Represents the various behind-the-scenes information needed by the model builder"""
    
    # codestores for code generation
    initinstrs=None
    """Codestore of ABAQUS initialization instructions"""
    
    assemblyinstrs=None
    """Codestore of ABAQUS assembly instructions"""

    bcinstrs=None
    """Codestore of ABAQUS boundary condition instructions"""

    meshinstrs=None
    """Codestore of ABAQUS meshing instructions"""

    fiberinstrs=None
    """Codestore of ABAQUS fiber orientation instructions"""

    runinstrs=None
    """Codestore of ABAQUS running instructions"""

    globals=None
    """Global variable dictionary providng context for main script and imported functions/methods, including functions from abqfuncs_*.py"""
    
    # Finite element wrapped classes/functions
    abq=None
    """Proxy object represents the "abq" module, for execution with initinstrs"""
    abq_assembly=None  # executed in assemblyinstrs
    """Proxy object represents the "abq" module, for execution with assemblyinstrs"""
    abqC=None
    """Proxy object represents the abaqusConstants module"""
    regionToolset=None
    """Proxy object represents the ABAQUS regionToolset module, for execution with assemblyinstrs"""
    
    mesh=None
    """Proxy object represents the ABAQUS "mesh" module, for execution with meshinstrs"""
    section=None
    """Proxy object represents the ABAQUS "section" module, for execution with initinstrs"""
    assembly=None
    """Proxy object represents the ABAQUS "assembly" module, for execution with assemblyinstrs"""
    connector=None
    """Proxy object represents the ABAQUS "connector" module, for execution with bcinstrs"""
    
    InitDict=None
    """A proxy object for a Python dictionary that is configured for execution with initinstrs"""

    autofiber=None
    """Abaqus reference to the auto fiber class"""

    # Finite element parameters presumed to be defined
    # in the abqparams initialization scripts
    FEModel=None        # AssemblyInstrs
    """Variable representing abq.mdb.models['FEModel'] which should be assigned in your initialization script"""
    
    LaminateAssembly=None # AssemblyInstrs
    """Variable representing abq.mdb.models['FEModel'].rootAssembly which should be assigned in your initialization script"""

    stepfile=None  # AssemblyInstrs
    """Variable wrapped for ABAQUS access representing name of generated CAD file to load"""
    
    # Body number database
    BodyNumDB_build=None
    """ModelBuilder BodyNumber database proxy used in build phase"""
    BodyNumDB=None
    """Actual body number database: dictionary indexed by body name of body number"""
    #BodyDB_build=None # BodyDB in build phase
    #BodyDB=None # dictionary indexed by layer name of list of body names


    # Geometry kernel parameters
    modelbuilder=None
    """delamo.OCCModelBuilder.OCCModelBuilder() object."""


    abqpointtolerance=None
    """Default point tolerance for finite element face and edge identification"""
    normaltolerance=None
    """Default normal tolerance for finite element face identification"""
    tangenttolerance=None
    """Default tangent tolerance for finite element edge identification"""
    

    # to_be_saved and uniquenumber are initialized in the 
    # constructor. Everything else is initialized in Initialize()
    to_be_saved=None
    """List of MBBody objects to be saved into the output CAD file"""
    to_be_saved_dependencies=None
    """MBBody objects in to_be_saved may be dependent on the existence of their parent objects. This is a list, a place to store such references so their parents don't get destroyed """
    
    uniquenumber=None
    """value for next unique number from get_unique()"""


    def __init__(self,**kwargs):
        self.uniquenumber=0

        # Put the CAD layers into a std::vector implementation
        # This  works like a regular Python list due to the appropriate definitions in Swig and C++ code
        self.to_be_saved=[] 
        self.to_be_saved_dependencies=[]

        for argname in kwargs:
            if not hasattr(self,argname):
                raise KeyError(argname)
            
            setattr(self,argname,kwargs[argname])
            pass
        pass


    def get_unique(self):
        """Get a unique number. The numbers from this method of a particular DelamoModeler object will never repeat."""
        self.uniquenumber+=1
        return self.uniquenumber-1


    def add_abaqus_code(self,instrs,func_code,globals):
        """Add the specified ABAQUS code func_code to the specified instrs codestore. Proxies for the left sides of any 
        assignment statements found in func_code will be added to the dictionary globals"""
        instrs.append_code(func_code)
        syntax_tree=ast.parse(func_code)
        newglobals={}
        _capture_assignments_as_variables(syntax_tree.body,instrs,globals,newglobals)
        globals.update(newglobals)

        pass
        # return instrs.wrap_function(None,func_name)

    def exec_abaqus_script(self,execinstrs,accessinstrs,filename,globals):
        """Add the ABAQUS code from the specified file to the specified execinstrs codestore. Proxies for the left sides of any 
        assignment statements found in func_code will be added to the dictionary globals, for access in the accessinstrs codestore."""
        fh=open(filename,"r")
        filetext=fh.read()
        fh.close()

        execinstrs.append_code(filetext)
        syntax_tree=ast.parse(filetext,filename=filename)

        newglobals={}
        _capture_assignments_as_variables(syntax_tree.body,accessinstrs,globals,newglobals)
        globals.update(newglobals)
        pass
        
    def abaqus_init_script(self,filename,globals):
        """Convenience routine for running a script in initinstrs
        Methods or attribute accesses for variables defined in this script will 
        be referenced in assemblyinstrs """
        self.exec_abaqus_script(self.initinstrs,self.assemblyinstrs,filename,globals)
        pass

    def abaqus_assembly_functions(self,filename,globals):
        """ Convenience routine for running a script in initinstrs
        Methods or attribute accesses for variables defined in this script will 
        be referenced in accessinstrs (for the moment this is identical to abaqus_init_script() """
        self.exec_abaqus_script(self.initinstrs,self.assemblyinstrs,filename,globals)
        pass

    def abaqus_bc_functions(self,filename,globals):
        """ Convenience routine for running a script in initinstrs
        Methods or attribute accesses for variables defined in this script will 
        be referenced in bcinstrs """
        self.exec_abaqus_script(self.initinstrs,self.bcinstrs,filename,globals)
        pass

    def abaqus_mesh_functions(self,filename,globals):
        """ Convenience routine for running a script in initinstrs
        Methods or attribute accesses for variables defined in this script will 
        be referenced in meshinstrs """
        self.exec_abaqus_script(self.initinstrs,self.meshinstrs,filename,globals)
        pass

    def abaqus_fiber_functions(self,filename,globals):
        """ Convenience routine for running a script in initinstrs
        Methods or attribute accesses for variables defined in this script will
        be referenced in fiberinstrs """
        self.exec_abaqus_script(self.initinstrs,self.fiberinstrs,filename,globals)
        pass

    def abaqus_run_functions(self,filename,globals):
        """ Convenience routine for running a script in initinstrs
        Methods or attribute accesses for variables defined in this script will 
        be referenced in runinstrs """
        self.exec_abaqus_script(self.initinstrs,self.runinstrs,filename,globals)
        pass


    def abaqus_functions(self,globals):
        """Add the delamo builtin ABAQUS scripts abqfuncs_assembly.py, abqfuncs_bc.py, abqfuncs_mesh.py, and abqfuncs_run.py to the 
        DelamoModeler codestores"""
        (api_dir,api_file) = os.path.split(self.__class__.__file__)

        abaqus_assembly_functions(os.path.join(api_dir,"abqfuncs_assembly.py"),globals)
        abaqus_bc_functions(os.path.join(api_dir,"abqfuncs_bc.py"),globals)
        abaqus_mesh_functions(os.path.join(api_dir,"abqfuncs_mesh.py"),globals)
        abaqus_run_functions(os.path.join(api_dir,"abqfuncs_run.py"),globals)
        pass
    
    @classmethod
    def Initialize(cls,globals,
                   pointtolerance=1e-7, # Geometry Kernel tolerance (mm)
                   pointtolerancefactor=100.0, # How much larger the tolerance should be in Abaqus
                   normaltolerance=100e-4,
                   tangenttolerance=100e-4,
                   GapWidth=0.3, # default gap of 0.3 mm
                   STLMeshSize=3.0,
                   Debug=False,
                   license_key=""):
        """Initialize the DelamoModeler, creating the codestores for initinstrs, assemblyinstrs, bcinstrs, meshinstrs, and runinstrs, and
        specifying the default tolerances"""
        
        # Create code stores
        initinstrs=codestore(globals=globals)

        assemblyinstrs=codestore(initinstrs.codevariables) # share same variable dictionary as init instructions
        bcinstrs=codestore(initinstrs.codevariables) # share same variable dictionary as init instructions
        meshinstrs=codestore(initinstrs.codevariables) # share same variable dictionary as init instructions
        fiberinstrs = codestore(initinstrs.codevariables)  # share same variable dictionary as init instructions
        runinstrs=codestore(initinstrs.codevariables) # share same variable dictionary as init instructions

        # imports from front matter in abqtemplate.py
        _abq=initinstrs.preexisting_variable("abq")        
        globals["abq"]=_abq

        _abq_assembly=assemblyinstrs.rewrapobj(_abq)
        globals["abq_assembly"]=_abq_assembly

        # abqC needs to be accessible when functions/classes
        # are defined because it is used for default
        # parameters, hence the namebinding_wrapper above
        _abqC=initinstrs.preexisting_variable("abqC")
        initinstrs.add_named_binding(abqC_wrapper,_abqC)
        globals["abqC"]=_abqC
        
        _regionToolset=assemblyinstrs.preexisting_variable("regionToolset") # Region operations wrapped into the assemblyinstrs pool
        globals["regionToolset"]=_regionToolset
        
        _mesh=meshinstrs.preexisting_variable("mesh")  # mesh operations wrapped into the meshinstrs pool
        globals["mesh"]=_mesh

        _section=initinstrs.preexisting_variable("section")  # section operations wrapped into the initinstrs pool
        globals["section"]=_section

        _assembly=assemblyinstrs.preexisting_variable("assembly")  # assembly operations wrapped into the assemblyinstrs pool
        globals["assembly"]=_assembly

        _connector=bcinstrs.preexisting_variable("connector")  # connector operations wrapped into the bcinstrs pool
        globals["connector"]=_connector


        
        # Wrap classes that will be used for code generation
        # Instructions go into the different codestores according to
        # which codestore was used to do the wrapping
        InitDict=initinstrs.wrap_class(dict)   # Dict that will operate during the initialization step

        autofiber=fiberinstrs.preexisting_variable("AutoFiber_abq")

        FEModel=assemblyinstrs.preexisting_variable("FEModel")
        LaminateAssembly=assemblyinstrs.preexisting_variable("LaminateAssembly")
        stepfile=assemblyinstrs.preexisting_variable("stepfile")
        
        # Body number database:
        # Two parts, first a dictionary filled out during initialization phase
        BodyNumDB_build = InitDict() # Define body number database as a dictionary for which operations will execute at initialization phase 
        
        # Second, BodyNumDB, a dictionary accessed later
        # Dictionary by bodyname of bodynum
        BodyNumDB=assemblyinstrs.rewrapobj(BodyNumDB_build)
        
        # Layer body database
        # Dictionary by layername of list of body names
        #BodyDB_build = InitDict() # define layer body database as a dictionary which will execute in initialization phase
        
        #BodyDB=assemblyinstrs.rewrapobj(BodyDB_build)
            
        # Initialize the OpenCascade Model Builder
        modelbuilder = OCCModelBuilder(PointTolerance=pointtolerance,NormalTolerance=normaltolerance,GapWidth=GapWidth,STLMeshSize=STLMeshSize,Debug=Debug)

        
        # Define default finite element tolerances
        # facepointtolerance = assemblyinstrs.assign_variable("facepointtolerance", facepointtolerancefactor*modelbuilder.tolerance())  # point positioning tolerance, in abaqus units (mm)g
        # normaltolerance = assemblyinstrs.assign_variable("normaltolerance", normaltolerance)
        abqpointtolerance=pointtolerancefactor*modelbuilder.PointTolerance
            
        DM=cls(initinstrs=initinstrs,
               assemblyinstrs=assemblyinstrs,
               bcinstrs=bcinstrs,
               meshinstrs=meshinstrs,
               fiberinstrs=fiberinstrs,
               runinstrs=runinstrs,
               globals=globals,
               abq=_abq,
               abq_assembly=_abq_assembly,
               abqC=_abqC,
               regionToolset=_regionToolset,
               mesh=_mesh,
               section=_section,
               assembly=_assembly,
               connector=_connector,
               InitDict=InitDict,
               autofiber=autofiber,
               FEModel=FEModel,
               LaminateAssembly=LaminateAssembly,
               stepfile=stepfile,
               BodyNumDB_build=BodyNumDB_build,
               BodyNumDB=BodyNumDB,
               #BodyDB_build=BodyDB_build,
               #BodyDB=BodyDB,
               modelbuilder=modelbuilder,
               abqpointtolerance=abqpointtolerance,
               normaltolerance=normaltolerance,
               tangenttolerance=tangenttolerance)
        
        # Insert abqfuncs libraries ...

        #   ... find installed directory
        thisfile=sys.modules[cls.__module__].__file__
        thisdir=os.path.split(thisfile)[0]

        # ...abqfuncs_assembly.py library
        DM.abaqus_assembly_functions(os.path.join(thisdir,"abqfuncs_assembly.py"),globals)

        # ...abqfuncs_bc.py library
        DM.abaqus_bc_functions(os.path.join(thisdir,"abqfuncs_bc.py"),globals)

        # ...abqfuncs_mesh.py library
        DM.abaqus_mesh_functions(os.path.join(thisdir,"abqfuncs_mesh.py"),globals)

        # ...abqfuncs_fiber.py library
        DM.abaqus_fiber_functions(os.path.join(thisdir,"abqfuncs_fiber.py"),globals)

        # ...abqfuncs_run.py library
        DM.abaqus_run_functions(os.path.join(thisdir,"abqfuncs_run.py"),globals)


        return DM

    #def _AppendAssemblyComponents(to_be_saved,assembly):
    #    for partname in assembly.parts:
    #        to_be_saved.append(assembly.parts[partname])
    #        pass
    #    for assemblyname in assembly.assemblies:
    #        self._AppendAssemblyComponents(to_be_saved,assembly.assemblies[assemblyname])
    #        pass
    #    pass
        
    def Finalize(self,script_to_generate,cad_file_path_from_script):
        """Write out the generated ABAQUS script and generated CAD file given the desired output filenames."""
        # target_cad_file_path, if given, will be joined with cad_file_name in the generated code. 
        # If target_cad_file_path is None, will encode full absolute path of cad_file_name into generated code  
        # Save SAT file using the std::vector implementation

        #to_be_saved=delamo.CADwrap.MBBodyList()
        #for part in parts:
        #    to_be_saved.append(part.gk_part)
        #    pass
        #for assembly in assemblies:
        #    self._AppendAssemblyComponents(to_be_saved,assembly)
        #    pass
        
        
        BodyNameList=[] 

        script_directory=os.path.split(script_to_generate)[0]
        cad_file_name=os.path.join(script_directory,cad_file_path_from_script)
        
        BodyNameList = self.modelbuilder.save(cad_file_name, self.to_be_saved)
        #self.modelbuilder.save(cad_file_name, layer_list, BodyNameList)
        # BodyNameList now contains a list of body names, ordered the same as the body numbers in the .sat file

        # Build BodyNumDB  ... this runs in immediate context
        # so that BodyNumDB is available during assembly
        for cnt in range(len(BodyNameList)):
            self.BodyNumDB_build[BodyNameList[cnt]]=cnt+1
            pass

        # Build BodyDB
        #for Body in self.to_be_saved:
        ##for layer in layer_list:
        #    self.BodyDB_build[Body.name()]=layer.bodynames(0)
        #    pass
            
            
        # assign .sat file name in initinstrs
        #if target_cad_file_path is None:
        #    acisfile=self.initinstrs.assign_variable("acisfile", os.path.abspath(cad_file_name))  # Primary ACIS geometry file to load
        #    pass
        #else:
        #    target_os_path_join = self.initinstrs.wrap_function(os.path.join)
        #    acisfile=self.initstrs.assign_variable("acisfile",target_os_path_join(target_cad_file_path,cad_file_name))
        #    pass
        stepfile=self.initinstrs.assign_variable("stepfile", cad_file_path_from_script)
        
        # Write out the generated python code
        # -----------------------------------------------

        write_abq_script(self.initinstrs,self.assemblyinstrs,self.bcinstrs,self.meshinstrs,self.fiberinstrs,self.runinstrs,script_to_generate)


        pass

    pass


class Part(object):
    """The Part class represents a CAD part and a parallel finite element (ABAQUS) part."""

    DM = None
    """DelamoModeler object"""
    
    name = None
    """Name of this part"""
    
    fe_part = None
    """Wrapper for ABAQUS part object, assembly phase"""
    
    fe_part_meshing = None
    """Wrapper for Abaqus part object, meshing phase"""
    
    fe_inst = None
    """Wrapper for Abaqus part instance, assembly phase"""
    
    shell = None
    """True/False: Is a shell (as opposed to a solid) part"""


    
    fe_datum_csys=None
    """Abaqus datum CSYS"""

    fe_materialorientation=None
    """Abaqus material orientation"""

    gk_layerbody=None
    """Geometry kernel representation of this part"""
    
    def __init__(self,**kwargs):
        self.shell=False
        for key in kwargs:
            assert(hasattr(self,key))
            setattr(self,key,kwargs[key])
            pass

        if self.fe_part_meshing is None and self.fe_part is not None:
            self.fe_part_meshing=self.DM.meshinstrs.rewrapobj(self.fe_part)

            pass
        
        pass


    
    def CreateInstance(self,dependent=abqC.ON):
        """ Create the finite element instance of this part """
        self.fe_inst=self.DM.LaminateAssembly.Instance(
            name='Instance-%s' % (self.name),
            part=self.fe_part,
            dependent=dependent)
        return self.fe_inst

    @classmethod
    def FromGeometryFile(cls,DM,filename,name,bodynum):
        """Define a finite element solid or shell part from a file on disk. This will not 
        have a geometry kernel representation"""
        # Always reopen file in same set of instructions (assemblyinstrs in this case)
        # immediately before calling PartFromGeometryFile() because
        # Abaqus mixes up files if you have more than one open. 
        geom_fh=DM.abq_assembly.mdb.openStep(filename,scaleFromFile=abqC.ON)
        fe_part=DM.FEModel.PartFromGeometryFile(
            geometryFile=geom_fh,
            bodyNum=bodynum, #M.globals["LookupBodies"](M.BodyNumDB,bodyname),
            combine=False,
            dimensionality=abqC.THREE_D,
            type=abqC.DEFORMABLE_BODY)
        
        newpart = cls(DM=DM,
                      name=name,fe_part=fe_part,gk_layerbody=gk_layerbody)

        DM.to_be_saved.append(gk_layerbody)
        # Create instance of new part
        newpart.CreateInstance(dependent=abqC.ON)
        
        return newpart
    

    @classmethod
    def FromGK3D(cls,DM,gk_layerbody,shell=False,no_FE_instance=False,omit_from_FE=False):
        """Define a solid or shell part from a 3D geometry kernel object. Ensures that 
        this geometry kernel object will be saved to the generated CAD file
        and programs ABAQUS to load it in and create an instance. """
        
        # Always reopen file in same set of instructions (assemblyinstrs in this case)
        # immediately before calling PartFromGeometryFile() because
        # Abaqus mixes up files if you have more than one open. 
        if not omit_from_FE:
            geom_fh=DM.abq_assembly.mdb.openStep(DM.stepfile,scaleFromFile=abqC.ON)
            fe_part=DM.FEModel.PartFromGeometryFile(
                name=gk_layerbody.Name,
                geometryFile=geom_fh,
                bodyNum=DM.BodyNumDB[gk_layerbody.Name], #M.globals["LookupBodies"](M.BodyNumDB,bodyname),
                combine=False,
                dimensionality=abqC.THREE_D,
                type=abqC.DEFORMABLE_BODY)
            pass
        else:
            fe_part=None
            pass

        newpart = cls(DM=DM,
                      name=gk_layerbody.Name,fe_part=fe_part,gk_layerbody=gk_layerbody,shell=shell)

        if not omit_from_FE:
            DM.to_be_saved.append(gk_layerbody)

            # Create instance of new part
            if not no_FE_instance:
                newpart.CreateInstance(dependent=abqC.ON)
                pass
            pass

        return newpart



    def AssignSection(self,Section,offset=0.0,offsetType=abqC.MIDDLE_SURFACE):
        """Assign ABAQUS finite element Section (which must be an Abaqus section object with a material applied to it) to this part. 
Calls SectionAssigment() method on the part"""
        
        # region is all cells or all faces from Lamina
        if self.shell:
            region=self.DM.regionToolset.Region(faces=self.fe_part.faces)
            pass
        else:            
            region=self.DM.regionToolset.Region(cells=self.fe_part.cells)
            pass

        self.fe_part.SectionAssignment(region=region,
                                       sectionName=Section.name,
                                       offset=offset,
                                       offsetType=offsetType,
                                       offsetField='',
                                       thicknessAssignment=abqC.FROM_SECTION)
        
        pass


    # def ReferenceSurfaceFromPart(self,surfname,surface_points_and_normals,facepointtolerance,normaltolerance):
    #    # Create surface associated with the part. Supply parameters to locate face(s)
    #    abqsurf=self.fe_part.Surface(name=surfname, side1Faces=self.DM.globals["GetFaces"](self.fe_part.faces,surface_points_and_normals,facepointtolerance,normaltolerance))
    #    return abqsurf
    


    
    def GetPartFace_point_normal(self,surface_point_and_normal,pointtolerance,normaltolerance):
        """Method to find a face of the part given a point and normal and tolerances"""
        return self.DM.globals["GetFace_point_normal"](self.fe_part.faces,
                                          surface_point_and_normal,
                                          pointtolerance,
                                          normaltolerance)

    def GetPartFace(self,surface_points,pointtolerance):
        """Method to find a face of the part given a point and tolerance"""
        return self.DM.globals["GetFace"](self.fe_part.faces,
                                          surface_points,
                                          pointtolerance)

    def GetMultiplePartFaces(self,surface_points,pointtolerance):
        """Method to find multiple faces of the part given a point and tolerance"""
        return self.DM.globals["GetMultipleFaces"](self.fe_part.faces,
                                                   surface_points,
                                                   pointtolerance)
    
    def GetPartEdge_point_tangent(self,edgepointtangent,pointtolerance,tangenttolerance):
        """Method to find an edge of the part given a point and tangent and tolerances"""
        return self.DM.globals["GetEdge_point_tangent"](self.fe_part.edges, self.fe_part.vertices, edgepointtangent,pointtolerance,tangenttolerance)

    def GetPartEdge(self,edgepoints,pointtolerance):
        """Method to find an edge of the part given a point and tolerance"""
        return self.DM.globals["GetEdge"](self.fe_part.edges, self.fe_part.vertices,edgepointtangent,pointtolerance)

    def GetMultiplePartEdges(self,edgepoints,pointtolerance):
        """Method to find multiple edges of the part given a point and tolerance"""
        return self.DM.globals["GetMultipleEdges"](self.fe_part.edges, self.fe_part.vertices,edgepoints,pointtolerance)

    
    #def GetPartNodes(self,nodepoints,pointtolerance):
    #    raise NotImplementedError('Use the GetPartVertices() method to find points associated with the part geometry. Use the GetInstanceNodes() method to find points associated with the instance mesh')
    #    
    #def GetPartVertices(self,nodepoints,pointtolerance):
    #    # Method to find vertices associated with the part
    #    # Vertices are associated with geometric features of the part
    #    # Nodes are associated with a mesh instance
    #    return GetVertices(self.fe_part.vertices,nodepoints,pointtolerance)
    
    
    #def GetPartEdges_ThreePoints(self,edgeandinteriorpoints,pointtolerance):
    #    # Method to find edges associated with the part
    #    return self.DM.globals["GetEdges_ThreePoints"](self.fe_part.edges, self.fe_part.vertices, edgeandinteriorpoints,pointtolerance)

    
    # Instance methods identical to part methods (do we need these as well?)
    def GetInstanceFaceRegion_point_normal(self,surface_point_and_normal,pointtolerance,normaltolerance):
        """Get an instance face region based on a point and normal and tolerances. Note that this creates a "Set-like" region, not the "Surface-like" regions used in the bonding functions(below) in api.py. See Abaqus Scripting Reference Guide section 45.3 for the distinction"""
        return self.DM.regionToolset.Region(faces=self.GetInstanceFace_point_normal(surface_point_and_normal,pointtolerance,normaltolerance))
    
    def GetInstanceFace_point_normal(self,surface_point_and_normal,pointtolerance,normaltolerance):
        """Get an instance face based on a point and normal and tolerances. """
        return self.DM.globals["GetFace_point_normal"](self.fe_inst.faces,surface_point_and_normal,pointtolerance,normaltolerance)

    def GetInstanceFaceRegion(self,surface_points,pointtolerance):
        """Get an instance face region based on a point and tolerance. Note that this creates a "Set-like" region, not the "Surface-like" regions used in the bonding functions(below) in api.py. See Abaqus Scripting Reference Guide section 45.3 for the distinction"""
        return self.DM.regionToolset.Region(faces=self.GetInstanceFace(surface_points,pointtolerance))

    def GetInstanceFaceRegionSurface(self,surface_points,pointtolerance):
        """Get an instance face region based on a point and tolerance. Note that this creates a "Surface-like" region used in the bonding functions(below) in api.py, not the "Set-like" regions used in some other contexts. See Abaqus Scripting Reference Guide section 45.3 for the distinction."""
        return self.DM.regionToolset.Region(side1Faces=self.GetInstanceFace(surface_points,pointtolerance))

    
    def GetInstanceFace(self,surface_points,pointtolerance):
        """Get an instance face based on a point and tolerance. """
        return self.DM.globals["GetFace"](self.fe_inst.faces,surface_points,pointtolerance)

    def GetMultipleInstanceFacesRegion(self,surface_points,pointtolerance):
        """Get an instance face region based on multile faces selected by surface points and tolerance. Note that this creates a "Set-like" region, not the "Surface-like" regions used in the bonding functions(below) in api.py. See Abaqus Scripting Reference Guide section 45.3 for the distinction"""
        return self.DM.regionToolset.Region(faces=self.GetMultipleInstanceFace(surface_points,pointtolerance))

    def GetMultipleInstanceFaces(self,surface_points,pointtolerance):
        """Get instance faces based on points and a tolerance. """
        return self.DM.globals["GetMultipleFace"](self.fe_inst.faces,surface_points,pointtolerance)
    
    def GetInstanceEdgeRegion_point_tangent(self,edgepointtangent,pointtolerance,tangenttolerance):
        """Get an instance edge region based on a point and tangent and tolerances. Note that this creates a "Set-like" region, not the "Surface-like" regions used in the bonding functions(below) in api.py. See Abaqus Scripting Reference Guide section 45.3 for the distinction"""
        return self.DM.regionToolset.Region(edges=self.GetInstanceEdge_point_tangent(edgepointtangent,pointtolerance,tangenttolerance))
    def GetInstanceEdge_point_tangent(self,edgepointtangents,pointtolerance,tangenttolerance):
        """Get an instance edge based on a point and tangent and tolerances. """
        return self.DM.globals["GetEdge_point_tangent"](self.fe_inst.edges, self.fe_inst.vertices,edgepointtangent,pointtolerance,tangenttolerance)

    def GetInstanceEdgeRegion(self,edgepoints,pointtolerance):
        """Get an instance edge region based on a point and tolerance. Note that this creates a "Set-like" region, not the "Surface-like" regions used in the bonding functions(below) in api.py. See Abaqus Scripting Reference Guide section 45.3 for the distinction"""
        return self.DM.regionToolset.Region(edges=self.GetInstanceEdge(edgepoints,pointtolerance))
    def GetInstanceEdgeRegionSurface(self,edgepoints,pointtolerance):
        """Get an instance edge region based on a point and tolerance. Note that this creates a "Surface-like" region with side1Edges= used in the bonding functions(below) in api.py, not the "Set-like" regions used in some other contexts. See Abaqus Scripting Reference Guide section 45.3 for the distinction."""
        # Note that this creates a "Surface-like" edge region using "side1Edges"
        return self.DM.regionToolset.Region(side1Edges=self.GetInstanceEdge(edgepoints,pointtolerance))
    def GetInstanceEdge(self,edgepoints,pointtolerance):
        """Get an instance edge based on a point and tolerance. """
        return self.DM.globals["GetEdge"](self.fe_inst.edges,self.fe_inst.vertices,edgepoints,pointtolerance)

    def GetMultipleInstanceEdgesRegion(self,edgepoints,pointtolerance):
        """Get instance edge regions based on a point and tolerance. Note that this creates a "Set-like" region, not the "Surface-like" regions used in the bonding functions(below) in api.py. See Abaqus Scripting Reference Guide section 45.3 for the distinction"""
        return self.DM.regionToolset.Region(edges=self.GetMultipleInstanceEdges(edgepoints,pointtolerance))
    def GetMultipleInstanceEdges(self,edgepoints,pointtolerance):
        """Get multiple instance edges based on points and tolerance. """
        return self.DM.globals["GetMultipleEdges"](self.fe_inst.edges,self.fe_inst.vertices,edgepoints,pointtolerance)

    #def GetInstanceNodes(self,nodepoints,pointtolerance):
    #    return GetNodes(self.fe_inst.nodes,nodepoints,pointtolerance)


    #def GetInstanceEdges_ThreePoints(self,edgeandinteriorpoints,pointtolerance):
    #    return self.DM.globals["GetEdges_ThreePoints"](self.fe_inst.edges, self.fe_inst.vertices,edgeandinteriorpoints,pointtolerance)

    def SeedPartEdgesByFaces(self,surface_points_and_normals,pointtolerance,normaltolerance,meshsize):
        """Seed the edges around a set of faces with a particular meshing size. This is used for localized mesh refinement. 
The actual implementation is the ABAQUS code in abqfuncs_mesh.py"""

        self.DM.globals["SeedPartEdgesByFaces"](self.fe_part,surface_points_and_normals,pointtolerance,normaltolerance,meshsize)
        pass
    
    def MeshSimple(self,ElemTypes,meshsize,ElemShape=None,ElemTechnique=None,refined_edges=[],pointtolerance=None,tangenttolerance=None,refinedmeshsize=None,DeviationFactor=None,MinSizeFactor=None):
        """Perform meshing of this Part. meshsize is nominal meshing size
        ElemShape should be abqC.HEX, abqC.HEX_DOMINATED (default), abqC.TET, etc. for solids or abqC.QUAD or abqC.QUAD_DOMINATED (default), or abqC.TRI for shells
        ElemTechnique can be abqC.SYSTEM_ASSIGN (default), abqC.FREE, or abqC.STRUCTURED
        refined edges is  a list of tuples: ((point1, tangent1),(point2,tangent2)) on each edge to be specially refined
        pointtolerance is tolerance size for finding refined edges
        refinedmeshsize is size of refined mesh regions
        DeviationFactor and MinSizeFactor are scaling factors for seeding (default 0.1) 
        """

        if self.shell:
            if ElemShape is None:
                ElemShape=abqC.QUAD_DOMINATED
                pass
            pass
        else:
            if ElemShape is None:
                ElemShape=abqC.HEX_DOMINATED
                pass
            pass

        if ElemTechnique is None:
            ElemTechnique=abqC.SYSTEM_ASSIGN
            pass

        if DeviationFactor is None:
            DeviationFactor=0.1
            pass
        
        if MinSizeFactor is None:
            MinSizeFactor=0.1
            pass

        if self.shell:
            self.fe_part_meshing.setElementType(regions=(self.fe_part_meshing.faces,),elemTypes=ElemTypes)

            
            self.fe_part_meshing.setMeshControls(regions=self.fe_part_meshing.faces,elemShape=ElemShape,technique=ElemTechnique)
            pass
        else:
            self.fe_part_meshing.setElementType(regions=(self.fe_part_meshing.cells,),elemTypes=ElemTypes)
            
        
            self.fe_part_meshing.setMeshControls(regions=self.fe_part_meshing.cells,elemShape=ElemShape,technique=ElemTechnique)
            pass

    
        self.fe_part_meshing.seedPart(size=meshsize,deviationFactor=DeviationFactor,minSizeFactor=MinSizeFactor)
        
        # refined_edges is a list of tuples ... endpoint1, interiorpoint, endpoint2
        # representing 3 points on the edge, with two of them being endpoints. 

        if len(refined_edges) > 0: 
            # edges defined for refinement
            #picked=self.GetPartEdges_ThreePoints(refined_edges,pointtolerance)
            picked=self.GetMultiplePartEdges(refined_edges,pointtolerance)
            self.fe_part_meshing.seedEdgeBySize(edges=picked,size=refinedmeshsize,deviationFactor=0.1,minSizeFactor=0.1,constraint=abqC.FINER)
            pass
        
        
        self.fe_part_meshing.generateMesh()
        pass
    



    def MeshCohesive(self,meshsize,ElemShape=None,Algorithm=None,ElemLibrary=None,refined_edges=[],pointtolerance=None,refinedmeshsize=None,DeviationFactor=None,MinSizeFactor=None,SweepSense=None):
        """Meshing routine for cohesive layers: 
        meshsize is nominal meshing size
        ElemShape should be abqC.HEX or abqC.HEX_DOMINATED (default)
        Algorithm can be abqC.ADVANCING_FRONT (default) or abqC.MEDIAL_AXIS
        ElemLibrary can be abqC.STANDARD (default) or abqC.EXPLICIT
        refined edges is a list of tuples: (endpoint1, interiorpoint, endpoint2) representing 3 points on each edge to be specially refined
        pointtolerance is tolerance size for finding refined edges
        refinedmeshsize is size of refined mesh regions
        DeviationFactor and MinSizeFactor are scaling factors for seeding (default 0.1) 
        SweepSense should be abqC.FORWARD (default) or abqC.REVERSE """
        
        assert(not self.shell)
        if ElemShape is None:
            ElemShape=abqC.HEX_DOMINATED
            pass
                
        if Algorithm is None:
            Algorithm=abqC.ADVANCING_FRONT
            pass
            
        if ElemLibrary is None:
            ElemLibrary=abqC.STANDARD
            pass
            
        if DeviationFactor is None:
            DeviationFactor=0.1
            pass
        
        if MinSizeFactor is None:
            MinSizeFactor=0.1
            pass
            
        if SweepSense is None:
            SweepSense=abqC.FORWARD
            pass

        cells=self.fe_part_meshing.cells # .getSequenceFromMask(...)

        #sys.stderr.write("WARNING: Assuming reference region for stacking direction is faces[0] (FIXME)\n")
        #refregion=self.fe_part_meshing.faces[0]

        # reference region should be the largest face
        refregion=self.DM.globals["GetLargestFace"](self.fe_part_meshing,self.fe_part_meshing.faces)
        
        self.fe_part_meshing.assignStackDirection(referenceRegion=refregion,cells=cells)
        
        self.fe_part_meshing.setMeshControls(regions=cells,
                                             elemShape=ElemShape,
                                             technique=abqC.SWEEP,
                                             algorithm=Algorithm)
        sys.stderr.write("WARNING: Assuming region for sweep path is cells[0] (FIXME) and edge is edges[1]\n")
        # !!!*** ALSO need to add contact model for delaminated region across cohesive layer ***!!!
        
        assembly_assert=self.DM.assemblyinstrs.preexisting_variable("assert")
        assembly_len=self.DM.assemblyinstrs.preexisting_variable("len")

        assembly_assert(assembly_len(self.fe_part_meshing.cells)==1)   # Given how we construct the cohesive layer it should never have anything but exactly one cell
        #self.fe_part_meshing.NEED_TO_DETERMINE_SWEEP_PATH

        (SweepPathEdgePoint, SweepPathEdgeTangent) = self.gk_layerbody.GetOffsetEdge()
        # Get proxied Abaqus edge object. 
        SweepPathEdge=self.GetPartEdge_point_tangent((SweepPathEdgePoint,SweepPathEdgeTangent),self.DM.abqpointtolerance,self.DM.tangenttolerance)

        self.fe_part_meshing.setSweepPath(region=self.fe_part_meshing.cells[0],edge=SweepPathEdge,sense=SweepSense)
        
        HexElemType = self.DM.mesh.ElemType(elemCode=abqC.COH3D8,elemLibrary=ElemLibrary)
        WedgeElemType = self.DM.mesh.ElemType(elemCode=abqC.COH3D6,elemLibrary=ElemLibrary)
        
        self.fe_part_meshing.setElementType(regions=(cells,),elemTypes=(HexElemType,WedgeElemType))
        self.fe_part_meshing.seedPart(size=meshsize,deviationFactor=DeviationFactor,minSizeFactor=MinSizeFactor)
        
        self.fe_part_meshing.generateMesh()
        pass
    
    pass

class LayerPart(Part):
    """A LayerPart is a Part that will be used in a Layer"""
    
    # So far no extra member variables
    
    def __init__(self,**kwargs):
        self.shell=False
        for key in kwargs:
            assert(hasattr(self,key))
            setattr(self,key,kwargs[key])
            pass
        
        if self.fe_part_meshing is None and self.fe_part is not None:
            self.fe_part_meshing=self.DM.meshinstrs.rewrapobj(self.fe_part)

            pass
        pass

    def ApplyLayup(self,coordsys,layupdirection,fiberorientation):
        """Assign self.fe_datum_csys and self.fe_materialorientation"""

        if fiberorientation is not None:
            self.fe_materialorientation = self.DM.fiberinstrs.rewrapobj(self.fe_part.MaterialOrientation)
            self.fe_materialorientation = self.fe_materialorientation(
                region=self.DM.regionToolset.Region(
                    cells=self.fe_part.cells),
                orientationType=abqC.FIELD, axis=abqC.AXIS_3,
                fieldName='%s_orientation' % self.name,
                localCsys=None,
                additionalRotationType=abqC.ROTATION_NONE,
                angle=0.0,
                additionalRotationField='',
                stackDirection=abqC.STACK_3)
        else:
            coordsys.ApplyLayup(self, layupdirection)
        pass


    pass
    

class Assembly(object):
    """The assembly class represents an assembly of parts and assemblies."""
    
    name=None
    """ The name of the assembly"""

    parts=None
    """An ordered dictionary by name of Part objects or Part subclasses"""

    
    assemblies=None
    """An ordered dictionary by name of Assembly objects or Assembly subclasses"""

    # There is also a synthetic attribute,
    # partlist = None  ... which is a list of the parts

    # and a synthetic attribute
    # singlepart = None  which only works if this is a one-part assembly
    def __init__(self,**kwargs):
        self.parts=collections.OrderedDict()
        self.assemblies=collections.OrderedDict()
        for key in kwargs:
            assert(hasattr(self,key))
            setattr(self,key,kwargs[key])
            pass
        
        pass

    @property
    def partlist(self):
        """This property gives a list of the parts in this assembly and 
    all subassemblies """
        parts = [ self.parts[partname] for partname in self.parts ]

        # also include anything within an assembly
        for assemblyname in self.assemblies:
            parts.extend(self.assemblies[assemblyname].partlist)
            pass
        return parts
    

    @property
    def singlepart(self):
        """Assuming this assembly contains only a single part, this property gives that part. Otherwise raises IndexError()."""
        partkeys=list(self.parts.keys())
        assemblykeys=list(self.assemblies.keys())
        if len(partkeys)!=1 or len(assemblykeys) != 0:
            raise IndexError("onlypart attribute is only valid for assemblies that contain a single part")
        return self.parts[partkeys[0]]
        
    def MeshSimple(self,*args,**kwargs):
        """This method iterates over the MeshSimple methods of each 
        part that is a child of this assembly (but not to subassemblies)"""
        # NOTE: Only applies to direct part children, not subassemblies
        for name in self.parts:
            self.parts[name].MeshSimple(*args,**kwargs)
            pass
        pass

    def MeshCohesive(self,*args,**kwargs):
        """This method iterates over the MeshCohesive methods of each 
        part that is a child of this assembly (but not to subassemblies)"""
        # NOTE: Only applies to direct part children, not subassemblies
        for name in self.parts:
            self.parts[name].MeshCohesive(*args,**kwargs)
            pass
        pass

    @classmethod
    def FromParts(cls,name,*args):
        """Creates an assembly given the name parameter, and 
        additional parameters representing the parts to include
        in the assembly"""
        assem=cls(name=name)
        
        for part in args:
            assem.parts[part.name]=part
            pass

        return assem
    @classmethod
    def FromAssemblies(cls,name,*args):
        """Creates an assembly given the name parameter, and 
        additional parameters representing the subassemblies to include
        in the assembly"""

        assem=cls(name=name)
        
        for assembly in args:
            assem.assemblies[assembly.name]=assembly
            pass

        return assem

    @classmethod
    def FromPartsAndAssemblies(cls,name,parts,assemblies):
        """Creates an assembly given the name parameter, and 
        lists of parts and subassemblies to include
        in the assembly"""

        assem=cls(name=name)
        
        for part in parts:
            assem.parts[part.name]=part
            pass
        for assembly in assemblies:
            assem.assemblies[assembly.name]=assembly
            pass
        
        return assem
    pass

class CoordSys(object):
    """ Abstract class representing a coordinate system. 
Each CoordSys subclass should implement the method ApplyLayup(self,layerpart,layerdirection)"""
    pass 

class SimpleCoordSys(CoordSys):
    """Concrete implementation of a CoordSys representing a fixed Cartesian coordinate frame"""
    
    # fibervec and crossfibervec correspond to a 0 deg ply
    fibervec=None
    """Unit vector along the fibers of a 0 degree ply"""
    
    crossfibervec=None
    """Unit vector along the fibers of a 90 degree ply"""

    outofplanevec=None
    """Out-of-plane unit vector, i.e. fibervec cross crossfibervec"""
    
    def __init__(self,fibervec,crossfibervec):
        self.fibervec=fibervec/np.linalg.norm(fibervec)
        # correct crossfibervec if necessary by Gram-Schmidt orthonormalization
        self.crossfibervec = crossfibervec - np.inner(self.fibervec,crossfibervec)*self.fibervec
        self.outofplanevec=np.cross(fibervec,crossfibervec)

        # multiply oriented_to_xyz_mat on right by coordinates in oriented frame to get coordinates in (x,y,z) space
        self.oriented_to_xyz_mat = np.array((self.fibervec,self.crossfibervec,self.outofplanevec),dtype='d').T
        # Multiply xyz_to_oriented_mat on right by (x,y,z) coordinates to get coordinates in oriented frame. 
        self.xyz_to_oriented_mat = self.oriented_to_xyz_mat.T 
        pass

    def ApplyLayup(self,layerpart,layupdirection):
        """Set the material orientation of the specified layerpart to layupdirection (in degrees) relative to this coordinate frame"""

        oriented_point1 = (np.cos(layupdirection*np.pi/180.0),
                           np.sin(layupdirection*np.pi/180.0),
                           0.0)
        oriented_point2 = (-np.sin(layupdirection*np.pi/180.0),
                           np.cos(layupdirection*np.pi/180.0),
                           0.0)

        xyz_point1=np.dot(self.oriented_to_xyz_mat,oriented_point1)
        xyz_point2=np.dot(self.oriented_to_xyz_mat,oriented_point2)
        layerpart.fe_datum_csys=layerpart.fe_part.DatumCsysByThreePoints(coordSysType=abqC.CARTESIAN,
                                                                         origin=(0.0,0.0,0.0),
                                                                         point1=tuple(xyz_point1),
                                                                         point2=tuple(xyz_point2),
                                                                         name="%s_csys" % (layerpart.name))
        
        
        # Can we simplify things by not defining a new CSYS for each, but just applying
        # angle, below??? ... would use DatumCsysByDefault() of Feature object
        
        layerpart.fe_materialorientation=layerpart.fe_part.MaterialOrientation(additionalRotationField='',
                                                                     additionalRotationType=abqC.ROTATION_NONE,
                                                                     angle=0.0,
                                                                     axis=abqC.AXIS_3,
                                                                     fieldName='',
                                                                     localCsys=layerpart.fe_part.datums[layerpart.fe_datum_csys.id],
                                                                     orientationType=abqC.SYSTEM,
                                                                     region=layerpart.DM.regionToolset.Region(cells=layerpart.fe_part.cells),
                                                                     stackDirection=abqC.STACK_3)
        
        pass
    
    pass
    

class Layer(Assembly):
    """Represents a layer of a composite material"""
    
    gk_layer=None
    """Underlying geometry kernel object"""
    
    layupdirection=None
    """Orientation of this layer, in degrees"""
    
    # fe_parts_builder=None
    LayerSection=None
    """ABAQUS Section used by this layer"""
    
    coordsys=None
    """Coordinate system used by this layer"""
    
    def __init__(self,**kwargs):
        self.parts=None
        self.assemblies=collections.OrderedDict()
        self.orientation=None
        for key in kwargs:
            assert(hasattr(self,key))
            setattr(self,key,kwargs[key])
            pass
        
        pass

    def Split(self,splitpath_filename,PointTolerance):
        self.gk_layer.Split(splitpath_filename,PointTolerance)
        pass

    def GetPartInstanceFaceRegionFromPoint(self,Point,PointTolerance):
        # May only be called after finalization
        gk_layerbodyname = self.gk_layer.FindLayerBodyNameByPoint(Point,PointTolerance)
        part = self.parts[gk_layerbodyname]
        region=part.GetInstanceFaceRegion(Point,PointTolerance) # get ABAQUS region
        return region
    
    @classmethod
    def CreateFromParams(cls,DM,create_params,name,LayerSection,layupdirection, meshsize=0.5, split=None, coordsys=None):
        """Create a layer given creation parameters to be passed to the geometry kernel, a name, ABAQUS Section, 
layup direction, etc."""
        # create_params should be a tuple with all parameters to create_layer, except for the gk_layer object at the end
        # e.g. mold, thickness, or layer, direction, thickness

        # ***!!!! Should this still be a layer object or a more generic part object???
        # ***!!!!!  OBSOLETE... NEEDS AN UPDATE 
        gk_layer=OCCLayer()

        actual_create_params=list(create_params)
        actual_create_params.append(gk_layer)
        # print(str(actual_create_params))
        DM.modelbuilder.create_layer(*actual_create_params)
        gk_layer.name(name)
        gk_layer.layup(layupdirection)  # is this really necessary???
        if split is not None:  # do layer splitting here, there seems no other option to change layer after creating it on the python side
            DM.modelbuilder.split_layer(gk_layer, split)
		
        # NOTE: Does not fill out "parts" member... that waits
        # until "Finalize" method when the part list is extracted from
        # the geometry kernel

        # Is this an appropriate place for this?
        # gk_layer.DMObj = gk_layer.CreateDMObject(gk_layer.RefMold.Shape, meshsize)

        return cls(name=name,gk_layer=gk_layer,layupdirection=layupdirection,LayerSection=LayerSection,coordsys=coordsys)

    def CreateFiberObject(self, DM, point, fibervec, normal, mp, fiberint=1.0, angle_error=0.01, final_plotting=False):
        """ Utilize Autofiber orientation package to calculate optimal fiber orientations at each mesh
         element centroid. """
        self.orientation = AutoFiber(self.gk_layer.DMObj,
                                     point, fibervec, normal,
                                     materialproperties=(
                                         [mp[0], mp[1], mp[2]],
                                         [mp[3], mp[4], mp[5]],
                                         [mp[6], mp[7], mp[8]]),
                                     fiberint=fiberint,
                                     angle_error=angle_error)

        self.LayupFiberObject(DM, self.layupdirection, final_plotting=final_plotting)

    def LayupFiberObject(self, DM, layupdirection, final_plotting=False):
        """ Use Autofiber object to determine layup orientation based on a fiber angle. """
        if self.orientation is not None:
            texcoord2inplane = self.orientation.layup(layupdirection, plotting=final_plotting)

            layerfiber = DM.fiberinstrs.assign_variable(("%s_%s_fiber" % (self.name, layupdirection)).replace("-", "n"),
                                                        (texcoord2inplane,
                                                         self.orientation.vertices,
                                                         self.orientation.vertexids,
                                                         self.orientation.inplanemat,
                                                         self.orientation.boxes,
                                                         self.orientation.boxpolys,
                                                         self.orientation.boxcoords,
                                                         self.orientation.facetnormals))
            for body in self.gk_layer.BodyList:
                body.fiberorientation = DM.autofiber.CreateFromParams(body.Name, DM.FEModel, layerfiber)
                body.fiberorientation.getMeshCenters()
                body.fiberorientation.getFiberOrientations()
                body.fiberorientation.CreateDiscreteField()
        else:
            raise ValueError("Fiber orientation is None. Must run CreateFiberObject on this layer first.")

    def Finalize(self,DM):
        """Build Python and Finite Element structure for this layer based on geometry kernel structure.
        It is not permissible to break a layer into pieces (multiple parts or bodies) after this 
        finalize call. It IS permissible to do surface imprints that operate in-place on the existing 
        layerbodies. """
        assert(self.parts is None) # this will fail if we try to finalize a second time
        self.parts=collections.OrderedDict()

        if self.coordsys is None:
            self.coordsys=SimpleCoordSys((1.0,0.0,0.0),(0.0,1.0,0.0))
            pass
        
        for gk_layerbody in self.gk_layer.BodyList:

            part=LayerPart.FromGK3D(DM,gk_layerbody)  # This also adds it to to_be_saved

            #part.CreateInstance(dependent=abqC.ON)

            part.ApplyLayup(self.coordsys,self.layupdirection,gk_layerbody.fiberorientation)
            part.AssignSection(self.LayerSection)
            self.parts[gk_layerbody.Name]=part
            
            DM.to_be_saved_dependencies.append(self.gk_layer) # prevent the gk_layer from being destroyed 
            
            #self.DM.to_be_saved.append(gk_part)  # mark the geometry layer object as to be saved to the .SAT file
            pass
    
        
        pass


    @classmethod
    def CreateFromMold(cls,DM,mold,direction,thickness,name,Section,layup,meshsize=0.5,coordsys=None):
        """Create a layer atop the specified mold. 
 * direction: "OFFSET" or "ORIG"
 * thickness: Thickness of layer (offsetting operation)
 * name: Unique name for layer
 * Section: ABAQUS section fo the layer
 * layup: Ply orientation in degrees
 * coordsys: Reference coordinate system for layup"""

        gk_layer=OCCLayer.CreateFromMold(name,mold,thickness,direction,DM.modelbuilder.PointTolerance,MeshSize=meshsize)
        
        return cls(name=name,gk_layer=gk_layer,layupdirection=layup,LayerSection=Section,coordsys=coordsys)


    
    pass





class LaminaContact(object):
    """LaminaContact represents the Contact/Adhesion/Cohesion between two layers. It includes 
storage of the various tie/contact/cohesive objects tying the two layers together"""
    
    DM=None
    """Reference to the DelamoModeler object"""
    
    bottomlamina=None
    """Bottom layer"""
    
    toplamina=None
    """Top layer"""
    
    cohesives=None
    """list of contact objects attaching these two layers"""
    
    contacts=None
    """list of contact objects between these two layers"""

    ties=None
    """list of Tie objects attaching these two layers"""


    # provide type of bottomlamina and toplamina to code generator so
    # we can simulate introspection

    def __init__(self,**kwargs):
        # LaminaContact object is created in bond_layers() 
        
        self.cohesives=[]
        self.contacts=[]
        self.ties=[]
        
        for key in kwargs:
            assert(hasattr(self,key))
            setattr(self,key,kwargs[key])
            pass
        
        
        pass


    def DefineContinuity(self,FEModel,bottom_laminapart,top_laminapart,point_normal,pointtolerance,normaltolerance,master=None):
        """Define a TIE (continuity) boundary condition between a face of a part of the bottom lamina and a corresponding face 
        of a part of the top lamina"""
        
        # Tie B.C.
        # Can specify master=top_bodynum or master=bottom_bodynum 
        # to force which element is master
        
        topface=top_laminapart.GetInstanceFace_point_normal(point_normal,
                                                            pointtolerance,
                                                            normaltolerance)
        bottomface=bottom_laminapart.GetInstanceFace_point_normal(point_normal,
                                                                  pointtolerance,
                                                                  normaltolerance)
        
        # Could create surfaces here
        # topsurface=Laminate.Surface(name="%sTop" % (self.toplamina.partname),side1Faces=topfaces)
        # bottomsurface=Laminate.Surface(name="%sBottom" % (self.bottomlamina.partname),side1Faces=bottmfaces)
        
        if master is not None: 
            # Note that this creates a "Surface-like" region, not the "Set-like" regions used in the Get...Region() functions above. See Abaqus Scripting Reference Guide section 45.3 for the distinction
            topregion=self.DM.regionToolset.Region(side1Faces=topface)
            bottomregion=self.DM.regionToolset.Region(side1Faces=bottomface)
            pass
            
                        
        # If the faces are subdivided on one side, 
        # have found that convergence is improved by using
        # the little subdivided faces as masters
        # and the big unified face as slave
        
        # note that right now, they really shouldn't be subdivided because they came from a single point and normal
        name="LaminaContact_%s_%s_Continuity_%d" % (bottom_laminapart.name,top_laminapart.name,self.DM.get_unique())


        if master is top_laminapart:
            tie=FEModel.Tie(name=name,
                            master=topregion,
                            slave=bottomregion,
                            positionToleranceMethod=abqC.COMPUTED,
                            adjust=abqC.ON,
                            tieRotations=abqC.ON,
                            thickness=abqC.ON)
            pass
        elif master is bottom_laminapart:
            # top is big unified face
            tie=FEModel.Tie(name=name,
                            master=bottomregion,
                            slave=topregion,
                            positionToleranceMethod=abqC.COMPUTED,
                            adjust=abqC.ON,
                            tieRotations=abqC.ON,
                            thickness=abqC.ON)
            pass
        else:
            assert(master is None)
            
            # Need to make decision on which one is master based on
            # numbers of faces -- This information is not available
            # until abaqus executes, so we have to make the decision
            # in Abaqus context
            tie=self.DM.globals["TieLayersWithMasterHavingMoreFaces"](FEModel,name,topface,bottomface)
            
            
            pass
        self.ties.append(tie)
        pass
    
        
    def DefineLamination(self,FEModel,CohesiveInteraction,bottom_laminapart,top_laminapart,point_normal,pointtolerance,normaltolerance,master=None):
        """Define a cohesive boundary condition between a face of a part of the bottom lamina and a corresponding face 
        of a part of the top lamina"""
        # Find top surfaces
        # Could also use a 'tie' constraint if we
        # don't want to allow this to delaminate
        
        topface=top_laminapart.GetInstanceFace_point_normal(point_normal,
                                                            pointtolerance,
                                                            normaltolerance)
        bottomface=bottom_laminapart.GetInstanceFace_point_normal(point_normal,
                                                                  pointtolerance,
                                                                  normaltolerance)
        
        nameindex="LaminaContact_%s_%s_Cohesive_%d" % (top_laminapart.name,bottom_laminapart.name,self.DM.get_unique())
        # Could create surfaces here
        # topsurface=Laminate.Surface(name="%sTop" % (self.toplamina.partname),side1Faces=topfaces)
        # bottomsurface=Laminate.Surface(name="%sBottom" % (self.bottomlamina.partname),side1Faces=bottmfaces)
        
        if master is not None:
            # Note that this creates a "Surface-like" region, not the "Set-like" regions used in the Get...Region() functions above. See Abaqus Scripting Reference Guide section 45.3 for the distinction
            topregion=self.DM.regionToolset.Region(side1Faces=topface)
            bottomregion=self.DM.regionToolset.Region(side1Faces=bottomface)
            pass
            
        
        #assert(numtopfaces==1 || numbottomfaces==1) # both sides can't be subdivided
        
        if master is top_laminapart:
            cohesive=FEModel.SurfaceToSurfaceContactStd(adjustMethod=abqC.NONE,
                                                        clearanceRegion=None,
                                                        createStepName="Initial",
                                                        datumAxis=None,
                                                        initialClearance=abqC.OMIT,
                                                        interactionProperty=CohesiveInteraction.name,
                                                        master=topregion,
                                                        slave=bottomregion,
                                                        name=nameindex,
                                                        sliding=abqC.SMALL,
                                                        thickness=abqC.ON)
            pass
        elif master is bottom_laminapart:
            # top is big unified face
            cohesive=FEModel.SurfaceToSurfaceContactStd(adjustMethod=abqC.NONE,
                                                        clearanceRegion=None,
                                                        createStepName="Initial",
                                                        datumAxis=None,
                                                        initialClearance=abqC.OMIT,
                                                        interactionProperty=CohesiveInteraction.name,
                                                        master=bottomregion,
                                                        slave=topregion,
                                                        name=nameindex,
                                                    sliding=abqC.SMALL,
                                                        thickness=abqC.ON)
            pass
        else:
            assert(master is None)
            # If the faces are subdivided on one side, 
            # have found that convergence is improved by using
            # the little subdivided faces as masters
            # and the big unified face as slave
        
            # Need to make decision on which one is master based on
            # numbers of faces -- This information is not available
            # until abaqus executes, so we have to make the decision
            # in Abaqus context
            
            cohesive=self.DM.globals["CohesiveWithMasterHavingMoreFaces"](FEModel,nameindex,CohesiveInteraction,topface,bottomface)
            pass
                
        self.cohesives.append(cohesive)
        
        pass



    def DefineDelamination(self,FEModel,ContactInteraction,bottom_laminapart,top_laminapart,point_normal,pointtolerance,normaltolerance,master=None):
        """Define a contact boundary condition between a face of a part of the bottom lamina and a corresponding face 
        of a part of the top lamina"""

        topface=top_laminapart.GetInstanceFace_point_normal(point_normal,
                                                            pointtolerance,
                                                            normaltolerance)
        bottomface=bottom_laminapart.GetInstanceFace_point_normal(point_normal,
                                                                   pointtolerance,
                                                                   normaltolerance)
        
        # Could create surfaces here
        # topsurface=Laminate.Surface(name="%sTop" % (self.toplamina.partname),side1Faces=topfaces)
        # bottomsurface=Laminate.Surface(name="%sBottom" % (self.bottomlamina.partname),side1Faces=bottmfaces)
        
        if master is None:
            topregion=self.DM.regionToolset.Region(side1Faces=topface)
            bottomregion=self.DM.regionToolset.Region(side1Faces=bottomface)
            pass
        
        
        nameindex="LaminaContact_%s_%s_Contact_%d" % (top_laminapart.name,bottom_laminapart.name,self.DM.get_unique())

        
        if master is top_laminapart:
            contact=FEModel.SurfaceToSurfaceContactStd(adjustMethod=abqC.NONE,
                                                       clearanceRegion=None,
                                                       createStepName='Initial',
                                                       datumAxis=None,
                                                       initialClearance=abqC.OMIT,
                                                       interactionProperty=ContactInteraction.name,
                                                       master=topregion,
                                                       slave=bottomregion,
                                                       name=nameindex,
                                                       #sliding=abqC.FINITE or abqC.SMALL
                                                       sliding=abqC.FINITE,
                                                       thickness=abqC.ON)
            
            pass
        elif master is bottom_laminapart:
            contact=FEModel.SurfaceToSurfaceContactStd(adjustMethod=abqC.NONE,
                                                       clearanceRegion=None,
                                                       createStepName='Initial',
                                                       datumAxis=None,
                                                       initialClearance=abqC.OMIT,
                                                       interactionProperty=ContactInteraction.name,
                                                       master=bottomregion,
                                                       slave=topregion,
                                                       name=nameindex,
                                                       #sliding=abqC.FINITE or abqC.SMALL
                                                       sliding=abqC.FINITE,
                                                       thickness=abqC.ON)
            pass
        else:
            assert(master is None)
            # If the faces are subdivided on one side, 
            # have found that convergence is improved by using
            # the little subdivided faces as masters
            # and the big unified face as slave
            
            # note that right now, they really shouldn't be subdivided because they came from a single point and normal

            contact=self.DM.globals["ContactWithMasterHavingMoreFaces"](FEModel,nameindex,ContactInteraction,topface,bottomface)
                
            pass
        self.contacts.append(contact)
        pass
        
        
    pass


def shell_and_cutout_from_shelltool(DM,shelltool_filename):
    """ Load in a .SAT file with two CAD objects. Both are presumed to be sheet bodies (shells)
    They should both be colored. Exactly one of them should have a red component greater than 
    0.75. The one without a red component greater than 0.75 is the "shell"; the one with 
    a red component greater than 0.75 is the "tool", which is assumed to intersect with the
    shell, yielding a closed curve. This function performs the intersection of shell and tool, 
    and makes a shell with a hole, plus a cut-out that fits in the hole. This function
    returns a Part reprepresenting the shell with hole, a second Part representing the 
    cut-out, and a list of (point,tangent) tuples for identifying the edges around 
    the hole (useful for mesh refinement)"""

    # ***!!!! Needs update
    layer_mold_list=delamo.CADwrap.LayerMoldList()
    DM.modelbuilder.create_shell_cutout(shelltool_filename,layer_mold_list)

    gk_shell = layer_mold_list[1]
    gk_cutout = layer_mold_list[0]


    # identify edges... e.g. for mesh refinement
    shell_edge_point_list = delamo.CADwrap.TPoint3dList()
    shell_edge_tangent_list = delamo.CADwrap.TPoint3dList()
    shell_edge_normal_list = delamo.CADwrap.TPoint3dList()

    # Note: Horribly load_shell_model() is horribly misnamed. It iterates around the hole in the shell, filling out the edge point, tangent, and normal lists
    # NOTE: It will eventually need a cutout parameter to distinguish which hole in the shell to iterate over. 
    DM.modelbuilder.load_shell_model(gk_shell, shell_edge_point_list, shell_edge_tangent_list, shell_edge_normal_list)

    #edgeinfo = [ (shell_edge_point_list[edgecnt],shell_edge_tangent_list[edgecnt]) for edgecnt in range(len(shell_edge_point_list)) ] 

    return (Part.FromGK3D(DM,gk_shell,shell=True),Part.FromGK3D(DM,gk_cutout,shell=True,omit_from_FE=True),list(shell_edge_point_list))


class shell_solid_coupling(object):
    """This class stores the various parameters required to implement shell-to-solid coupling"""
    
    shell=None
    """Shell representation with hole, the first Part returned from shell_and_cutout_from_shelltool()"""
    
    cutout=None
    """Cutout representation, the second part returned from shell_and_cutout_from_shelltool()"""
    
    shell_edge_point_list=None
    """List of points on the edges of the shell cutout"""
    
    shell_edge_tangent_list=None
    """List of tangents on the edges of the shell cutout"""
    
    shell_edge_normal_list=None
    """List of normals on the edges of the shell cutout (in the plane of the shell, will be normal to the edge of the layer).""" 

    surfacefaces=None
    """List of lists of ABAQUS side1faces type face sets, from all shell/cutout edges"""

    def __init__(self,**kwargs):
        for key in kwargs:
            assert(hasattr(self,key))
            setattr(self,key,kwargs[key])
            pass
        pass

    def bond_layer(self,DM,layer,layeroffset):
        """Bond a layer to the shell"""
        # Apply layeroffset to get the layer boundary location

        # ***!!! Needs rework ***!!!
        LayerEdgePointList=delamo.CADwrap.TPoint3dList()
        DM.modelbuilder.translate_shell_edge_points(self.shell_edge_point_list,self.shell_edge_normal_list,layeroffset,LayerEdgePointList)
        
        layerSideFacePointList = delamo.CADwrap.TPoint3dList()
        layerSideFaceNormalList = delamo.CADwrap.TPoint3dList()
        DM.modelbuilder.find_side_faces(layer.gk_layer,LayerEdgePointList,layerSideFacePointList,layerSideFaceNormalList)
        
        
        for edgecnt in range(len(self.shell_edge_point_list)):
            solidface=DM.globals["GetFace_point_normal"](layer.singlepart.fe_inst.faces,(layerSideFacePointList[edgecnt],layerSideFaceNormalList[edgecnt]),DM.pointtolerance,DM.normaltolerance)
            self.surfacefaces[edgecnt].append(solidface)
            pass
        pass

    def Finalize(self,DM,influenceDistance=None,positionTolerance=None):
        """Create the ABAQUS ShellSolidCoupling object."""
        # Do we need to create a surface from 
        # DM.Surface(side1faces=...,name=...)?

        for edgecnt in range(len(self.shell_edge_point_list)):

            # Add together (union) all of the solid faces
            # corresponding to this edge region
            solidfaces=self.surfacefaces[edgecnt][0]
            for solidface in self.surfacefaces[edgecnt][1:]:
                solidfaces=solidfaces+solidface
                pass

            SolidFaceRegion=DM.regionToolset.Region(side1Faces=solidfaces)
            #SolidSurf=Laminate.Surface(name="SolidSurf",side1Faces=solidfaces)

            #ShellGeomPartBC=DM.bcinstrs.rewrapobj(ShellGeomPart)

            shelledge=DM.globals["GetEdge_point_tangent"](self.shell.fe_inst.edges,self.shell.fe_part.vertices,(tuple(self.shell_edge_point_list[edgecnt]),self.shell_edge_tangent_list[edgecnt]),0.1,DM.normaltolerance)
            ShellRegion=DM.regionToolset.Region(side1Edges=shelledge)
        
            #ShellSurf=Laminate.Surface(name="ShellSurf",side1Edges=shelledges)
            ssc_params={
                "name": 'ssc%d-%d' % (DM.get_unique(),edgecnt),
                "shellEdge": ShellRegion,
                "solidFace": SolidFaceRegion,
                "positionToleranceMethod": abqC.COMPUTED,
                "influenceDistanceMethod": abqC.DEFAULT
            }
        
            #DM.FEModel.ShellSolidCoupling(name='ssc%d-%d' % (DM.get_unique(),edgecnt),
            #                              shellEdge=ShellRegion,
            #                              solidFace=SolidFaceRegion,
            #                              positionToleranceMethod=abqC.COMPUTED,
            #                              influenceDistanceMethod=abqC.DEFAULT)
            
            if influenceDistance is not None:
                ssc_params["influenceDistanceMethod"]=abqC.SPECIFIED
                ssc_params["influenceDistance"]=influenceDistance
                pass
                
            if positionTolerance is not None:
                ssc_params["positionToleranceMethod"]=abqC.SPECIFIED
                ssc_params["positionTolerance"]=positionTolerance
                pass
            DM.FEModel.ShellSolidCoupling(**ssc_params)
            pass
        pass
    
    @classmethod
    def from_shell_and_cutout(cls,DM,shell,cutout):
        """Create shell_solid_coupling object from a shell and cutout, for example as returned by shell_and_cutout_from_shelltool()"""

        # ***!!! NEEDS UPDATE
        shell_edge_point_list = delamo.CADwrap.TPoint3dList()
        shell_edge_tangent_list = delamo.CADwrap.TPoint3dList()
        shell_edge_normal_list = delamo.CADwrap.TPoint3dList()
    
        # Note:  load_shell_model() is horribly misnamed. It iterates around the hole in the shell, filling out the edge point, tangent, and normal lists
        # NOTE: It will eventually need a cutout parameter to distinguish which hole in the shell to iterate over. 
        DM.modelbuilder.load_shell_model(shell.gk_part, shell_edge_point_list, shell_edge_tangent_list, shell_edge_normal_list)
        
        surfacefaces=[]

        # Fill surfacefaces with empty lists
        for edgecnt in range(len(shell_edge_point_list)):
            surfacefaces.append([])
            pass
            
        ssc = cls(shell=shell,
                  cutout=cutout,
                  shell_edge_point_list=shell_edge_point_list,
                  shell_edge_tangent_list=shell_edge_tangent_list,
                  shell_edge_normal_list=shell_edge_normal_list,
                  surfacefaces=surfacefaces)
        

        return ssc
    pass
        
# !!!*** the parameters of the bond_layers() call must be kept in sync
# with the bond_layers_params and bond_layers_default_params
# variables in processor.py
# *** ALSO NEED TO CHANGE CODE in processor.py/annotate_bond_layers_calls
# *** WHERE "basename" and "phase" parameters
def bond_layers(DM,layer1,layer2,defaultBC="TIE",delamBC="CONTACT",delamRingBC="NONE",CohesiveInteraction=None,ContactInteraction=None,delaminationlist=None,master_layer=None,cohesive_layer=None,delamo_sourceline=None,delamo_phase=None,delamo_basename=None):
    """Bond two layers together. Parameters:
* DM: DelamoModeler object
* layer1: First layer
* layer2: Second layer
* defaultBC: The boundary condition for the bonded zone: Generally "TIE", "COHESIVE", or "COHESIVE_LAYER"
* delamBC: The boundary condition for the bulk of the delaminated region(s). Generally  "CONTACT"
* delamRingBC: The boundary condition for the outer zone of the delaminated region: Generally "NONE"
* CohesiveInteraction: The ABAQUS interaction property for any cohesive portions of the bond
* ContactInteraction: The ABAQUS interaction property for any contact portions of the bond
* delaminationlist: A list of files with delamination outlines (loops of 3D coordinates)
* master_layer: A preference for which layer should be the master layer in ABAQUS
* cohesive_layer: If defaultBC=="COHESIVE_LAYER", this layer will be used to achieve the bond. The cohesive_layer should have finite thickness but should NOT be finalized (this routine will finalize the cohesive_layer). You will also need to mesh the cohesive_layer afterward using MeshCohesive()
* delamo_sourceline: An index of the line number of this bond_layers() call, used to distinguish between mulitiple bonding steps. 
* delamo_phase: Phase of this multi-step defect insertion process. 
* delamo_basename: Base directory name for creation of layer boundary (STL) meshes. 
"""
    # *** Should we be able to specify how fine the .stl mesh should be? ***!!!

    
    
    if defaultBC=="COHESIVE_LAYER":
        # This is really two bonding operations one between layer1 and cohesive_layer,
        # and one between layer2 and cohesive_layer

        if cohesive_layer is None:
            raise ValueError("When bonding layers with a cohesive layer, you must have explicitly created the cohesive layer and passed it as the cohesive_layer parameter")

        if delaminationlist is not None:
        
            # Split cohesive layer according to the delaminations
            for cnt in range(len(delaminationlist)):
                cohesive_layer.Split(delaminationlist[cnt],DM.abqpointtolerance)
                pass
            pass

        # Finalize the cohesive layer now that it has been split
        cohesive_layer.Finalize(DM)

        #import pdb
        #pdb.set_trace()

        #if delaminationlist is not None:
        #    DM.modelbuilder.apply_delaminations(layer1.gk_layer,cohesive_layer.gk_layer,delaminationlist) # Imprint faces on both sides, 
        #    pass
        
        # Now that delaminations have been used to break the cohesive layer into pieces, imprint the broken-down cohesive layer onto layer1
        DM.modelbuilder.imprint_layers(layer1.gk_layer,cohesive_layer.gk_layer)

        # return adjacent layers in face_adjacency_list
        face_adjacency_list = DM.modelbuilder.adjacent_layer_boundary_conditions(layer1.gk_layer,cohesive_layer.gk_layer,bc_map={ "TIE": defaultBC, "CONTACT": delamBC, "NONE": delamRingBC })  
        
        # Find delamination region from FAL so cohesive_layer is removed in this region
        for face_adjacency in face_adjacency_list:
            if face_adjacency["bcType"]=="NONE" or face_adjacency["bcType"]=="CONTACT":
                name_to_remove=face_adjacency["name2"] # name2 because cohesive_layer was 2nd parameter to adjacent_layer_boundary_conditions
                cohesive_layer_bodynames=[ layerbody.Name for layerbody in cohesive_layer.gk_layer.BodyList ]
                for cnt in range(len(cohesive_layer_bodynames)):
                    if cohesive_layer_bodynames[cnt]==name_to_remove:
                        # Remove body from layer. This is OK because Layers are mutable
                        del cohesive_layer.gk_layer.BodyList[cnt] 
                        break
                    pass
                pass
            pass
        
        # Now call ourselves to create TIE bond between layer1 and cohesive_layer, except in the delaminated region

        # BCTypes have already been written to layer1 and cohesive_layer by adjacent_layer_boundary_conditions() call above
        bond_layers(DM,layer1,cohesive_layer,defaultBC="TIE",delamBC=delamBC,delamRingBC=delamRingBC,CohesiveInteraction=CohesiveInteraction,ContactInteraction=ContactInteraction,delaminationlist=None,master_layer=layer1,cohesive_layer=None,delamo_sourceline=None,delamo_phase=None,delamo_basename=None)

        
        
        # Now call ourselves to create TIE bond between cohesive_layer and layer2, even in the delaminated region
        # because CONTACT/NONE will exist between layer1 and cohesive_layer
        bond_layers(DM,cohesive_layer,layer2,defaultBC="TIE",delamBC="TIE",delamRingBC="TIE",CohesiveInteraction=CohesiveInteraction,ContactInteraction=ContactInteraction,delaminationlist=None,master_layer=layer2,cohesive_layer=None,delamo_sourceline=None,delamo_phase=None,delamo_basename=None)
        
        return

    if cohesive_layer is not None:
        raise ValueError("cohesive_layer may only be set when defaultBC==BC_DEFAULT_COHESIVE_LAYER")
    
        
    #face_adjacency_list = delamo.CADwrap.FAL()
    #if delaminationlist is None:
    #face_adjacency_list = DM.modelbuilder.adjacent_layer_boundary_conditions(layer1.gk_layer,layer2.gk_layer,defaultBC) # Imprint faces on both sides, return adjacent layers in face_adjacency_list
    #    pass
    #else:

    if delaminationlist is not None:
        DM.modelbuilder.apply_delaminations(layer1.gk_layer,layer2.gk_layer,delaminationlist)
        pass
    

    face_adjacency_list = DM.modelbuilder.adjacent_layer_boundary_conditions(layer1.gk_layer,layer2.gk_layer,bc_map={ "TIE": defaultBC, "CONTACT": delamBC, "NONE": delamRingBC })
    

    # Can not bond non-existant objects
    ThisContact = LaminaContact(DM=DM,bottomlamina=layer1, toplamina=layer2)
    
    print("Face adjacency: Layer %s and Layer %s:" % (layer1.gk_layer.Name,layer2.gk_layer.Name))
    
    for face_adjacency in face_adjacency_list:
        # NOTE: We only use point1, vector1 because
        # the faces on both sides HAVE to line upt, therefore 
        # a single (point, normal) must work for both
        
        if master_layer is not None:
            if master_layer is layer1:
                master_part=layer1.parts[face_adjacency['name1']]
                pass
            else:
                master_part=layer2.parts[face_adjacency['name2']]
                pass
            pass
        else:
            master_part=None
            pass
        
        if face_adjacency['bcType'] == "COHESIVE":
            print("    Cohesive, Body %s to %s" %(face_adjacency['name1'],face_adjacency['name2']))
            assert(CohesiveInteraction is not None) # Must provide CohesiveInteraction parameter if such interaction is found between the two surfaces
            ThisContact.DefineLamination(DM.FEModel,
                                         CohesiveInteraction,
                                         layer1.parts[face_adjacency['name1']],
                                         layer2.parts[face_adjacency['name2']],
                                         (face_adjacency['point1'],
                                          face_adjacency['normal1']),
                                         DM.abqpointtolerance,
                                         DM.normaltolerance,
                                         master=master_part)
            
            pass
        
        elif face_adjacency['bcType'] == "CONTACT":
            print("    Contact, Body %s to %s" %(face_adjacency['name1'],face_adjacency['name2']))
            # Must provide CohesiveInteraction parameter if such interaction is found between the two surfaces
            assert(ContactInteraction is not None)
            ThisContact.DefineDelamination(DM.FEModel,
                                           ContactInteraction,
                                           layer1.parts[face_adjacency['name1']],
                                           layer2.parts[face_adjacency['name2']],
                                           (face_adjacency['point1'],
                                            face_adjacency['normal1']),
                                           DM.abqpointtolerance,
                                           DM.normaltolerance,
                                           master=master_part)
            pass
        elif face_adjacency['bcType'] == "TIE":
            print("    Tie, Body %s to %s" %(face_adjacency['name1'],face_adjacency['name2']))
            
            ThisContact.DefineContinuity(DM.FEModel,
                                         layer1.parts[face_adjacency['name1']],
                                         layer2.parts[face_adjacency['name2']],
                                         (face_adjacency['point1'],
                                           face_adjacency['normal1']),
                                         DM.abqpointtolerance,
                                         DM.normaltolerance,
                                         master=master_part)
            pass
        
        elif face_adjacency['bcType'] == "NONE":
            print("    Nomodel, Body %s to %s" %(face_adjacency['name1'],face_adjacency['name2']))
            pass
        else:
            print("    Other; bcType=%d, Body %s to %s" %(face_adjacency.bcType,face_adjacency['name1'],face_adjacency['name2']))
            raise ValueError(face_adjacency.bcType)

        if delamo_sourceline is not None:
            # import statement not at top of file
            # to avoid an import loop
            from delamo.process import output_dir
            
            output_directory=delamo.process.output_dir(delamo_basename)
            delamo_fname=os.path.join(output_directory,"layerboundary_%s_%5.5d.stl" % (delamo_phase,delamo_sourceline))
            DM.modelbuilder.save_layer_surface_stl(delamo_fname,layer1.gk_layer,layer2.gk_layer)
            pass
        
        
        pass

    pass
