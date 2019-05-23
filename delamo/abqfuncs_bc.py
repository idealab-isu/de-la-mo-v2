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


def UniqueNameForDict(dct,name_with_percent_d):
    cnt=0
    while (name_with_percent_d % (cnt)) in dct:
        cnt+=1
        pass
    return name_with_percent_d % (cnt)




def TieLayersWithMasterHavingMoreFaces(FEModel,nameindex,topfaces,bottomfaces):
    numtopfaces=CountFaces(topfaces)
    numbottomfaces=CountFaces(bottomfaces)
    
    # Note that this creates a "Surface-like" region, not the "Set-like" regions used in the Get...Region() functions in api.py. See Abaqus Scripting Reference Guide section 45.3 for the distinction
    topregion=regionToolset.Region(side1Faces=topfaces)
    bottomregion=regionToolset.Region(side1Faces=bottomfaces)

    
    if numtopfaces > numbottomfaces:
        # bottom is big unified face
        tie=FEModel.Tie(name=nameindex,
                        master=topregion,
                        slave=bottomregion,
                        positionToleranceMethod=abqC.COMPUTED,
                        adjust=abqC.ON,
                        tieRotations=abqC.ON,
                        thickness=abqC.ON)
        pass
    else:
        # top is big unified face
        tie=FEModel.Tie(name=nameindex,
                        master=bottomregion,
                        slave=topregion,
                        positionToleranceMethod=abqC.COMPUTED,
                        adjust=abqC.ON,
                        tieRotations=abqC.ON,
                        thickness=abqC.ON)
        
        pass
    return tie

def CohesiveWithMasterHavingMoreFaces(FEModel,nameindex,CohesiveInteraction,topfaces,bottomfaces):
    
    numtopfaces=CountFaces(topfaces)
    numbottomfaces=CountFaces(bottomfaces)
    
    topregion=regionToolset.Region(side1Faces=topfaces)
    bottomregion=regionToolset.Region(side1Faces=bottomfaces)
    
    if numtopfaces > numbottomfaces:
        # bottom is big unified face
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
    else:
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
    return cohesive



def ContactWithMasterHavingMoreFaces(FEModel,nameindex,ContactInteraction,topfaces,bottomfaces):
    
    numtopfaces=CountFaces(topfaces)
    numbottomfaces=CountFaces(bottomfaces)
    
    topregion=regionToolset.Region(side1Faces=topfaces)
    bottomregion=regionToolset.Region(side1Faces=bottomfaces)
    if (numtopfaces > numbottomfaces):
        # bottom is big unified face
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
    else:
        
        # top is big unified face
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
    return contact



def BCFaces_OBSOLETE(Parts,BodyNum_if_not_merged,points_normals,facepointtolerance,normaltolerance):
    """ Given Parts (could be a Lamina) and a BodyNumber (unless contains a single body number), a tuple of (surface point, normal)s, and tolerances to identify
    the faces, return the face object to use for a boundayr condition"""

    if BodyNum_if_not_merged is not None:
        instancefaces=Parts.bodynumbers[BodyNum_if_not_merged].instance.faces
        pass
    else:
        assert(len(Parts.parts)==1) # Need to provide a body number if lamina elements not merged
        instancefaces=Parts.parts[0].instance.faces
        pass
    
    return GetFaces(instancefaces,
                    points_normals,
                    facepointtolerance,
                    normaltolerance)

def BCEdges_OBSOLETE(Parts,BodyNum_if_not_merged,points_tangents,facepointtolerance,normaltolerance):
    """ Given Parts (could be a Lamina) and a BodyNumber (unless contains a single body number), a tuple of (surface point, normal)s, and tolerances to identify
    the faces, return the face object to use for a boundayr condition"""

    if BodyNum_if_not_merged is not None:
        instanceedges=Parts.bodynumbers[BodyNum_if_not_merged].instance.edges
        instancevertices=Parts.bodynumbers[BodyNum_if_not_merged].instance.vertices
        pass
    else:
        assert(len(Parts.parts)==1) # Need to provide a body number if lamina elements not merged
        instanceedges=Parts.parts[0].instance.edges
        instancevertices=Parts.parts[0].instance.vertices
        pass
    
    return GetEdges(instanceedges,
                    instancevertices,
                    points_tangents,
                    facepointtolerance,
                    normaltolerance)

