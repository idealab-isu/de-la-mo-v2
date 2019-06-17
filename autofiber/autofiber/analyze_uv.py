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
import numpy as np

try:
    from cStringIO import StringIO  # python 2.x
    pass
except ImportError:
    from io import StringIO # python 3.x
    pass


from spatialnde.coordframes import concrete_affine
from spatialnde.coordframes import coordframe
from spatialnde.ndeobj import ndepart,ndeassembly
from spatialnde.cadpart.rectangular_plate import rectangular_plate
from spatialnde.cadpart.appearance import simple_material
from spatialnde.cadpart.appearance import texture_url
from spatialnde.cadpart.polygonalsurface_texcoordparameterization import polygonalsurface_texcoordparameterization
from spatialnde.exporters.x3d import X3DSerialization
from spatialnde.exporters.vrml import VRMLSerialization


# *** NOTE: We may want to either eliminate the "surface" abstraction
# for meshed boundary representations or develop a similar set of
# routines that cross surface boundaries so that overlapping textures
# can exist across surface boundaries ***

def IdentifyTexMaps(part,surfaceparameterizationmapping=None):    
    """ Trace through and identify all texture maps, given a part
    (ndepart instance)

    If texture_urls are given through an appearance node, those are used. If
    the appearance node is missing or does not provide a texture_url, then
    a numbered name is used for that surface, of the form _unnamed_surface_%d

    Returns a tuple of two dictionaries:
       surface_texurl is indexed by the id of the surface object and contains the the texture url as a string. 
       surfaces_bytexurl is indexed by texture url strings and contains a list of surface objects that share that texture url
"""

    
    surface_texurl={}  # texture url, indexed by id(surface)

    surfaces_bytexurl={} # list of surfaces, indexed by texurl

    surfacecnt=0
    for surface in part.implpart.surfaces:
        if surfaceparameterizationmapping is None:
            surfaceparameterization=surface.intrinsicparameterization
            pass
        else:
            surfaceparameterization=surfaceparameterizationmapping[id(surface)]
            pass
        
        if isinstance(surfaceparameterization,polygonalsurface_texcoordparameterization): # only extras texture coordinates
            
            # Currently only extras texture url
            if surface.appearance is not None and hasattr(surface.appearance,"texture_url"):
                surface_texurl[id(surface)]=surface.appearance.texture_url
                pass
            else:
                surface_texurl[id(surface)]="_unnamed_surface_%d" % (surfacecnt)
                pass

            if surface_texurl[id(surface)] not in surfaces_bytexurl:
                surfaces_bytexurl[surface_texurl[id(surface)]]=[]
                pass
            surfaces_bytexurl[surface_texurl[id(surface)]].append(surface)
            
            pass
        pass
    
    return (surface_texurl,surfaces_bytexurl)


def BuildEdgeDict(surface):
    """ Create edge dictionary from a surface.
    The edge dictionary is indexed by a tuple 
    (vertexindex1,vertexindex2) of indices into surface.vertexes. 
    It contains a list of polygon ids that have an edge that shares
    these two vertices. 
    
    This function assumes that identical vertices in the surface 
    have been merged, so the vertexindex uniquely identifies the
    vertex. 

    This function returns the edge dictionary
"""

    edges={}  # lists of polygon numbers, indexed by a tuple of the  vertex indices connected by the edge
    for polynum in range(surface.vertexidx_indices.shape[0]):
        firstidx=surface.vertexidx_indices[polynum]
                
        numvertices=surface.numvertices[polynum]
        for firstvertex in range(numvertices):
            nextvertex=(firstvertex+1) % numvertices
            
            firstvertexidx=surface.vertexidx[firstidx+firstvertex]
            nextvertexidx=surface.vertexidx[firstidx+nextvertex]
            
            if (firstvertexidx,nextvertexidx) in edges:
                edges[(firstvertexidx,nextvertexidx)].append((polynum,firstvertex,nextvertex))
                pass
            elif (nextvertexidx,firstvertexidx) in edges:
                edges[(nextvertexidx,firstvertexidx)].append((polynum,nextvertex,firstvertex))
                pass
            else:
                edges[(firstvertexidx,nextvertexidx)] = [ (polynum,firstvertex,nextvertex) ]
                pass
            pass
        pass
    return edges


def DetermineAdjacency(surface,edges,surfaceparameterization=None,texture=False):
    """ Build an adjacency index for the given surface 
        with the given edgedict. The adjacency index
        has the same layout as surface.vertexidx, but
        additional entries may be -1 as facets may have
        fewer adjacencies than vertices. 

    The adjacency index contains the polygon numbers adjacent
    to the given polygon. 

    If texture=True is given as a parameter, then the adjacency
    index built up will only show adjacencies both in the polygon 
    mesh and in the texture. In that case, the polygon indices
    will include both original polygons and redundant copies, 
    so the polygon indices will range from 0 to 
    vertexidx_indices.shape[0]+texcoordredundant_polystartindexs.shape[0]

"""
    if surfaceparameterization is None:
        surfaceparameterization=surface.intrinsicparameterization
        pass

    if texture:
        idxsize = surface.vertexidx.shape[0]
        if surfaceparameterization.texcoordredundant_texcoordidx is not None:
            idxsize += surfaceparameterization.texcoordredundant_texcoordidx.shape[0]
        pass
        adjacencyidx = np.ones(idxsize, dtype=np.int32) * -1
        pass
    else:
        adjacencyidx = np.ones(surface.vertexidx.shape[0],dtype=np.int32)*-1
        pass
    
    
    for (vertex1,vertex2) in edges.keys():
        numedges=len(edges[(vertex1,vertex2)])
        assert(numedges <= 2) # If this assertion fails, then there was an edge between polygons that was shared between more than two polygons (?)
        if numedges == 2:

            ((polynum1,polynum1v1,polynum1v2),(polynum2,polynum2v1,polynum2v2)) = edges[(vertex1,vertex2)]
            #textureadjacent=False
            
            if texture: 
                # Also check if the texcoords are adjacent
                # for these edges
                # Texture coordinates for polygon #1
                firstidx1=surface.vertexidx_indices[polynum1]
                # 1st vertex (a) of polygon 1
                texcoord1a=[ surfaceparameterization.texcoord[surfaceparameterization.texcoordidx[firstidx1+polynum1v1],:] ]
                # 2nd vertex (b) of polygon 1
                texcoord1b=[ surfaceparameterization.texcoord[surfaceparameterization.texcoordidx[firstidx1+polynum1v2],:] ]
                extpolynum1=[ polynum1 ]
                
                # redundant texture mappings of polygon 1
                if surfaceparameterization.texcoordredundant_firstpolynum is not None:
                    firstpolynum1=surfaceparameterization.texcoordredundant_firstpolynum[polynum1]
                    numcopies1=surfaceparameterization.texcoordredundant_numcopies[polynum1]
                    for texcnt1 in range(numcopies1):
                        startindex1=surfaceparameterization.texcoordredundant_polystartindexes[firstpolynum1+texcnt1]
                        # identify the polygon index of this redundant texture copy
                        
                        extpolynum1.append(surface.vertexidx_indices.shape[0] + startindex1)
                        # 1st vertex (a)
                        texcoord1a.append(surfaceparameterization.texcoord[surfaceparameterization.texcoordredundant_texcoordidx[startindex1+polynum1v1],:])
                        # 2nd vertex (b) of polygon 1
                        texcoord1b.append(surfaceparameterization.texcoord[surfaceparameterization.texcoordredundant_texcoordidx[startindex1+polynum1v2],:])

                        pass
                    pass
                
                
                # Texture coordinates for polygon #2
                firstidx2=surface.vertexidx_indices[polynum2]
                # 1st vertex (a) of polygon 2
                texcoord2a=[ surfaceparameterization.texcoord[surfaceparameterization.texcoordidx[firstidx2+polynum2v1],:] ]
                # 1st vertex (b) of polygon 2
                texcoord2b=[ surfaceparameterization.texcoord[surfaceparameterization.texcoordidx[firstidx2+polynum2v2],:] ]
                extpolynum2=[ polynum2 ]


                # redundant texture mappings of polygon #2
                if surfaceparameterization.texcoordredundant_firstpolynum is not None:
                    firstpolynum2=surfaceparameterization.texcoordredundant_firstpolynum[polynum2]
                    numcopies2=surfaceparameterization.texcoordredundant_numcopies[polynum2]
                    for texcnt2 in range(numcopies2):
                        startindex2=surfaceparameterization.texcoordredundant_polystartindexes[firstpolynum2+texcnt2]
                        # identify the polygon index of this redundant texture copy
                        extpolynum2.append(surface.vertexidx_indices.shape[0] + startindex2)

                        # 1st vertex (a)
                        texcoord2a.append(surfaceparameterization.texcoord[surfaceparameterization.texcoordredundant_texcoordidx[startindex2+polynum2v1],:])
                        # 2nd vertex (b) of polygon 2
                        texcoord2b.append(surfaceparameterization.texcoord[surfaceparameterization.texcoordredundant_texcoordidx[startindex2+polynum2v2],:])
                        pass
                    pass
                

                
                #if (((texcoord1a==texcoord2a).all() and (texcoord1b==texcoord2b).all()) or
                #    ((texcoord1a==texcoord2b).all() and (texcoord1b==texcoord2a).all())):
                #    # texcoords match on both sides of the edge
                #    textureadjacent=True
                #    pass
                for cnt1 in range(len(texcoord1a)):
                    for cnt2 in range(len(texcoord2a)):
                        if (((texcoord1a[cnt1] == texcoord2a[cnt2]).all() and (texcoord1b[cnt1] == texcoord2b[cnt2]).all()) or 
                            ((texcoord1a[cnt1] == texcoord2b[cnt2]).all() and (texcoord1b[cnt1] == texcoord2a[cnt2]).all())):
                            # ... Then there is an adjacency between
                            # extpolynum1[cnt1] and
                            # extpolynum2[cnt2]
                            # textureadjacent=True
                            
                            # Store adjacency in adjacencyidx
                            for (polynuma,polynumb) in ((extpolynum1[cnt1],extpolynum2[cnt2]),(extpolynum2[cnt2],extpolynum1[cnt1])):
                                adjacencycnt=0
                                if polynuma < surface.vertexidx_indices.shape[0]:
                                    adjacencyidxidx=surface.vertexidx_indices[polynuma]
                                    pass
                                else:
                                    adjacencyidxidx=surfaceparameterization.texcoordredundant_polystartindexes[polynuma-surface.vertexidx_indices.shape[0]]
                                    pass
                                
                                while adjacencyidx[adjacencyidxidx+adjacencycnt] >= 0:
                                    adjacencycnt += 1
                                    pass
                                adjacencyidx[adjacencyidxidx+adjacencycnt]=polynumb
                                pass

                            pass
                        pass
                    pass
                
                
                
                pass

            if not(texture):
                # Store adjacency in adjacencyidx
                for (polynuma,polynumb) in ((polynum1,polynum2),(polynum2,polynum1)):
                    adjacencycnt=0
                    while (adjacencyidx[surface.vertexidx_indices[polynuma]+adjacencycnt] >= 0):
                        adjacencycnt+=1
                        pass
                    adjacencyidx[surface.vertexidx_indices[polynuma]+adjacencycnt]=polynumb
                    pass
                pass
            pass
        pass
    
    return adjacencyidx

def FindTexPatches(surface,texadjacencyidx,surfaceparameterization=None):
    """ Given a surface and a texture adjacency index (which contains the texture polygon numbers adjacent, to the given texture polygon, with polygons findable from surface.vertexidx_indices and surfaceparameterization.texcoordredundant...), separate the polygons of the surface into groups that have adjacent texture. Return a list of lists of polygon numbers."""

    if surfaceparameterization is None:
        surfaceparameterization=surface.intrinsicparameterization
        pass
    
    # Now we can go through and find adjacent polygons
    # for which the edge texcoords also match
    UnusedPoly=np.ones(surface.vertexidx_indices.shape[0],dtype=np.bool)
    WhereUnused=np.arange(surface.vertexidx_indices.shape[0],dtype=np.int32)

    Patches=[]
    while len(WhereUnused) > 0:
        # Take first unused element
        polynum=WhereUnused[0]
        
        ThisPatch=[]
        BuildPatch(ThisPatch,UnusedPoly,surface.vertexidx_indices,texadjacencyidx, polynum )
        
        Patches.append(ThisPatch)
        WhereUnused=np.where(UnusedPoly)[0]
        pass
    
    

    return Patches

def BuildPatch(ThisPatch,UnusedPoly,vertexidx_indices,texadjacencyidx,startpolynum):
    # Build a patch (list of adjacent polygons) from the pool of
    # unused polygons, the indexes into the texadjacencyidx, and
    # the texadjacencyidx.
    # this is a formerly recursive function. Call it once with a single polynum
    # and ThisPatch as an empty list,
    # and it will fill up ThisPatch, removing the polynums from UnusedPoly

    polynums = [ startpolynum ]

    while len(polynums) > 0:
        polynum = polynums.pop()
        
        firstidx=vertexidx_indices[polynum]
    
        ThisPatch.append(polynum)
        UnusedPoly[polynum]=False
                
        neighborcnt=0
        while texadjacencyidx[firstidx+neighborcnt] >= 0:

            if UnusedPoly[texadjacencyidx[firstidx+neighborcnt]]:
                polynums.append(texadjacencyidx[firstidx+neighborcnt])
                pass
            neighborcnt+=1
            pass
        
        pass
    pass

def FindEdgeLoops(boundaryvertexpairlist):
    # boundaryvertexpairlist is a list of boundary vertex pairs for
    # this patch boundary.

    # Resolve into a list of sorted edge loops

    edgeloops=[]
    edgeloopcnt=0
    
    # Build dictionary by vertexnum of
    # lists of two vertexpairs
    vertexpairbyvertexdict={}
    for vertexpair in boundaryvertexpairlist:
        if vertexpair[0] not in vertexpairbyvertexdict:
            vertexpairbyvertexdict[vertexpair[0]]=[ vertexpair ]
            pass
        else:
            vertexpairbyvertexdict[vertexpair[0]].append(vertexpair)
            pass

        if vertexpair[1] not in vertexpairbyvertexdict:
            vertexpairbyvertexdict[vertexpair[1]]=[ vertexpair ]
            pass
        else:
            vertexpairbyvertexdict[vertexpair[1]].append(vertexpair)
            pass

        pass
    
    # Create a set with a copy of the pair list
    # pull a vertex from it, then walk the list
    # until it is empty
    boundaryvertexpairset=set(boundaryvertexpairlist)

    while len(boundaryvertexpairset) > 0:
        curpair=boundaryvertexpairset.pop()
        edgeloops.append([ curpair ])        

        firstvertex=curpair[0]
        lastvertex=curpair[1]
        while lastvertex != firstvertex:
            # Find a pair with lastvertex
            vertexpairs=vertexpairbyvertexdict[lastvertex]
            assert(len(vertexpairs)==2)
            if vertexpairs[0]==curpair:
                newpair=vertexpairs[1]
                pass
            else:
                assert(vertexpairs[1]==curpair)
                newpair=vertexpairs[0]
                pass

            # update lastvertex one step around the loop
            if newpair[0]==lastvertex:
                lastvertex=newpair[1]
                pass
            else:
                assert(newpair[1]==lastvertex)
                lastvertex=newpair[0]
                pass
            # update curpair

            curpair=newpair

            # Store this edge as part of the loop
            edgeloops[edgeloopcnt].append(curpair)

            # Loop terminates when lastvertex matches firstvertex
            # (i.e. when we close the loop)
            pass
        edgeloopcnt+=1
        pass
    return edgeloops

def PatchFindBoundaries(Patches,edges):
    # Patches is a list of patches (each patch is a list of polygons)
    # Determine the boundaries of a patch

    # Find all of the edges that are included by exactly one
    # of the polygons
    polygonsets = [ set(ThisPatch) for ThisPatch in Patches ]

    #edgepolycnt = {}  # dictionary by (vertexidx1,vertexidx2) (like edgedict) of array of per-patch counts
    boundaryvertexpairs = [ [] for ThisPatch in Patches ] # list by patch of lists of vertexidpairs for the boundary of a patch
    
    
    for vertexidxpair in edges.keys():
        EdgePolys=edges[vertexidxpair]
        assert(len(EdgePolys) > 0 and len(EdgePolys) <= 2)

        for cnt in range(len(Patches)):
            if EdgePolys[0] in polygonsets[cnt] and (len(EdgePolys) < 2 or EdgePolys[1] not in polygonsets[cnt]):
                # On the boundary of this patch
                boundaryvertexpairs[cnt].append(vertexidxpair)
                pass
            elif len(EdgePolys) == 2 and EdgePolys[1] in polygonsets[cnt] and EdgePolys[0] not in polygonsets[cnt]:
                # On the boundary of this patch
                boundaryvertexpairs[cnt].append(vertexidxpair)
                pass
                
            pass
        pass
    
    boundaryloops = [ FindEdgeLoops(boundaryvertexpairlist) for boundaryvertexpairlist in boundaryvertexpairs ] # list by patch of lists of vertexidpairs loops for the boundary of a patch
    
if __name__=="__main__":

    # define lab frame
    labframe = coordframe()

    # define camera frame
    camframe = coordframe()

    # define object frame
    objframe = coordframe()

    specimen=ndepart.fromx3d(objframe,None,"C17-NASAMODELCURVED_UV.x3d",tol=1e-6)
    (surface_texurl,surfaces_bytexurl)=IdentifyTexMaps(specimen)

    Patches={}
    # Now iterate over all of the unique textures
    for texurl in surfaces_bytexurl.keys():
        for surface in surfaces_bytexurl[texurl]:

            edges=BuildEdgeDict(surface)

            texadjacencyidx=DetermineAdjacency(surface,edges,texture=True)
            
            # now texadjacencyidx is like vertexidx but gives adjacent polygons that also have adjacent texture

            Patches[id(surface)]=FindTexPatches(surface,texadjacencyidx)
            
        
            pass
        
        pass
    
    pass

