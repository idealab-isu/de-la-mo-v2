# Name of results surfaces --done--
# Single surface for loads and deflections --done--

# Does getNodes work on vertices and nodes
    # Nodes are the discrete mesh points, while vertices are geometric objects
    # This depends on what we might be using the set for
    # If we are extracting results we likely want to find nodes
    # If we are applying boundary conditions, we probably want to find vertices
    # Vertices use vert.pointOn, while nodes use node.coordinates

# Are there differences between node arrays and vertex arrays
    # Nodes are different than vertices and have different attributes

# Not create surface, reference surface (ReferenceSurfaceFromAssembly)
# Create set at part level or assembly level?
    # If the instance is made independent of the assembly, only geometric sets can be created
    # If the instance is dependent on the assembly, then you can create both geometric and mesh sets
    # Part sets apply only to the part, while assembly sets can span instances



# Where do we use vertices and where do we use nodes?
# Split create set into separate functions?
# Get faces from assembly using assembly.faces or assembly.instance.faces

# A set represents a collection of geometric or mesh objects
# A set can contain vertices, edges, faces, nodes, or elements

# If a set is defined by a vertex (or vertices), both the vertices and the mesh nodes at those
# vertices can be accessed

# If a set is defined by an edge (or edges), the edges can be accessed as well as the mesh nodes
# that fall on the edge

# If a set is defined by a face (or faces), the face can be accessed as well as the mesh nodes
# and elements that include the specified face

# If a set is defined by a node (or nodes), the mesh nodes can be accessed
# If a set is defined by an element (or elements), the mesh elements can be accessed as well
# as the nodes that are included in the elements.

# Any combination of these objects can make up a set, and the output will be a set that is the
# combination of these specific cases

# The same applies to surfaces, except surfaces cannot contain vertices or edges
# A geometric surface will contain the geometric faces, as well as meshed elements and nodes
# A mesh surface will contain only the elements and nodes

# For all intents and purposes, a "surface" seems to be a subset of a "set"
# i.e. anything we can describe with a surface we can describe with a set
# One difference is that surfaces can specify which side of a face it is referring to


# 3/18/16 Cheryl Rose suggests using Continuum shell element
# to consist of multiple lamina.
#  * Also suggests having option to define cohesive element rather than
#    cohesive surface
#  * disbonded regions should use contact model

MaskSubStrParseObj=re.compile(r"""#([0-9A-Fa-f]+)(?::(\d+))?""")
def CountFaces(facearray):
    Mask=facearray.getMask()
    #print("mask=%s\n\n\n" % (str(Mask)))
    assert(len(Mask)==1)  # should be a length-1 tuple
    MaskStr=Mask[0]
    assert(MaskStr.startswith("(\'[")) # check for brackets
    assert(MaskStr.endswith("]\',),"))
    MaskStr=MaskStr[3:-5] # strip off brackets

    numfaces=0
    for MaskSubStr in MaskStr.split():
        MaskSubStr=MaskSubStr.strip()
        if len(MaskSubStr)==0:
            continue
        MatchObj=MaskSubStrParseObj.match(MaskSubStr)
        if MatchObj is None:
            raise ValueError("Error matching mask substring \"%s\"" % (MaskSubStr))
        # MatchObj.group(1) is hex number, MatchObj.group(2) is repetition count
        LongMask=long(MatchObj.group(1),base=16)
        nreps=1  # number of repetitions
        if MatchObj.group(2) is not None:
            nreps=long(MatchObj.group(2)) # !!!*** If this errors out, may need to do conversion from base 16, not base 10 in long()
            pass

        # Convert to binary string representation; count 1's
        numfaces += bin(LongMask).count("1")*nreps
        
        pass
    return numfaces

    
def GetMaskStr(masknum):
    # Note: mask should be a long

    # masks are hexadecimal numbers. The bits of the
    # hexadecimal number represent the face indices.
    # The mask is broken into groups of 32 bits.
    # The first group represents the lowest-order 32 bits
    # The second group represents the second-lowest-order 32 bits
    # etc. (a shorthand notation for repeated groups uses
    # #0badfeed:repetitionnumber to repeat a 32 bit group more
    # than once

    #maskstr="[ #%x ]" % (masknum)
    
    maskstr="[ "
    while masknum != 0:
        maskstr += "#%x " % (masknum & 0xffffffff)
        masknum >>= 32
        pass
    
    maskstr+="]"
    return maskstr


def _get_face_indices_point_normal(facearray,facepointnormal,pointtolerance,normaltolerance):
    # IMPORTANT: Normals must be unit vectors!!!
    # facepointnormal is a (points,normal) tuple
    gotfaces=[]
    point=facepointnormal[0]
    normal=facepointnormal[1]
    #print("len(facearray)=%d" % (len(facearray)))
    #print("facearray=%s" % (str(facearray)))
    
    for facecnt in range(len(facearray)):
        # if facepointnormal[0][2]==1.5:
        #print("facecnt=%d" % (facecnt))
        #     pass
        actualnormal=None
        try: 
            actualnormal=facearray[facecnt].getNormal(point=point)
        except abq.AbaqusException:
            pass
        # if facepointnormal[0][2]==1.5:
        #print("actualnormal=%s" % (str(actualnormal)))
        #print("point=%s" % (str(point)))
        #print("normal=%s" % (str(normal)))
        #     pass
        if actualnormal is not None and np.abs(np.dot(actualnormal,normal))-1.0 > -normaltolerance:
            # normal matches... does point?
            # Create a mask containing this face alone, so we can
            # get a FaceArray containing only this face
            # masks are hexadecimal numbers. The bits of the
            # hexadecimal number represent the face indices.
            # The mask is broken into groups of 32 bits.
            # The first group represents the lowest-order 32 bits
            # The second group represents the second-lowest-order 32 bits
            # etc. (a shorthand notation for repeated groups uses
            # #0badfeed:repetitionnumber to repeat a 32 bit group more
            # than once
            # Our face index is facecnt
            masknum=(long(1) << facecnt)
            #maskstr="[ #%x ]" % (masknum)
            maskstr=GetMaskStr(masknum)
                
            #print("facecnt=%d; maskstr=%s" % (facecnt,maskstr))
            reduced_facearray=facearray.getSequenceFromMask(mask=(maskstr,))
            closest=reduced_facearray.getClosest(coordinates=(point,),searchTolerance=pointtolerance)
            # if facepointnormal[0][2]==1.5:
            #     print("closest=%s" % (repr(closest)))
            #     pass
            # closest is a dictionary, indexed by coordinate number (just 0)
            assert(len(list(closest.keys()))<=1)  # Should get no more than one result, since there is only one face involved
            if len(list(closest.keys())) > 0:
                # got a result
                assert(0 in closest) # only coordinate number is 0
                closestface=closest[0][0] # first element is face
                closestpoint=closest[0][1] # second element is point
                assert(closestface.index==facecnt) # should only get the face we are working on, unless getClosest() is like findAt() and might find the wrong face
                assert(np.linalg.norm(np.array(closestpoint)-np.array(point)) <= pointtolerance) # distance should be within tolerance
                gotfaces.append(closestface)
                pass
            
            pass
        #actualnormals.append(actualnormal)
        
        
        pass
        
    return gotfaces



def _get_face_indices_point(facearray,point,pointtolerance):
    # IMPORTANT: Normals must be unit vectors!!!
    # point is an (x,y,z) tuple
    gotfaces=[]
    #print("len(facearray)=%d" % (len(facearray)))
    #print("facearray=%s" % (str(facearray)))
    
    for facecnt in range(len(facearray)):

        # does point match?
        # Create a mask containing this face alone, so we can
        # get a FaceArray containing only this face
        # masks are hexadecimal numbers. The bits of the
        # hexadecimal number represent the face indices.
        # The mask is broken into groups of 32 bits.
        # The first group represents the lowest-order 32 bits
        # The second group represents the second-lowest-order 32 bits
        # etc. (a shorthand notation for repeated groups uses
        # #0badfeed:repetitionnumber to repeat a 32 bit group more
        # than once
        # Our face index is facecnt
        masknum=(long(1) << facecnt)
        #maskstr="[ #%x ]" % (masknum)
        maskstr=GetMaskStr(masknum)
                
        #print("facecnt=%d; maskstr=%s" % (facecnt,maskstr))
        reduced_facearray=facearray.getSequenceFromMask(mask=(maskstr,))
        closest=reduced_facearray.getClosest(coordinates=(point,),searchTolerance=pointtolerance)
        # if facepointnormal[0][2]==1.5:
        #     print("closest=%s" % (repr(closest)))
        #     pass
        # closest is a dictionary, indexed by coordinate number (just 0)
        assert(len(list(closest.keys()))<=1)  # Should get no more than one result, since there is only one face involved
        if len(list(closest.keys())) > 0:
            # got a result
            assert(0 in closest) # only coordinate number is 0
            closestface=closest[0][0] # first element is face
            closestpoint=closest[0][1] # second element is point
            assert(closestface.index==facecnt) # should only get the face we are working on, unless getClosest() is like findAt() and might find the wrong face
            assert(np.linalg.norm(np.array(closestpoint)-np.array(point)) <= pointtolerance) # distance should be within tolerance
            gotfaces.append(closestface)
            pass
        
        
        pass
        
    return gotfaces




def GetFaces_points_normals(facearray,facepointsnormals,pointtolerance,normaltolerance):
    # facepointsnormals is a sequence
    # of (point,normal) tuples
    # IMPORTANT: Normals must be unit vectors!!!
    
    allfaces=[]
    
    for pointcnt in range(len(facepointsnormals)):
        gotfaces=_get_face_indices_point_normal(facearray,facepointsnormals[pointcnt],pointtolerance,normaltolerance)
        if len(gotfaces)==0:
            print("facearray=%s; facepointnormal=%s" % (repr(facearray),str(facepointsnormals[pointcnt])))
            pass

        # Should be able to find exactly one face from each (point, normal)
        assert(len(gotfaces) > 0) # should be able to find faces
        assert(len(gotfaces) < 2) # should be able to find unique face

        allfaces.extend(gotfaces)
        pass

    
    # Create a resulting facearray by creating a mask from allfaces
    masknum=long(0)
    for gotface in allfaces:
        masknum |= (long(1) << gotface.index)
        pass
    
    maskstr=GetMaskStr(masknum)
    result_facearray=facearray.getSequenceFromMask(mask=(maskstr,))
    
    return result_facearray


def GetFace_point_normal(facearray,facepointnormal,pointtolerance,normaltolerance):
    # facepointnormal is a (point,normal) tuples
    # IMPORTANT: Normals must be unit vectors!!!
    
    gotfaces=_get_face_indices_point_normal(facearray,facepointnormal,pointtolerance,normaltolerance)
    if len(gotfaces)==0:
        print("facearray=%s; facepointnormal=%s" % (repr(facearray),str(facepointnormal)))
        pass
    
    # Should be able to find exactly one face from each (point, normal)
    assert(len(gotfaces) > 0) # should be able to find faces
    assert(len(gotfaces) < 2) # should be able to find unique face
    
    
    # Create a resulting facearray by creating a mask from gotfaces
    masknum=long(0)
    masknum |= (long(1) << gotfaces[0].index)
    
    maskstr=GetMaskStr(masknum)
    result_facearray=facearray.getSequenceFromMask(mask=(maskstr,))
    
    return result_facearray



def GetMultipleFaces(facearray,points,pointtolerance):
    # points is either an (x,y,z) tuple or an array of (x,y,z) tuples
    # Attempts to get all faces within the tolerance of any of the given points
    
    gotfaces=[]
    if isinstance(points[0],collections.Sequence):
        for point in points:
            gotfaces.extend(_get_face_indices_point(facearray,point,pointtolerance))
            pass
        pass
    else:
        gotfaces.extend(_get_face_indices_point(facearray,points,pointtolerance))
        pass

    if len(gotfaces)==0:
        print("facearray=%s; points=%s" % (repr(facearray),str(points)))
        pass
    
    # Should be able to find at least one face from each point
    assert(len(gotfaces) > 0) # should be able to find faces
    
    
    # Create a resulting facearray by creating a mask from gotfaces
    masknum=long(0)
    masknum |= (long(1) << gotfaces[0].index)
    
    maskstr=GetMaskStr(masknum)
    result_facearray=facearray.getSequenceFromMask(mask=(maskstr,))
    
    return result_facearray

def GetFace(facearray,points,pointtolerance):
    # points is either an (x,y,z) tuple or an array of (x,y,z) tuples
    # Attempts to get a single face within the tolerance of all of the given points
    
    allfaces=[]
    if not isinstance(points[0],collections.Sequence):
        points=(points,)
        pass
    
    for point in points:
        gotfaces=_get_face_indices_point(facearray,point,pointtolerance)
        if len(gotfaces)==0:
            print("facearray=%s; facepoint=%s" % (repr(facearray),str(point)))
            pass
        assert(len(gotfaces) > 0) # should be able to find faces
        allfaces.extend(gotfaces)
        pass

    
    # Should be able to find at least one face from each point
    assert(len(set([gotface.index for gotface in allfaces])) < 2) # should be exactly one unique face
    
    
    # Create a resulting facearray by creating a mask from gotfaces
    masknum=long(0)
    masknum |= (long(1) << gotfaces[0].index)
    
    maskstr=GetMaskStr(masknum)
    result_facearray=facearray.getSequenceFromMask(mask=(maskstr,))
    
    return result_facearray



def OnEdge(vertexarray,point, edge, pointtolerance):
    
    def bisect_distance(start, end, point):
        try: 
            # Use edge.getCurvature() as a way to find the closest point on the curve
            evalpt = edge.getCurvature(parameter=(end+start)/2.0)['evaluationPoint']
            pass
        except abq.AbaqusException as abqexc:
            if abqexc.message == 'Radius of curvature is not defined for a straight edge':
                (v1, v2)=edge.getVertices()
                p1=vertexarray[v1].pointOn[0] # Corresponds to normalized parameter=0.0
                p2=vertexarray[v2].pointOn[0] # corresponds to normalized parameter=1.0
                evalpt = np.array(p1) + (end+start)/2.0 *(np.array(p2)-np.array(p1))
                pass
            else: 
                raise
            pass
            
        #d = (evalpt[0]-point[0], evalpt[1]-point[1], evalpt[2]-point[2])
        d = np.array(evalpt)-np.array(point) 
        #sys.modules["__main__"].__dict__.update(globals())
        #sys.modules["__main__"].__dict__.update(locals())

        return np.sqrt(d[0]**2.0 + d[1]**2.0 + d[2]**2.0)
    
    edgelen = edge.getSize()
    
    n = 0
    start = 0.0
    end = 1.0
    
    while n <= math.log((pointtolerance/edgelen)**(-1), 2):
        # The maximum error after n iterations is: 1/(2**n)*edgelen
        
        prevstart = start
        prevend = end
        
        d1 = bisect_distance(prevstart, (prevend+prevstart)/2.0, point)
        d2 = bisect_distance((prevend+prevstart)/2.0, prevend, point)
        
        if d1<d2:
            if d1<=pointtolerance:
                return True
                
            start = prevstart
            end = (prevend+prevstart)/2.0
            pass
        else : #  d2<=d1:
            if d2<=pointtolerance:
                return True
                
            start = (prevend+prevstart)/2.0
            end = prevend
            pass
        n+=1
        pass
        
    return False


def _get_edge_indices_point_tangent(edgearray,vertexarray,edgepointtangent,pointtolerance,tangenttolerance):

    point=edgepointtangent[0]
    tangent=edgepointtangent[1]

    gotedges=[]
    for edgecnt in range(len(edgearray)):
        actualtangent=None
        
        # First check if point lies on edge
        if OnEdge(vertexarray,point, edgearray[edgecnt], pointtolerance):
            try:
                actualtangent=edgearray[edgecnt].getCurvature(point=point)['tangent']
                pass
            except abq.AbaqusException as abqexc:
                if abqexc.message == 'Radius of curvature is not defined for a straight edge':
                    (v1, v2)=edgearray[edgecnt].getVertices()
                    p1=vertexarray[v1].pointOn[0]
                    p2=vertexarray[v2].pointOn[0]
                    t=np.array(p2)-np.array(p1)
                    m=np.linalg.norm(t)
                    actualtangent=(t[0]/m,t[1]/m,t[2]/m)
                    pass
                else: 
                    raise
                pass
            pass
        else:
            continue    # Point must be on edge

        if np.abs(np.dot(actualtangent,tangent))-1.0 > -tangenttolerance:
            masknum=(long(1) << edgecnt)
            maskstr=GetMaskStr(masknum)
            
            #if Troubleshoot:
            #    #print("Got tangent",file=sys.stderr)
            #    sys.modules["__main__"].__dict__.update(globals())
            #    sys.modules["__main__"].__dict__.update(locals())
            #    pass
            
            reduced_edgearray=edgearray.getSequenceFromMask(mask=(maskstr,))
            closest=reduced_edgearray.getClosest(coordinates=(point,),searchTolerance=pointtolerance)
            
            assert(len(list(closest.keys()))<=1)  # Should get no more than one result, since there is only one edge involved
            #if Troubleshoot:
            #    raise ValueError("Foo!")
            if len(list(closest.keys())) > 0:
                #if Troubleshoot: 
                #    #print("Got point",file=sys.stderr)
                #    pass
                assert(0 in closest) # only coordinate number is 0
                closestedge=closest[0][0] # first element is edge 
                closestpoint=closest[0][1] # second element is point
                assert(closestedge.index==edgecnt) # should only get the edge we are working on, unless getClosest() is like findAt() and might find the wrong edge 
                assert(np.linalg.norm(np.array(closestpoint)-np.array(point)) <= pointtolerance) # distance should be within tolerance
                gotedges.append(closestedge)
                pass
            pass
        pass
    return gotedges


def _get_edge_indices_point(edgearray,vertexarray,point,pointtolerance):


    gotedges=[]
    for edgecnt in range(len(edgearray)):
        actualtangent=None
        
        # First check if point lies on edge
        if not OnEdge(vertexarray,point, edgearray[edgecnt], pointtolerance):
            continue
            
        masknum=(long(1) << edgecnt)
        maskstr=GetMaskStr(masknum)

        reduced_edgearray=edgearray.getSequenceFromMask(mask=(maskstr,))
        closest=reduced_edgearray.getClosest(coordinates=(point,),searchTolerance=pointtolerance)
            
        assert(len(list(closest.keys()))<=1)  # Should get no more than one result, since there is only one edge involved

        if len(list(closest.keys())) > 0:
            #if Troubleshoot: 
            #    #print("Got point",file=sys.stderr)
            #    pass
            assert(0 in closest) # only coordinate number is 0
            closestedge=closest[0][0] # first element is edge 
            closestpoint=closest[0][1] # second element is point
            assert(closestedge.index==edgecnt) # should only get the edge we are working on, unless getClosest() is like findAt() and might find the wrong edge 
            assert(np.linalg.norm(np.array(closestpoint)-np.array(point)) <= pointtolerance) # distance should be within tolerance
            gotedges.append(closestedge)
            pass
        pass
    return gotedges

def GetEdges_points_tangents(edgearray,vertexarray,edgepointstangents,pointtolerance,tangenttolerance):
    # edgepointtangent is a sequence
    # of (points,tangent) sequences
    # IMPORTANT: Tangents must be unit vectors!!!

    # Troubleshooting:
    #if edgepointtangent[0][1][1]==1.0:
    #    sys.modules["__main__"].__dict__.update(globals())
    #    sys.modules["__main__"].__dict__.update(locals())
    #    Troubleshoot=True
    #    print("Troubleshooting!",file=sys.stderr)
    #    pass
    #else:
    Troubleshoot=False

    alledges=[]
    for pointcnt in range(len(edgepointstangents)):
        gotedges=_get_edge_indices_point_tangent(edgearray,vertexarray,edgepointstangents[pointcnt],pointtolerance,tangenttolerance)
        if len(gotedges)==0:
            print("edgearray=%s; edgepointnormal=%s" % (repr(edgearray),str(edgepointnormal)))
            pass

        # Should be able to find exactly one face from each (point, tangent)
        assert(len(gotedges) > 0) # should be able to find edge
        assert(len(gotedges) < 2) # should be able to find unique edge

        alledges.extend(gotedges)
        pass

    
    # Create a resulting edgearray by creating a mask from alledges
    masknum=long(0)
    for gotedge in alledges:
        masknum |= (long(1) << gotedge.index)
        pass

    maskstr=GetMaskStr(masknum)
    result_edgearray=edgearray.getSequenceFromMask(mask=(maskstr,))

    return result_edgearray


def GetEdge_point_tangent(edgearray,vertexarray,edgepointtangent,pointtolerance,tangenttolerance):
    # edgepointtangent is a (points,tangent) tuple
    # IMPORTANT: Tangents must be unit vectors!!!

    # Troubleshooting:
    #if edgepointtangent[0][1][1]==1.0:
    #    sys.modules["__main__"].__dict__.update(globals())
    #    sys.modules["__main__"].__dict__.update(locals())
    #    Troubleshoot=True
    #    print("Troubleshooting!",file=sys.stderr)
    #    pass
    #else:
    #Troubleshoot=False

    gotedges=_get_edge_indices_point_tangent(edgearray,vertexarray,edgepointtangent,pointtolerance,tangenttolerance)
    if len(gotedges)==0:
        print("edgearray=%s; edgepointnormal=%s" % (repr(edgearray),str(edgepointnormal)))
        pass
    
    # Should be able to find exactly one face from each (point, tangent)
    assert(len(gotedges) > 0) # should be able to find edge
    assert(len(gotedges) < 2) # should be able to find unique edge

    
    # Create a resulting edgearray by creating a mask from alledges
    masknum=long(0)
    masknum |= (long(1) << gotedges[0].index)
    maskstr=GetMaskStr(masknum)
    result_edgearray=edgearray.getSequenceFromMask(mask=(maskstr,))

    return result_edgearray




def GetMultipleEdges(edgearray,vertexarray,points,pointtolerance):
    # points is either an (x,y,z) tuple or an array of (x,y,z) tuples
    # Attempts to get all edges within the tolerance of any of the given points
    
    gotedges=[]
    if isinstance(points[0],collections.Sequence):
        for point in points:
            gotedges.extend(_get_edge_indices_point(edgearray,vertexarray,point,pointtolerance))
            pass
        pass
    else:
        gotedges.extend(_get_edge_indices_point(edgearray,vertexarray,points,pointtolerance))
        pass

    if len(gotedges)==0:
        print("edgearray=%s; points=%s" % (repr(edgearray),str(points)))
        pass
    
    # Should be able to find at least one face from each point
    assert(len(gotedges) > 0) # should be able to find faces
    
    
    # Create a resulting facearray by creating a mask from gotfaces
    masknum=long(0)
    masknum |= (long(1) << gotedges[0].index)
    
    maskstr=GetMaskStr(masknum)
    result_edgearray=edgearray.getSequenceFromMask(mask=(maskstr,))
    
    return result_edgearray



def GetEdge(edgearray,verticesarray,points,pointtolerance):
    # points is either an (x,y,z) tuple or an array of (x,y,z) tuples
    # Attempts to get a single edge within the tolerance of all of the given points
    
    alledges=[]
    if not isinstance(points[0],collections.Sequence):
        points=(points,)
        pass
    
    for point in points:
        gotedges=_get_edge_indices_point(edgearray,verticesarray,point,pointtolerance)
        if len(gotedges)==0:
            print("edgearray=%s; point=%s" % (repr(edgearray),str(point)))
            pass
        assert(len(gotedges) > 0) # should be able to find faces
        alledges.extend(gotedges)
        pass

    
    # Should be able to find at least one edge from each point
    assert(len(set([gotedge.index for gotedge in alledges])) < 2) # should be exactly one unique face
    
    
    # Create a resulting facearray by creating a mask from gotedges
    masknum=long(0)
    masknum |= (long(1) << gotedges[0].index)
    
    maskstr=GetMaskStr(masknum)
    result_edgearray=edgearray.getSequenceFromMask(mask=(maskstr,))
    
    return result_edgearray




def GetVertices(vertexarray,vertexpoints,pointtolerance):
    # vertexpoints is a sequence of points (geometric points)

    gotvertices=[]
    for pointcnt in range(len(vertexpoints)):
        point=vertexpoints[pointcnt]

        for vertexcnt in range(len(vertexarray)):
            d = np.linalg.norm(np.array(vertexarray[vertexcnt].pointOn)-np.array(vertexpoints[pointcnt]))
            #closest=nodearray.getClosest(coordinates=(point,),searchTolerance=pointtolerance)

            if np.abs(d)<pointtolerance:
                gotvertices.append(vertexarray[vertexcnt])
                pass
            pass
        pass

    assert(len(gotvertices) > 0) # should be able to find vertices

    masknum=long(0)
    for gotvertex in gotvertices:
        masknum |= (long(1) << gotvertex.index)
        pass
    
    maskstr=GetMaskStr(masknum)
    result_vertexarray=vertexarray.getSequenceFromMask(mask=(maskstr,))
    
    return result_vertexarray
    

    


def GetEdges_ThreePoints_OBSOLETE(edgearray,vertexarray,edgeendandinteriorpoints,pointtolerance):
    # edgepointtangent is a sequence
    # of (points,tangent) sequences
    # IMPORTANT: Tangents must be unit vectors!!!

    gotedges=[]
    for pointcnt in range(len(edgeandinteriorpoints)):
        endpoint1=edgeandinteriorpoints[pointcnt][0]
        interiorpoint=edgeandinteriorpoinrts[pointcnt][1]
        endpoint2=edgeandinteriorpoints[pointcnt][2]

        for edgecnt in range(len(edgearray)):
            try: 
                (v1, v2)=edgearray[edgecnt].getVertices()
                p1=vertexarray[v1].pointOn[0]
                p2=vertexarray[v2].pointOn[0]
                dist11=np.linalg.norm(np.array(p1)-endpoint1)
                dist12=np.linalg.norm(np.array(p1)-endpoint2)
                dist21=np.linalg.norm(np.array(p2)-endpoint1)
                dist22=np.linalg.norm(np.array(p2)-endpoint2)


                if ((dist11 < pointtolerance and dist22 < pointtolerance) or # endpoints match 1 <-> 1 and 2 <-> 2
                    (dist12 < pointtolerance and dist21 < pointtolerance)): # endpoints match 1 <-> 2
                    
                    masknum=(long(1) << edgecnt)
                    maskstr=GetMaskStr(masknum)
                    reduced_edgearray=edgearray.getSequenceFromMask(mask=(maskstr,))
                    closest=reduced_edgearray.getClosest(coordinates=(interiorpoint,),searchTolerance=pointtolerance)
                    assert(len(list(closest.keys()))<=1)  # Should get no more than one result, since there is only one face involved
                    if len(list(closest.keys())) > 0:                    
                        # got a match
                        assert(0 in closest) # only coordinate number is 0
                        closestedge=closest[0][0] # first element is edge
                        closestpoint=closest[0][1] # second element is point
                        assert(closestedge.index==edgecnt) # should only get the face we are working on, unless getClosest() is like findAt() and might find the wrong face (it doesn't... checked!) 

                        assert(np.linalg.norm(np.array(closestpoint)-np.array(point)) <= pointtolerance) # distance should be within tolerance
                    
                        gotedges.append(closestedge)
                    
                        pass
                    pass
                pass
            
            #except abq.AbaqusException:
            #    pass
            finally:
                pass

            pass

    if len(gotedges) < len(edgeandinteriorpoints):
        print("GetEdges_ThreePoints: Unable to find edges... edgearray=%s; edgeandinteriorpoints=%s" % (repr(edgearray),str(edgeandinteriorpoints)))
        raise ValueError("GetEdges_ThreePoints: Unable to find edges... edgearray=%s; edgeandinteriorpoints=%s" % (repr(edgearray),str(edgeandinteriorpoints)))

    
    assert(len(gotedges) > 0) # should be able to find edges

    # Create a resulting edgearray by creating a mask from gotedges
    masknum=long(0)
    for gotedge in gotedges:
        masknum |= (long(1) << gotedge.index)
        pass
    
    maskstr=GetMaskStr(masknum)
    result_edgearray=edgearray.getSequenceFromMask(mask=(maskstr,))
    
    return result_edgearray


def LookupBody(BodyNumDB,name):
    # Wrapped function to do body number lookup in assembly phase
    return BodyNumDB[name] 



    
