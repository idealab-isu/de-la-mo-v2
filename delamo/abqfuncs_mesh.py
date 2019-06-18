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


def GetNodes(nodearray,nodepoints,pointtolerance):
    # nodepoints is a sequence of points (meshed points)

    gotnodes=[]
    for pointcnt in range(len(nodepoints)):
        point=nodepoints[pointcnt]

        for nodecnt in range(len(nodearray)):
            d = np.linalg.norm(np.array(nodearray[nodecnt].coordinates)-np.array(nodepoints[pointcnt]))
            #closest=nodearray.getClosest(coordinates=(point,),searchTolerance=pointtolerance)

            if np.abs(d)<pointtolerance:
                gotnodes.append(nodearray[nodecnt])
                pass
            pass
        pass

    assert(len(gotnodes) > 0) # should be able to find edges

    masknum=long(0)
    for gotnode in gotnodes:
        masknum |= (long(1) << gotnode.index)
        pass
    
    maskstr=GetMaskStr(masknum)
    result_nodearray=nodearray.getSequenceFromMask(mask=(maskstr,))
    
    return result_nodearray

def SeedPartEdgesByFaces(fe_part,surface_points_and_normals,facepointtolerance,normaltolerance,meshsize):
    # Seed the edges around a set of faces with a particular meshing size
    faces=GetFaces(fe_part.faces,surface_points_and_normals,facepointtolerance,normaltolerance)
    
    edges=()
    for face in faces:
        edges += face.getEdges()
        pass 
    masknum=long(0)
    for edgeidx in edges:
        masknum |= (long(1) << edgeidx)
        pass
    
    maskstr=GetMaskStr(masknum)
    result_edgearray=fe_part.edges.getSequenceFromMask(mask=(maskstr,))
    fe_part.seedEdgeBySize(size=meshsize,edges=result_edgearray)
    pass



def GetLargestFace(part,faces):
    areas = [ part.getArea(faces.getSequenceFromMask(mask=(GetMaskStr(long(1) << facecnt),))) for facecnt in range(len(faces)) ]
    maxidx=np.argmax(areas)
    return faces[maxidx]

    
