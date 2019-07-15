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


class AutoFiber_abq(object):

    __name__ = None
    model = None
    vertices = None
    vertexids = None
    facetnormals = None
    inplanemat = None
    boxes = None
    boxpolys = None
    boxcoords = None

    texcoord2inplane = None

    meshcoords = None
    orientations = None

    def __init__(self, **kwargs):
        for key in kwargs:
            assert (hasattr(self, key))
            setattr(self, key, kwargs[key])
            pass

        pass

    @classmethod
    def CreateFromParams(cls, name, model, fiber):
        return cls(__name__=name, model=model, vertices=fiber[1], vertexids=fiber[2], inplanemat=fiber[3], boxes=fiber[4], boxpolys=fiber[5], boxcoords=fiber[6], facetnormals=fiber[7], texcoord2inplane=fiber[0])

    def getMeshCenters(self):
        self.meshcoords = np.empty((0, 3))
        elements = self.model.parts[self.__name__].elements
        for i in range(0, len(elements)):
            centroid = np.empty((0, 3))
            for node in elements[i].getNodes():
                centroid = np.vstack((centroid, node.coordinates))

            centroid = np.sum(centroid, axis=0) / centroid.shape[0]

            self.meshcoords = np.vstack((self.meshcoords, centroid))
        pass

    def getFiberOrientations(self):
        if self.meshcoords is not None:
            elements = []
            directions = []

            rangex = np.array([self.boxcoords[:, 0], self.boxcoords[:, 3]]).T
            rangey = np.array([self.boxcoords[:, 1], self.boxcoords[:, 4]]).T
            rangez = np.array([self.boxcoords[:, 2], self.boxcoords[:, 5]]).T

            def in_box(p):
                containers = np.where(np.logical_and(np.logical_and(
                    np.logical_and(p[0] > rangex[:, 0], p[0] < rangex[:, 1]),
                    np.logical_and(p[1] > rangey[:, 0], p[1] < rangey[:, 1])),
                    np.logical_and(p[2] > rangez[:, 0], p[2] < rangez[:, 1])))[0]
                bindx = self.boxes[containers][:, -1][self.boxes[containers][:, -1] != -1][0]

                polyslist = np.array([self.boxpolys[bindx]])
                count = 1
                while True:
                    if self.boxpolys[bindx + count] == -1:
                        break
                    polyslist = np.append(polyslist, self.boxpolys[bindx + count])
                    count += 1

                return bindx, polyslist

            for i in range(0, self.meshcoords.shape[0]):
                idx, element = None, None

                point = self.meshcoords[i]

                try:
                    box_idx, polys = in_box(point)
                except IndexError:
                    distances = np.sqrt(np.einsum('ij,ij->i', self.vertices - point, self.vertices - point))
                    idx = np.where(distances == np.min(distances))
                    vert = self.vertices[idx][0]
                    box_idx, polys = in_box(vert)

                for j in polys:
                    check, projpnt = self.check_proj_inplane_pnt(point, self.vertices[self.vertexids[j]])
                    if self.point_in_polygon_3d(self.vertices[self.vertexids][j], projpnt, self.inplanemat[j]):
                        element = j
                        break

                if element is None:
                    if idx is None:
                        distances = np.sqrt(np.einsum('ij,ij->i', self.vertices - point, self.vertices - point))
                        idx = np.where(distances == np.min(distances))
                    vertneighbors = np.unique(np.where((self.vertexids == idx))[0])
                    element = vertneighbors[0]

                if element is not None:
                    red_texcoords2inplane = self.texcoord2inplane[element][:2, :2]
                    texutexbasis = np.array([1.0, 0.0])
                    texu2dbasis = np.dot(red_texcoords2inplane, texutexbasis)
                    u3D = np.dot(self.inplanemat[element].T, texu2dbasis)
                    elements.append(i+1)
                    directions += self.calcunitvector(u3D).tolist() + np.cross(self.facetnormals[element], self.calcunitvector(u3D)).tolist()
                else:
                    print("Failed to find point on surface: %s" % point)
            orientations = (('', 6, tuple(elements), tuple(directions)),)
            self.orientations = orientations
        else:
            print("No orientation points defined.")

    def CreateDiscreteField(self):
        self.model.DiscreteField(name='%s_orientation' % self.__name__, location=abqC.ELEMENTS, fieldType=abqC.ORIENTATION,
                                 dataWidth=6, defaultValues=(0, 1, 0, 1, 0, 0),
                                 data=self.orientations,
                                 orientationType=abqC.CARTESIAN, partLevelOrientation=True)

    def calcunitvector(self, vector):
        """ Returns the unit vector of the vector.  """
        if len(vector.shape) >= 2:
            return vector / np.linalg.norm(vector, axis=1)[:, np.newaxis]
        else:
            return vector / np.linalg.norm(vector)

    def check_proj_inplane_pnt(self, point, element_vertices):
        """
        https://math.stackexchange.com/questions/544946/determine-if-projection-of-3d-point-onto-plane-is-within-a-triangle
        :param point: Test point to check
        :param element_vertices: Vertices of current element
        :return: True or False, depending on if the projected point is inside or outside
        """
        u = element_vertices[1] - element_vertices[0]
        v = element_vertices[2] - element_vertices[0]
        normal = np.cross(u, v)
        w = point - element_vertices[0]

        gamma = np.dot(np.cross(u, w), normal) / np.dot(normal, normal)
        beta = np.dot(np.cross(w, v), normal) / np.dot(normal, normal)
        alpha = 1 - gamma - beta

        projpoint = alpha * element_vertices[0] + beta * element_vertices[1] + gamma * element_vertices[2]
        if 0 <= alpha <= 1 and 0 <= beta <= 1 and 0 <= gamma <= 1:
            return True, projpoint
        else:
            return False, projpoint

    def point_in_polygon_2d(self, vertices_rel_point):
        import sys
        # Apply winding number algorithm.
        # This algorithm is selected -- in its most simple form --
        # because it is so  simple and robust in the case of the
        # intersect point being on or near the edge. It may well
        # be much slower than optimal. It tries to return True
        # in the edge case.

        # Should probably implement a faster algorithm then drop
        # down to this for the special cases.

        # See Hormann and Agathos, The point in polygon problem
        # for arbitrary polygons, Computational Geometry 20(3) 131-144 (2001)
        # http://dx.doi.org/10.1016/S0925-7721(01)00012-8
        # https://pdfs.semanticscholar.org/e90b/d8865ddb7c7af2b159d413115050d8e5d297.pdf

        # Winding number is sum over segments of
        # acos((point_to_vertex1 dot point_to_vertex2)/(magn(point_to_vertex1)*magn(point_to_vertex_2))) * sign(det([ point_to_vertex1  point_to_vertex2 ]))
        # where sign(det) is really: What is the sign of the z
        # component of (point_to_vertex1 cross point_to_vertex2)

        # Special cases: magn(point_to_vertex1)==0 or
        #  magn_point_to_vertex2   -> point is on edge
        # det([ point_to_vertex1  point_to_vertex2 ]) = 0 -> point may be on edge

        windingnum = 0.0
        numvertices = vertices_rel_point.shape[0]

        for VertexCnt in range(numvertices):
            NextVertex = VertexCnt + 1
            if NextVertex == numvertices:
                # final vertex... loop back to the start
                NextVertex = 0
                pass

            # calculate (thisvertex - ourpoint) -> vec1
            vec1 = vertices_rel_point[VertexCnt, :]
            magn1 = np.linalg.norm(vec1)

            # calculate (nextvertex - ourpoint) -> vec2
            vec2 = vertices_rel_point[NextVertex, :]
            magn2 = np.linalg.norm(vec2)

            if magn1 == 0.0 or magn2 == 0.0:
                # Got it!!!
                return True

            vec1 = vec1 / magn1
            vec2 = vec2 / magn2

            det = vec1[0] * vec2[1] - vec2[0] * vec1[1]  # matrix determinant

            cosparam = (vec1[0] * vec2[0] + vec1[1] * vec2[1])  # /(magn1*magn2);

            if cosparam < -1.0:
                # Shouldn't be possible...just in case of weird roundoff
                cosparam = -1.0

            if cosparam > 1.0:
                # Shouldn't be possible...just in case of weird roundoff
                cosparam = 1.0

            if det > 0:
                windingnum += np.arccos(cosparam)
            elif det < 0:
                windingnum -= np.arccos(cosparam)
            else:
                # det==0.0
                # Vectors parallel or anti-parallel

                if cosparam > 0.9:
                    # Vectors parallel. We are OUTSIDE. Do Nothing
                    pass
                elif cosparam < -0.9:
                    # Vectors anti-parallel. We are ON EDGE */
                    return True
                else:
                    assert 0  # Should only be able to get cosparam = +/- 1.0 if abs(det) > 0.0 */
                    pass
                pass
            pass

        windingnum = abs(windingnum) * (
                    1.0 / (2.0 * np.pi))  # divide out radians to number of winds; don't care about clockwise vs. ccw
        if windingnum > .999 and windingnum < 1.001:
            # Almost exactly one loop... got it!
            return True
        elif windingnum >= .001:
            #
            sys.stderr.write(
                "spatialnde.geometry.point_in_polygon_2d() Got weird winding number of %e; assuming inaccurate calculation on polygon edge\n" % (
                    windingnum))
            # Could also be self intersecting polygon
            # got it !!!
            return True

        # If we got this far, the search failed
        return False

    def point_in_polygon_3d(self, vertices, point, inplanemat):
        """ assumes vertices are coplanar, with given orthonormal 2D basis inplanemat.  """
        vert3d_rel_point = vertices-point[np.newaxis, :]
        vert2d_rel_point = np.inner(vert3d_rel_point, inplanemat)

        return self.point_in_polygon_2d(vert2d_rel_point)
