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


import numpy as np
np.seterr(divide='ignore', invalid='ignore')
from spatialnde import geometry


class EdgeError(Exception):
    pass


def calcunitvector(vector):
    """ Returns the unit vector of the vector.  """
    if len(vector.shape) >= 2:
        return vector / np.linalg.norm(vector, axis=1)[:, np.newaxis]
    else:
        if np.linalg.norm(vector) > 0.0:
            return vector / np.linalg.norm(vector)
        else:
            return vector


def calcnormal(points):
    """ Returns the normal for the given 2d points"""
    v1 = points[2] - points[0]
    v2 = points[1] - points[0]
    return np.cross(v1, v2)


def angle_between_vectors(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'
        https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
    """
    v1_u = calcunitvector(v1)
    v2_u = calcunitvector(v2)
    detarray = np.vstack((v1_u, v2_u))
    detarray = np.vstack((detarray, np.ones(detarray.shape[1]))).T
    det = np.linalg.det(detarray)
    if det < 0:
        return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
    else:
        return 2*np.pi - np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


def vector_inbetween(v1, v2, v3, error=1e-10):
    """
    Determines if a vector (v1) is between v2 and v3
    https://stackoverflow.com/questions/13640931/how-to-determine-if-a-vector-is-between-two-other-vectors
    :param v1: Test vector
    :param v2: Given vector
    :param v3: Given vector
    :return: True if vector is between v2 and v3, false if not between
    """
    if np.dot(np.cross(v2, v1), np.cross(v2, v3)) >= -error and np.dot(np.cross(v3, v1), np.cross(v3, v2)) >= -error:
        return True
    else:
        return False


def rot_vector(oldnormal, newnormal, vector, force=False):
    """ Rotate a vector given an axis and an angle of rotation
        Returns: Vector reoriented from an old element face to a new element
        https://en.wikipedia.org/wiki/Rodrigues'_rotation_formula
    """
    if np.linalg.norm(oldnormal - newnormal) < 1e-10:
        return vector
    else:
        vector_a = np.cross(oldnormal, newnormal)
        sinphi = np.linalg.norm(vector_a)
        cosphi = np.dot(oldnormal, newnormal)
        a_hat = vector_a/sinphi
        if np.arccos(cosphi) >= np.deg2rad(85.0) and not force:
            # print("Edge detected...geodesic path completed")
            raise EdgeError
        else:
            return calcunitvector(vector * cosphi - np.cross(vector, a_hat) * sinphi + a_hat * np.dot(vector, a_hat) * (1 - cosphi))


def rot_vector_angle(vector, normal, angle):
    return calcunitvector(vector * np.cos(np.deg2rad(angle)) + np.cross(normal, vector) * np.sin(np.deg2rad(angle)) + normal * np.dot(normal, vector) * (1 - np.cos(np.deg2rad(angle))))


def check_proj_inplane_pnt(point, element_vertices):
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


def proj_vector(vector, newnormal):
    """
    Project a vector onto a surface defined by newnormal
    :param vector: Vector to be projected
    :param newnormal: Normal of projected surface
    :return: Vector projected on surface defined by newnormal
    """
    return calcunitvector(vector - np.dot(vector, newnormal) * newnormal)


def check_intersection(p1, q1, p2, q2):
    """
    Check for an intersection between (p1, p2) and (q1, q2)
    https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
    """
    def orientation(p, q, r):
        val = (q[1] - p[1]) * (r[0] - q[0]) - (q[0] - p[0])*(r[1] - q[1])
        if val == 0.0:
            return 0
        elif val > 0:
            return 1
        else:
            return 2

    def onSegment(p, q, r):
        if q[0] <= np.max([p[0], r[0]]) and q[0] >= np.min([p[0], r[0]]) and q[1] <= np.max([p[1], r[1]]) and q[1] >= np.min([p[1], r[1]]):
            return True
        return False

    if (p1 == p2).all():
        return False

    o1 = orientation(p1, q1, p2)
    o2 = orientation(p1, q1, q2)
    o3 = orientation(p2, q2, p1)
    o4 = orientation(p2, q2, q1)

    if o1 != o2 and o3 != o4:
        return True

    # Special Cases
    # p1, q1 and p2 are colinear and p2 lies on segment p1q1
    if o1 == 0 and onSegment(p1, p2, q1):
        return True
    # p1, q1 and q2 are colinear and q2 lies on segment p1q1
    if o2 == 0 and onSegment(p1, q2, q1):
        return True
    # p2, q2 and p1 are colinear and p1 lies on segment p2q2
    if o3 == 0 and onSegment(p2, p1, q2):
        return True
    # p2, q2 and q1 are colinear and q1 lies on segment p2q2
    if o4 == 0 and onSegment(p2, q1, q2):
        return True

    # Doesn't fall in any of the above cases
    return False


def find_intpnt(P1, P2, P3, P4):
    """ Line-Line intersection method
        Returns: A point in 2d that intersects line P1P2 and P3P4
        https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection
    """
    return np.array([((P1[0] * P2[1] - P1[1] * P2[0]) * (P3[0] - P4[0]) - (P1[0] - P2[0]) * (P3[0] * P4[1] - P3[1] * P4[0])) /
                     ((P1[0] - P2[0]) * (P3[1] - P4[1]) - (P1[1] - P2[1]) * (P3[0] - P4[0])),
                     ((P1[0] * P2[1] - P1[1] * P2[0]) * (P3[1] - P4[1]) - (P1[1] - P2[1]) * (P3[0] * P4[1] - P3[1] * P4[0])) /
                     ((P1[0] - P2[0]) * (P3[1] - P4[1]) - (P1[1] - P2[1]) * (P3[0] - P4[0]))])


def find_edge(point, direction, bary, error):
    """
    Determines which edge number is intersected first (0, 1, 2) -> (d12, d23, d31)
    https://math.stackexchange.com/questions/2292895/walking-on-the-surface-of-a-triangular-mesh
    :param point: Start point
    :param direction: Current fiber direction
    :param error: Numerical tolerance
    :return: Edge number (0, 1, 2) or -1 if on an edge
    """
    point = point + calcunitvector(np.sum(bary, axis=0) / 3 - point) * 1e-8

    if direction[1] != 0.0:
        d0 = -point[1] / direction[1]
    else:
        d0 = -1

    if direction[0] + direction[1] != 0.0:
        d1 = (1 - point[0] - point[1]) / (direction[0] + direction[1])
    else:
        d1 = -1

    if direction[0] != 0.0:
        d2 = -point[0] / direction[0]
    else:
        d2 = -1

    if d0 >= error and (d0 <= d1 or d1 <= error) and (d0 <= d2 or d2 <= error):
        return 0
    elif d1 >= error and (d1 <= d0 or d0 <= error) and (d1 <= d2 or d2 <= error):
        return 1
    elif d2 >= error and (d2 <= d0 or d0 <= error) and (d2 <= d1 or d1 <= error):
        return 2
    elif d0 == -1 or d1 == -1 or d2 == -1:
        # print("Following an edge...")
        raise EdgeError
    elif d0 < 0 and d1 < 0 and d2 < 0:
        return find_edge(point, -direction, bary, error)
    else:
        # print("Edge is uncertain...terminating geodesic")
        raise EdgeError


def find_neighbors(element, vertexids_indices, adjacencyidx):
    """
    Finds neighboring elements
    :param element: Current element
    :param vertexids_indices: Indices of the mesh indices
    :param adjacencyidx: Built from spatialnde, index of element adjacency
    :return: An array of element numbers that neighbor the current element
    """
    firstidx = vertexids_indices[element]
    neighbors = np.array([element, adjacencyidx[firstidx]])
    indx = 1
    counter = 1
    while indx > 0:
        indx = adjacencyidx[firstidx + counter]
        if np.isscalar(indx):
            if indx != -1:
                neighbors = np.append(neighbors, indx)
            counter += 1
        else:
            print("### Edge encountered ###")
            break
    return neighbors


def calcdistance(unitvector, oldvertex, meshpoints):
    """
    Calculate perpendicular distance between a ray and a point
    :param unitvector: Reference vector to calculate distance from
    :param oldvertex: Start point for unitvector
    :param meshpoints: Test points
    :return: Perpendicular and parallel distance to each mesh point
    """
    perpvectors = -1*((oldvertex - meshpoints) - np.multiply(np.dot((oldvertex - meshpoints), unitvector[:, np.newaxis]), unitvector[np.newaxis, :]))
    paraldistance = np.dot(meshpoints - oldvertex, unitvector)
    return perpvectors, paraldistance


def calcclosestpoint(unitvector, oldpoint, meshpoints, normal):
    """
    Find closest mesh vertex defined by the distances calculated in calcdistance
    :param unitvector: Reference direction vector
    :param oldpoint: Start point for unitvector
    :param meshpoints: All test points
    :return: Closest point relative to unitvector
    """
    trimedmeshpoints = np.delete(meshpoints, np.where((meshpoints == oldpoint).all(axis=1)), axis=0)
    perpdistances, paraldistances = calcdistance(unitvector, oldpoint, trimedmeshpoints)
    point_idx = np.argmin(np.linalg.norm(perpdistances, axis=1))

    fpointu = paraldistances[point_idx]
    point_3d = trimedmeshpoints[point_idx]
    vector2pnt = perpdistances[point_idx]
    testval = np.dot(calcunitvector(np.cross(unitvector, vector2pnt)), normal)

    fpointv = testval * np.linalg.norm(vector2pnt)
    if np.isnan(fpointv):
        fpointv = 0.0
    fpoint = np.array([fpointu, fpointv])
    return fpoint, point_3d


def calcbarycentric(point, element_vertices):
    """
    Convert 3d point to barycenteric coordinates
    https://en.wikipedia.org/wiki/Barycentric_coordinate_system
    https://math.stackexchange.com/questions/2292895/walking-on-the-surface-of-a-triangular-mesh
    :param point: 3d point
    :param element_vertices: Vertices of current element
    :return: point in barycenteric coordinates
    """
    d1 = (element_vertices[0, 0] * (element_vertices[1, 1] - element_vertices[2, 1]) + element_vertices[1, 0] * (element_vertices[2, 1] - element_vertices[0, 1]) + element_vertices[2, 0] * (element_vertices[0, 1] - element_vertices[1, 1]))
    d2 = (element_vertices[0, 0] * (element_vertices[1, 2] - element_vertices[2, 2]) + element_vertices[1, 0] * (element_vertices[2, 2] - element_vertices[0, 2]) + element_vertices[2, 0] * (element_vertices[0, 2] - element_vertices[1, 2]))
    d3 = (element_vertices[0, 1] * (element_vertices[1, 2] - element_vertices[2, 2]) + element_vertices[1, 1] * (element_vertices[2, 2] - element_vertices[0, 2]) + element_vertices[2, 1] * (element_vertices[0, 2] - element_vertices[1, 2]))

    a1 = np.abs(d1)
    a2 = np.abs(d2)
    a3 = np.abs(d3)

    if a1 >= a2 and a1 >= a3:
        u = (point[0] * (element_vertices[2, 1] - element_vertices[0, 1]) + element_vertices[0, 0] * (
            point[1] - element_vertices[2, 1]) + element_vertices[2, 0] * (element_vertices[0, 1] - point[1])) / d1
        v = (point[0] * (element_vertices[0, 1] - element_vertices[1, 1]) + element_vertices[0, 0] * (
            element_vertices[1, 1] - point[1]) + element_vertices[1, 0] * (point[1] - element_vertices[0, 1])) / d1
    elif a2 >= a1 and a2 >= a3:
        u = (point[0] * (element_vertices[2, 2] - element_vertices[0, 2]) + element_vertices[0, 0] * (
            point[2] - element_vertices[2, 2]) + element_vertices[2, 0] * (element_vertices[0, 2] - point[2])) / d2
        v = (point[0] * (element_vertices[0, 2] - element_vertices[1, 2]) + element_vertices[0, 0] * (
            element_vertices[1, 2] - point[2]) + element_vertices[1, 0] * (point[2] - element_vertices[0, 2])) / d2
    else:
        u = (point[1] * (element_vertices[2, 2] - element_vertices[0, 2]) + element_vertices[0, 1] * (
            point[2] - element_vertices[2, 2]) + element_vertices[2, 1] * (element_vertices[0, 2] - point[2])) / d3
        v = (point[1] * (element_vertices[0, 2] - element_vertices[1, 2]) + element_vertices[0, 1] * (
            element_vertices[1, 2] - point[2]) + element_vertices[1, 1] * (point[2] - element_vertices[0, 2])) / d3
    return u, v


def invcalcbarycentric(pointuv, element_vertices):
    """
    Convert barycenteric coordinates into 3d
    https://en.wikipedia.org/wiki/Barycentric_coordinate_system
    https://math.stackexchange.com/questions/2292895/walking-on-the-surface-of-a-triangular-mesh
    :param pointuv: Point in barycenteric coordinates (u, v)
    :param element_vertices: Vertices of current element
    :return: pointuv in 3d coordinates (x, y, z)
    """
    return element_vertices[0] + pointuv[0] * (element_vertices[1] - element_vertices[0]) + pointuv[1] * (element_vertices[2] - element_vertices[0])


def calcbarycentricdirection(vector, element_vertices):
    """
    Convert a direction vector from 3d to barycenteric coordinates
    https://en.wikipedia.org/wiki/Barycentric_coordinate_system
    https://math.stackexchange.com/questions/2292895/walking-on-the-surface-of-a-triangular-mesh
    :param vector: Direction vector in 3d
    :param element_vertices: Vertices of current element
    :return: Vector in barycenteric coordinates (du, dv)
    """
    d1 = (element_vertices[0, 0] * (element_vertices[1, 1] - element_vertices[2, 1]) + element_vertices[1, 0] * (element_vertices[2, 1] - element_vertices[0, 1]) + element_vertices[2, 0] * (element_vertices[0, 1] - element_vertices[1, 1]))
    d2 = (element_vertices[0, 0] * (element_vertices[1, 2] - element_vertices[2, 2]) + element_vertices[1, 0] * (element_vertices[2, 2] - element_vertices[0, 2]) + element_vertices[2, 0] * (element_vertices[0, 2] - element_vertices[1, 2]))
    d3 = (element_vertices[0, 1] * (element_vertices[1, 2] - element_vertices[2, 2]) + element_vertices[1, 1] * (element_vertices[2, 2] - element_vertices[0, 2]) + element_vertices[2, 1] * (element_vertices[0, 2] - element_vertices[1, 2]))

    a1 = np.abs(d1)
    a2 = np.abs(d2)
    a3 = np.abs(d3)

    if a1 >= a2 and a1 >= a3:
        du = (vector[0] * (element_vertices[2, 1] - element_vertices[0, 1]) + vector[1] * (
        element_vertices[0, 0] - element_vertices[2, 0])) / d1
        dv = (vector[0] * (element_vertices[0, 1] - element_vertices[1, 1]) + vector[1] * (
        element_vertices[1, 0] - element_vertices[0, 0])) / d1
    elif a2 >= a1 and a2 >= a3:
        du = (vector[0] * (element_vertices[2, 2] - element_vertices[0, 2]) + vector[2] * (
        element_vertices[0, 0] - element_vertices[2, 0])) / d2
        dv = (vector[0] * (element_vertices[0, 2] - element_vertices[1, 2]) + vector[2] * (
        element_vertices[1, 0] - element_vertices[0, 0])) / d2
    else:
        du = (vector[1] * (element_vertices[2, 2] - element_vertices[0, 2]) + vector[2] * (
        element_vertices[0, 1] - element_vertices[2, 1])) / d3
        dv = (vector[1] * (element_vertices[0, 2] - element_vertices[1, 2]) + vector[2] * (
        element_vertices[1, 1] - element_vertices[0, 1])) / d3
    return du, dv


def invcalcbarycentricdirection(vectoruv, element_vertices):
    """
    Convert vector in barycenteric coordinates into a 3d vector
    https://en.wikipedia.org/wiki/Barycentric_coordinate_system
    https://math.stackexchange.com/questions/2292895/walking-on-the-surface-of-a-triangular-mesh
    :param vectoruv: Vector in barycenteric coordinate (du, dv)
    :param element_vertices: Vertices of current element
    :return: Vectoruv in 3d space (dx, dy, dz)
    """
    # element_vertices = vertices[vertexids[element]]
    return vectoruv[0] * (element_vertices[1] - element_vertices[0]) + vectoruv[1] * (element_vertices[2] - element_vertices[0])


def check_inplane_pnt(point, element_vertices):
    """
    Determines if a point is within the plane of the current element face
    :param point: A point within or on the edge of the current element
    :param element_vertices: Vertices of current element
    :return: True if the point is within in the plane, or False if otherwise
    """
    detarray = np.vstack((element_vertices, point)).T
    detarray = np.vstack((detarray, np.ones(detarray.shape[1])))
    det = np.linalg.det(detarray)
    err = 1e-10
    if det == 0 or np.abs(det) < err:
        return True
    else:
        return False


def check_inplane_vector(vector, normal):
    """
    Determines if a vector is in plane with the current element
    :param vector: Test vector
    :param normal: Normal of element
    :return:
    """
    checkdot = np.dot(normal, vector)
    err = 1e-9
    if -err < checkdot < err:
        return True
    else:
        return False


def find_element_vertex(point, unitvector, curnormal, vertices, vertexids, facetnormals):
    """
    Determines which element is next given a vertex and an angle
    :param point: Vertex in the mesh
    :param unitvector: Fiber direction vector
    :param curnormal: Current element normal direction vector
    :param vertices: Mesh vertices
    :param vertexids: Id's of mesh element vertices
    :param facetnormals: Normals of each element in mesh
    :return: The element in which the fiber direction vector resides
    """
    newvector = None
    element = None
    # Find neighboring elements:
    idx = np.where(np.linalg.norm(vertices - point, axis=1) == np.min(np.linalg.norm(vertices - point, axis=1)))
    neighbors = np.unique(np.where((vertexids == idx))[0])

    # For each neighbor and the given direction determine which will include the next point
    nvectors = vertices[vertexids[neighbors]] - point
    for i in neighbors:
        # Calculate fiber direction at current point
        try:
            newvector = rot_vector(curnormal, facetnormals[i], unitvector)
        except EdgeError:
            continue
        if check_inplane_vector(newvector, facetnormals[i]):
            nvector = nvectors[np.where((neighbors == i))]
            evectors = calcunitvector(nvector[0][np.unique(np.nonzero(nvector)[1])].reshape(-1, 3))
            if vector_inbetween(newvector, evectors[0], evectors[1]):
                element = i
                break
    return element, newvector


def find_element_within(point, unitvector, normal, vertices, vertexids, facetnormals, inplanemat):
    """
    Determines which element a point is within
    :param point: Vertex in the mesh
    :param unitvector: Fiber direction vector
    :param normal: Current element normal direction vector
    :return: The element that the point is within
    """
    vertexid = np.where(np.linalg.norm(vertices - point, axis=1) == np.min(np.linalg.norm(vertices - point, axis=1)))
    neighbors = np.where(vertexids == vertexid)[0]

    for i in neighbors:
        if check_inplane_pnt(point, vertices[vertexids[i, :]]):
            if geometry.point_in_polygon_3d(vertices[vertexids[i]], point, inplanemat[i]):
                if check_inplane_vector(unitvector, facetnormals[i]):
                    return i, unitvector, None
                else:
                    try:
                        newvector = rot_vector(normal, facetnormals[i], unitvector)
                    except EdgeError:
                        continue
                    if check_inplane_vector(newvector, facetnormals[i]):
                        return i, newvector, None
        else:
            test, projpnt = check_proj_inplane_pnt(point, vertices[vertexids[i]])
            if test:
                if geometry.point_in_polygon_3d(vertices[vertexids[i]], projpnt, inplanemat[i]):
                    if check_inplane_vector(unitvector, facetnormals[i]):
                        return i, unitvector, projpnt
                    else:
                        try:
                            newvector = rot_vector(normal, facetnormals[i], unitvector)
                        except EdgeError:
                            continue
                        if check_inplane_vector(newvector, facetnormals[i]):
                            return i, newvector, projpnt
    return None, None, None


def traverse_element(af, element, point, unitfiberdirection, length, uv_start, direction=1, parameterization=True):
    """
    Traverse a triangular element
    :param af: Autofiber object
    :param element: Current triangular element
    :param point: Current point
    :param unitfiberdirection: Current direction vector
    :param length: Current length of the geodesic
    :param uv_start: Start point of geodesic in uv space
    :param direction: (1) for positive geodesic direction, (-1) for negative geodesic direction
    :param parameterization: Are we going to record geodesic details for use in parameterization calculation?
    :return: next intersection point, next unitfiberdirection based on next element, next element
    """
    # Determine the elements surrounding the current element
    neighbors = find_neighbors(element, af.vertexids_indices, af.adjacencyidx)

    # Retrieve current element vertices
    element_vertices = af.vertices[af.vertexids[element]]

    # Calculate the barycentric coordinates for each vertex of the current element
    element_vertices_bary = np.zeros((element_vertices.shape[0], 2))
    for i in range(0, element_vertices.shape[0]):
        element_vertices_bary[i, :] = calcbarycentric(element_vertices[i], element_vertices)

    # Calculate the barycentric coordinates for the current 3d point
    pointuv = np.array(calcbarycentric(point, element_vertices))

    # Calculate the barycentric direction for the current 3d fiber direction vector
    duv = calcunitvector(np.array(calcbarycentricdirection(unitfiberdirection, element_vertices)))

    # Calculate another point in the direction of the fiber in order to calculate intersection point later
    lnpoint = pointuv + duv * 1.1

    # Signed distance to each edge ([d12, d23, d31])
    edges_dict = np.array([[0, 1], [1, 2], [2, 0]])
    # Determine which edge will be intersected first
    edge_num = find_edge(pointuv, duv, element_vertices_bary, 0.0000000005)
    # Retrieve the corresponding vertex indices to the intersected edge
    edge = edges_dict[edge_num]
    nextedge = element_vertices_bary[edge]

    # Find the point of intersection between the edge and a line in the fiber direction
    int_pnt = find_intpnt(pointuv, lnpoint, nextedge[0], nextedge[1])

    intersected = []
    prev_lines = af.georecord.get(element, [[], None])[0]
    for line in prev_lines:
        if check_intersection(pointuv, int_pnt, line[0], line[1]):
            # print("Intersection detected...terminating current geodesic")
            int_pnt_3d = invcalcbarycentric(find_intpnt(pointuv, int_pnt, line[0], line[1]), element_vertices)
            distance = np.linalg.norm(int_pnt_3d - point)
            intersected.append((distance, int_pnt_3d))
    if len(intersected) > 0:
        return min(intersected, key=lambda t: t[0])[1], None, None

    if parameterization:
        af.fiberdirections[element] = direction * unitfiberdirection
        if element not in list(af.georecord.keys()):
            af.georecord[element] = [[], None]

        fpoint, closest_point = calcclosestpoint(direction * unitfiberdirection, point, element_vertices, af.facetnormals[element])
        closest_point_idx = np.where((af.vertices == closest_point).all(axis=1))[0][0]
        af.georecord[element][0].append(
            (pointuv, int_pnt, point, unitfiberdirection, closest_point_idx, uv_start, length, direction))

    # Retrieve the 3d coordinates of the edge vertices
    nextedgec = af.vertices[af.vertexids[element, edge]]
    # Remove the current element from the neighbors array
    neighbors = np.delete(neighbors, np.where((neighbors == element)))
    # Redefine the all_vertices variable to not include the current element
    all_vertices = af.vertices[af.vertexids[neighbors]]
    # Find which element is across the intersected edge
    nextelement = None
    for i in range(0, all_vertices.shape[0]):
        # Check for the first vertex of the edge
        test = np.linalg.norm(all_vertices[i] - nextedgec[0], axis=1)
        if 0 in test:
            # If the first vertex is in this element check for the second vertex of the edge
            test2 = np.linalg.norm(all_vertices[i] - nextedgec[1], axis=1)
            if 0 in test2:
                # If both vertices are in this element set nextelement to the current element
                nextelement = neighbors[i]
                break
    # Convert the intersection point into 3d coordinates
    int_pnt_3d = invcalcbarycentric(int_pnt, element_vertices)

    if nextelement is None:
        return int_pnt_3d, None, nextelement

    # Rotate current 3d fiber vector to match the plane of the next element
    nextunitvector = rot_vector(af.facetnormals[element], af.facetnormals[nextelement], unitfiberdirection)

    return int_pnt_3d, nextunitvector, nextelement
