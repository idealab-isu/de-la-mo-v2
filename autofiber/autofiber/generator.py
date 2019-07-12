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


import os, sys, time
import numpy as np
np.seterr(divide='ignore', invalid='ignore')

from spatialnde.coordframes import coordframe
from spatialnde.ndeobj import ndepart
from spatialnde.cadpart.polygonalsurface_texcoordparameterization import polygonalsurface_texcoordparameterization

from . import geodesic as GEO
from . import analyze_uv as AUV
from . import optimization as OP


def calcunitvector(vector):
    """ Returns the unit vector of the given vector.  """
    if len(vector.shape) >= 2:
        return vector / np.linalg.norm(vector, axis=1)[:, np.newaxis]
    else:
        if np.linalg.norm(vector) > 0.0:
            return vector / np.linalg.norm(vector)
        else:
            return vector


class AutoFiber:

    def __init__(self, cadfile, initpoint, initdirection, initnormal, materialproperties=(228.0, 0.2, None), fiberint=0.1, angle_error=0.01, accel=False):
        """
        Calculate geodesic based parameterization of a triangular meshed CAD model

        :param cadfile: Path to CAD file (currently supports x3d, stl, and De-La-Mo DMObjects)
        :param initpoint: 3D point closest to the center of the surface we would like to work on
        :param initdirection: 3D unit vector representing the fiber direction at initpoint
        :param initnormal: 3D unit vector indicating the surface normal at initpoint (used to determine which surface to operate on)
        :param materialproperties: Composite fiber material properties, if empty then a set of default properties will be used.
        materialproperties is setup as follows: (E, poisson's ratio, G) if E is a list: for anisotropic [E1, E2, E3],
        nu is [nu12, nu13, nu23], and the shear modulus G is [G12, G13, G23]. If E is not a list then an isotropic material
        is used and G will be computed from E and nu.
        :param fiberint: Perpendicular distance between generated geodesics
        :param angle_error: Error incorporated into the initdirection, any error defined here is reversed during optimization
        :param accel: Utilize OpenCL parallel geodesic generator (WIP - not functioning)
        """
        # Get CAD file of part that we would like to parameterize
        self.cadfile = cadfile
        # Get the angle error which will be applied to initdirection
        self.error = angle_error

        # Gather options, set to default if option doesn't exist
        # accel: activate opencl optimization features (WIP)
        self.accel = accel

        # Init spatialnde objects
        self.obj = None
        self.objframe = None

        # load model into spatialnde obj
        self.loadobj()

        # Init model variables
        self.vertices = None
        self.vertexids_indices = None
        self.vertexids = None
        self.facetnormals = None
        self.refpoints = None
        self.inplanemat = None
        self.edges = None
        self.adjacencyidx = None
        self.surfaces = None
        self.surface_vertexids = None
        self.boxes = None
        self.boxpolys = None
        self.boxcoords = None

        # load relevant model variables
        self.loadvars()

        # Init fiber material properties
        # Defaults are just basic isotropic material properties
        self.fiberint = fiberint
        self.E = materialproperties[0]
        self.nu = materialproperties[1]

        # Init geodesic variables
        self.startpoints = np.empty((0, 3))
        self.startuv = np.empty((0, 2))
        self.startelements = np.empty(0, dtype=int)
        self.sfiberdirections = np.empty((0, 3))
        self.fiberdirections = np.empty((self.vertexids.shape[0], 3)) * np.nan

        self.surfacenormal = initnormal
        # Init geodesic path record variables
        self.georecord = {}
        self.geoints = []

        # Init uv parameterization parameters
        self.geoparameterization = np.empty((self.vertices.shape[0], 2)) * np.nan
        # Init optimization parameterization
        self.optimizedparameterization = None

        # Calculate compliance tensor
        if isinstance(self.E, list):
            # Orthotropic
            G = materialproperties[2]
            if G is None:
                raise ValueError("G property is not defined.")
            self.compliance_tensor = np.array([[1/self.E[0], -self.nu[0]/self.E[1], 0],
                                               [-self.nu[0]/self.E[0], 1/self.E[1], 0],
                                               [0, 0, 1/(2*G[0])]])
        else:
            # Isotropic
            G = self.E / (2 * (1 + self.nu))
            self.compliance_tensor = np.array([[1/self.E, -self.nu/self.E, 0],
                                               [-self.nu/self.E, 1/self.E, 0],
                                               [0, 0, 1/G]])
        # Calculate stiffness tensor
        self.stiffness_tensor = np.linalg.inv(self.compliance_tensor)
        # Calculate 2D normalized points for each element
        self.normalized_2d = OP.calc2d(self.obj, self.vertices[self.vertexids])

        # Apply angle error to initdirection
        # New initdirection corresponds to a vector rotated around initnormal by angle_error
        initdirection = GEO.rot_vector_angle(initdirection, initnormal, angle_error)

        # Determine if initpoint is on the surface or not
        # If not then we will set initpoint to a projected point on the surface
        t1, t2, projpnt = GEO.find_element_within(initpoint, initdirection, initnormal, self.vertices, self.vertexids, self.facetnormals, self.inplanemat)
        if t1 is None:
            closestvertex_ind = np.where(np.linalg.norm(self.vertices - initpoint, axis=1) == np.min(np.linalg.norm(
                self.vertices - initpoint, axis=1)))
            initpoint = self.vertices[closestvertex_ind[0], :][0]
        elif projpnt is not None:
            initpoint = projpnt
        self.initpoint = initpoint

        # Determine which surface we will be laying the geodesics upon
        self.determine_surface(initpoint, initdirection)

        # Find start points for all geodesics
        # Start from initpoint and go in the positive initdirection and the negative initdirection
        # Starting parameterization point is [0, 0]
        self.find_startpoints(initpoint, initdirection, initnormal, np.array([0.0, 0.0]))

        print("Calculating parameterization...")
        start_time = time.time()

        # Calculate each geodesic defined in self.find_startpoints
        self.calc_geodesics(0)
        # Given the calculated geodesics attempt to create a parameterization
        # We may need to fill in missing elements or elements with few geodesics to get a full parameterization
        self.create_parameterization()

        stop_time = time.time()
        elapsed = stop_time - start_time
        print("\r\nTime to calculate geodesic parameterization: %f seconds" % elapsed)

    def loadobj(self):
        """
        Load a given CAD model using SpatialNDE
        Currently supported models are X3D, STL, and De-La-Mo/DMObjects
        """
        if isinstance(self.cadfile, str):
            if os.path.splitext(self.cadfile)[1] == ".x3d":
                self.objframe = coordframe()
                self.obj = ndepart.fromx3d(self.objframe, None, self.cadfile, tol=1e-6)
            elif os.path.splitext(self.cadfile)[1] in [".stl", ".STL"]:
                self.objframe = coordframe()
                self.obj = ndepart.fromstl(self.objframe, None, self.cadfile, tol=1e-6)
            else:
                raise Exception("Unsupported file type.")
        elif isinstance(self.cadfile, object):
            # print("Loading %s data type." % self.cadfile.__class__.__name__)
            if self.cadfile.__class__.__name__ is "DMObject":
                self.objframe = coordframe()
                self.obj = ndepart.fromDMobject(self.objframe, None, self.cadfile, recalcnormals=False, tol=1e-6)
            else:
                raise Exception("Unsupported object type.")
        else:
            raise Exception("Unsupported data type.")

    def loadvars(self):
        """ Load spatialnde data into the corresponding model variables"""
        self.vertices = self.obj.implpart.surfaces[0].vertices
        self.vertexids_indices = self.obj.implpart.surfaces[0].vertexidx_indices
        self.vertexids = self.obj.implpart.surfaces[0].vertexidx.reshape(self.vertexids_indices.shape[0], 4)[:, 0:3]
        self.facetnormals = self.obj.implpart.surfaces[0].facetnormals
        self.refpoints = self.obj.implpart.surfaces[0].refpoints
        self.inplanemat = self.obj.implpart.surfaces[0].inplanemats
        self.edges = AUV.BuildEdgeDict(self.obj.implpart.surfaces[0])
        self.adjacencyidx = AUV.DetermineAdjacency(self.obj.implpart.surfaces[0], self.edges)
        self.surfaces = AUV.FindTexPatches(self.obj.implpart.surfaces[0], self.adjacencyidx)
        self.obj.implpart.surfaces[0].intrinsicparameterization = None
        self.boxes = self.obj.implpart.surfaces[0].boxes
        self.boxpolys = self.obj.implpart.surfaces[0].boxpolys
        self.boxcoords = self.obj.implpart.surfaces[0].boxcoords

    def find_close_geodesic(self, elements, point):
        """
        Find the closest geodesic in self.georecord and elements

        :param elements: A given set of search elements for a close geodesic
        :param point: The point we want to find a close geodesic to
        :return: tuple(georecord, element closest geodesic is in,
        [vector perpendicular to u direction norm(vector) = v distance, u distance])
        """
        geodets = None
        for element in elements:
            if element in self.georecord.keys():
                # georecord[element] = [pointuv (bary), int_pnt (bary), point (3D), unitfiberdirection (3D), closest_point_idx (idx), uv_start, length, direction]
                geodesics = self.georecord[element][0]
                for g in range(0, len(geodesics)):
                    d2left = GEO.calcdistance(geodesics[g][7]*geodesics[g][3], geodesics[g][2], point)
                    if geodets is None or np.linalg.norm(d2left[0]) < np.linalg.norm(geodets[2][0]):
                        geodets = (geodesics[g], element, d2left)
        if geodets:
            return geodets
        else:
            raise IndexError("Cannot find a close geodesic")

    def determine_surface(self, initpoint, initdirection):
        """
        Determine which surface of a 3D model we should operate on.
        This should single out a solo surface without sharp (90 degree) edges

        :param initpoint: Starting point
        :param initdirection: Starting direction
        :return: Sets self.surface_vertexids
        """
        # Determine if start point is on a vertex or not then perform the necessary calculation to find the next element
        if 0 in np.linalg.norm(self.vertices - initpoint, axis=1):
            # Current point is a vertex:
            element, newvector = GEO.find_element_vertex(initpoint, initdirection, self.surfacenormal,
                                                         self.vertices, self.vertexids, self.facetnormals)
        else:
            # Current point is not on a vertex but within a polygon:
            element, newvector, _ = GEO.find_element_within(initpoint, initdirection, self.surfacenormal,
                                                            self.vertices, self.vertexids, self.facetnormals,
                                                            self.inplanemat)

        surface = [(i, surface.index(element)) for i, surface in enumerate(self.surfaces) if element in surface][0][0]

        surface_polygons = self.surfaces[surface]
        self.surface_vertexids = self.vertexids[surface_polygons]

        # self.vertexids = self.surface_vertexids
        # self.vertices = self.vertices[np.unique(self.surface_vertexids)].reshape(-1, 3)

    def find_startpoints(self, initpoint, initdirection, normal, cfpoint):
        """
        Determines starting location, direction, and element for each geodesic
        Spawns geodesics perpendicular to initdirection
        Geodesic start points are dropped in self.fiberint intervals

        :param initpoint: Starting location (should be close to the center of the model)
        :param initdirection: Starting direction vector
        :param normal: Surface normal vector
        :param cfpoint: initpoint location in parameterization space
        :return: Appends new geodesic start information to the relevant lists
        """
        # We want to travel perpendicular to the fiber direction
        pointdirections = np.array([np.cross(normal, initdirection)])

        # For both positive and negative directions
        directions = [1, -1]
        for i in range(0, pointdirections.shape[0]):
            for direction in directions:
                pointdirection = direction * pointdirections[i]
                point = initpoint

                # Determine if start point is on a vertex or not then perform
                # the necessary calculation to find the next element for the seed geodesic
                if 0 in np.linalg.norm(self.vertices - point, axis=1):
                    # Current point is a vertex:
                    element, newvector = GEO.find_element_vertex(point, pointdirection, normal, self.vertices,
                                                                 self.vertexids, self.facetnormals)
                else:
                    # Current point is not on a vertex but within a polygon:
                    element, newvector, _ = GEO.find_element_within(point, pointdirection, normal, self.vertices,
                                                                    self.vertexids, self.facetnormals, self.inplanemat)

                if 0 in np.linalg.norm(self.vertices - point, axis=1):
                    # Current point is a vertex:
                    selement, _ = GEO.find_element_vertex(point, initdirection, normal, self.vertices, self.vertexids,
                                                          self.facetnormals)
                else:
                    # Current point is not on a vertex but within a polygon:
                    selement, _, _ = GEO.find_element_within(point, initdirection, normal, self.vertices, self.vertexids, self.facetnormals, self.inplanemat)

                if element is None or newvector is None:
                    # print("Point: %s , can't find element or newvector. Element: %s, Newvector: %s" % (point, element, newvector))
                    continue

                if not GEO.check_proj_inplane_pnt(point, self.vertices[self.vertexids[element]]):
                    # print("Point: %s , proj_inplane_pnt failed." % point)
                    continue

                if selement is not None:
                    # Rotate the given fiber vector to be in plane with the start element
                    try:
                        newfibervector = GEO.rot_vector(normal, self.facetnormals[selement], initdirection)
                    except GEO.EdgeError:
                        continue
                    # Set initial point values
                    self.startpoints = np.vstack((self.startpoints, point))
                    self.startuv = np.vstack((self.startuv, cfpoint))
                    self.startelements = np.append(self.startelements, selement)
                    self.sfiberdirections = np.vstack((self.sfiberdirections, newfibervector))
                else:
                    # Can't find a starting element
                    continue

                # Make sure the direction we are going is in plane with the starting element we just found
                pointdirection = GEO.rot_vector(normal, self.facetnormals[element], pointdirection)

                p = 1
                while True:
                    # Set the distance left to travel to the fiber interval
                    dleft = self.fiberint
                    try:
                        while True:
                            # Traverse the current element
                            int_pnt_3d, nextunitvector, nextelement = GEO.traverse_element(self, element, point, pointdirection, None, None, direction=direction, parameterization=False)

                            # Determine how far we are from the calculated intersection point on the other side of the current triangle
                            d2int = np.linalg.norm(int_pnt_3d - point)

                            # Update the distance left to reflect our movement to the edge of the triangle
                            dleft = dleft - d2int

                            # The intersected edge is too far
                            if dleft <= 0:
                                # So our next point is in the current direction but only as far as how much distance we had left to travel
                                point = int_pnt_3d + pointdirection * dleft

                                # Drop down a geodesic here
                                # Rotate the starting direction vector to be in plane with this element
                                selementvec = GEO.rot_vector(self.facetnormals[self.startelements[-1]],
                                                             self.facetnormals[element], self.sfiberdirections[-1],
                                                             force=True)

                                # Append geodesic start point, element, and direction to the corresponding array
                                self.startpoints = np.vstack((self.startpoints, point))
                                self.startelements = np.append(self.startelements, element)
                                self.sfiberdirections = np.vstack((self.sfiberdirections, selementvec))
                                break
                            else:
                                # We still have distance to cover past this intersection point
                                if nextelement is None:
                                    # print("Couldn't find next element.")
                                    raise GEO.EdgeError
                                # Update point, element, and direction to the next triangles details
                                point = int_pnt_3d
                                element = nextelement
                                pointdirection = nextunitvector
                    except GEO.EdgeError:
                        # print("End of geodesic detected. %s" % direction)
                        break
                    # We added a geodesic start point so we must also add its starting u, v point as well
                    self.startuv = np.vstack((self.startuv, np.array([cfpoint[0], direction * self.fiberint * p + cfpoint[1]])))
                    p += 1
        # We couldn't find any start points so we can display the model and useful variables to figure out what we did wrong
        if self.startpoints.shape[0] < 1:
            import matplotlib.pyplot as plt
            from mpl_toolkits.mplot3d import axes3d

            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(self.vertices[:, 0], self.vertices[:, 1], self.vertices[:, 2])
            ax.scatter(initpoint[0], initpoint[1], initpoint[2])
            ax.quiver(initpoint[0], initpoint[1], initpoint[2], initdirection[0], initdirection[1], initdirection[2])
            plt.show()
            raise IndexError("No geodesic start points found.")

    def calc_geodesics(self, startidx):
        """
        Computes the a geodesic path in the positive and negative direction for each point define in start_points beginning
        at the index startidx
        :param: startidx: Which index in self.startpoints to begin calculating geodesics at
        :return: Adds the relevant details for each geodesic to self.georecord
        """
        for i in range(startidx, self.startpoints.shape[0]):
            sys.stdout.write('\r')
            percent_complete = ((i + 1) / float(self.startpoints.shape[0])) * 100.0
            sys.stdout.write("[%-50s] %d%%" % ('=' * int(percent_complete/2), percent_complete))
            sys.stdout.flush()
            self.calc_geodesic(self.startpoints[i], self.startelements[i], self.sfiberdirections[i],
                               self.startuv[i], direction=1, parameterization=True)
            self.calc_geodesic(self.startpoints[i], self.startelements[i], self.sfiberdirections[i],
                               self.startuv[i], direction=-1, parameterization=True)

    def calc_geodesic(self, point, element, unitfiberdirection, uv_start, direction=1, parameterization=False, save_ints=True):
        """
        Calculates the geodesic path from a point in a given direction

        :param point: The geodesic's starting point
        :param element: The first element we will traverse
        :param unitfiberdirection: The desired direction the geodesic will propogate
        :param uv_start: Starting point in uv space
        :param direction: Positive (1) or negative (-1) direction of unitfiberdirection?
        :param parameterization: Do we want to save this geodesic into self.georecord for use in the parameterization?
        :param save_ints: Do we want to save intersection points for plotting purposes?
        :return: Adds the resulting geodesic path to the self.geoints and/or the self.georecord
        """
        int_pnt_3d = point
        unitfiberdirection = direction * unitfiberdirection

        # Create an empty array of intersection points to visualize geodesics
        int_pnts = np.array([point])

        length = 0.0
        p = 0
        while True:
            try:
                int_pnt_3d, nextunitvector, nextelement = GEO.traverse_element(self, element, point,
                                                                                                  unitfiberdirection,
                                                                                                  length, uv_start,
                                                                                                  parameterization=parameterization,
                                                                                                  direction=direction)
            except GEO.EdgeError:
                break

            # Update and store the calculated fiber points and the intersection points
            int_pnts = np.vstack((int_pnts, int_pnt_3d))

            # Calculate the new length of the fiber
            length = direction * np.linalg.norm(int_pnt_3d - point) + length

            if nextelement is None:
                break

            # Update previous iteration values with the next iteration values
            point = int_pnt_3d
            unitfiberdirection = nextunitvector
            element = nextelement
            p += 1

        if save_ints:
            self.geoints.append(int_pnts)
        return length, int_pnt_3d, element

    def check_negative_area(self, record):
        """
        Check to see if we have a triangle with a flipped normal

        :param record: Parameterization we want to check
        :return: True if there is a flipped triangle, False otherwise
        """
        rel_uvw = np.pad(record[self.vertexids], [(0, 0), (0, 0), (0, 1)], "constant", constant_values=1)
        vdir = 0.5 * np.linalg.det(rel_uvw)
        vdir[np.isnan(vdir)] = 0
        if (vdir < 0).any():
            return True
        else:
            return False

    def interpolate_geodesic(self, point, element, minassigned):
        """
        Determine a geodesics starting direction and uv parameterization location that is between two other geodesics

        :param point: A point between two geodesics
        :param element: The element point is within or a vertex of
        :param minassigned: Minimum number of neighbor elements that have geodesics within them for the used starting element
        :return: The direction and uv parameterization location of the geodesic at point
        """
        fiberdirection, cfpoint, shared_cg1, shared_cg2, v = None, None, None, None, None

        # Find neighbors that contain geodesics to the given element
        neighbors = GEO.find_neighbors(element, self.vertexids_indices, self.adjacencyidx)
        neighbors = np.intersect1d(neighbors, self.georecord.keys())

        # If the number of neighbors with geodesics is more than minassigned
        if neighbors.shape[0] > minassigned:
            for neighbor in neighbors:
                # Find a shared vertex between element and neighbor
                sharedvertex = self.vertices[np.intersect1d(self.vertexids[neighbor], self.vertexids[element])][1]
                # Find the closest geodesic in neighbor
                shared_cg, _, distance = self.find_close_geodesic([neighbor], sharedvertex)
                check_dir = np.cross(shared_cg[7] * shared_cg[3], self.facetnormals[element])
                # Spin off two geodesics in opposite directions that are perpendicular to the closest goedeisc in neighbor
                distance1, int_pnt_3d1, element1 = self.calc_geodesic(point, element, check_dir, None, parameterization=False, save_ints=False)
                distance2, int_pnt_3d2, element2 = self.calc_geodesic(point, element, check_dir, None, direction=-1, parameterization=False, save_ints=False)

                ulist = []
                flist = []
                dlist = []
                try:
                    # Determine the u, v distance and direction of the geodesic in the positive direction
                    shared_cg1, _, d1 = self.find_close_geodesic([element1], int_pnt_3d1)
                    u1 = d1[1] + shared_cg1[6] + shared_cg1[5][0]
                    f1 = shared_cg1[7] * shared_cg1[3]
                    ulist.append(u1)
                    flist.append(f1)
                    dlist.append(distance1)
                except IndexError:
                    pass

                try:
                    # Determine the u, v distance and direction of the geodesic in the negative direction
                    shared_cg2, _, d2 = self.find_close_geodesic([element2], int_pnt_3d2)
                    u2 = d2[1] + shared_cg2[6] + shared_cg2[5][0]
                    f2 = shared_cg2[7] * shared_cg2[3]
                    v2 = shared_cg2[5][1]
                    ulist.append(u2)
                    flist.append(f2)
                    dlist.append(distance2)
                except IndexError:
                    v2 = shared_cg1[5][1]
                    pass

                if shared_cg1 is not None and shared_cg2 is not None:
                    num_v = int(np.abs(shared_cg1[5][1] - shared_cg2[5][1]) / self.fiberint)
                    if num_v > 1:
                        v = dlist[-1]
                    else:
                        v = sum(dlist) / len(dlist)
                elif shared_cg1 is not None or shared_cg2 is not None:
                    v = dlist[-1]

                cfpoint = np.array([sum(ulist) / len(ulist), v + v2])

                fiberdirection = sum(flist) / len(flist)
                fiberdirection = GEO.proj_vector(fiberdirection, self.facetnormals[element])
        return fiberdirection, cfpoint

    def interpolate_point(self, vertex):
        """
        Determine the uv coordinates of a point using the same method as interpolate_geodesic
        (Currently not being used as direct assignment is faster and more accurate than spinning off more geodesics)

        :param vertex: Starting vertex
        :return: uv coordinates of vertex
        """
        fiberdirection, cfpoint, shared_cg1, shared_cg2, v = None, None, None, None, None

        vidx = np.where(np.linalg.norm(self.vertices - vertex, axis=1) == np.min(np.linalg.norm(self.vertices - vertex, axis=1)))
        neighbors = np.unique(np.where((self.vertexids == vidx))[0])
        neighbors = np.intersect1d(neighbors, self.georecord.keys())

        shared_cg, element, distance = self.find_close_geodesic(neighbors, vertex)
        check_dir = np.cross(shared_cg[7] * shared_cg[3], self.facetnormals[element])
        distance1, int_pnt_3d1, element1 = self.calc_geodesic(vertex, element, check_dir, None, parameterization=False, save_ints=False)
        distance2, int_pnt_3d2, element2 = self.calc_geodesic(vertex, element, check_dir, None, direction=-1,
                                                              parameterization=False, save_ints=False)

        ulist = []
        dlist = []
        try:
            shared_cg1, _, d1 = self.find_close_geodesic([element1], int_pnt_3d1)
            u1 = d1[1] + shared_cg1[6] + shared_cg1[5][0]
            ulist.append(u1)
            dlist.append(distance1)
        except IndexError:
            pass

        try:
            shared_cg2, _, d2 = self.find_close_geodesic([element2], int_pnt_3d2)
            u2 = d2[1] + shared_cg2[6] + shared_cg2[5][0]
            v2 = shared_cg2[5][1]
            ulist.append(u2)
        except IndexError:
            v2 = shared_cg1[5][1]
            pass

        if shared_cg1 is not None and shared_cg2 is not None:
            num_v = int(np.abs(shared_cg1[5][1] - shared_cg2[5][1]) / self.fiberint)
            if num_v > 1:
                v = dlist[-1]
            else:
                v = sum(dlist) / len(dlist)
        elif shared_cg1 is not None or shared_cg2 is not None:
            v = dlist[-1]

        cfpoint = np.array([sum(ulist) / len(ulist), v + v2])

        return cfpoint

    def fill_missing_geodesics(self, elements, minassigned):
        """
        Spin off more geodesics in elements that contain no geodesics (i.e. fill holes in the initially spawned geodesics)

        :param elements: Elements without an geodesics
        :param minassigned: Minimum number of neighbors that contain geodesics
        :return: Nothing
        """
        # print("Filling missing elements")
        loc = self.startpoints.shape[0]

        for i in elements:
            centroid = self.vertices[self.vertexids[i]].sum(axis=0) / 3
            try:
                num_geos = len(self.georecord[i][0])
            except KeyError:
                # print("No geodesics found in element %s." % i)
                fiberdirection, cfpoint = self.interpolate_geodesic(centroid, i, minassigned)
                if fiberdirection is not None and cfpoint is not None:
                    self.find_startpoints(centroid, fiberdirection, self.facetnormals[i], cfpoint)
                    self.calc_geodesics(loc)
            loc = self.startpoints.shape[0]
        return np.setdiff1d(range(0, self.vertexids.shape[0]), self.georecord.keys()).size

    def fill_low_density_geodesics(self, minassigned):
        """
        Spin off more geodesics in elements that contain a low number of geodesics per area
        (This method is a little unreliable because the geodesics/area threshold isn't well defined for all models)

        :param minassigned: Minimum number of neighbors that contain geodesics
        :return: Nothing
        """
        # print("Filling low density elements")
        loc = self.startpoints.shape[0]

        for i in range(0, self.vertexids.shape[0]):
            centroid = self.vertices[self.vertexids[i]].sum(axis=0) / 3
            try:
                num_geos = len(self.georecord[i][0])
            except KeyError:
                # print("Still can't find a geodesic in element %s." % i)
                continue
            area = 0.5 * np.linalg.det(self.vertices[self.vertexids[i]])
            geos2area = num_geos / area
            if geos2area < 400.0:
                fiberdirection, cfpoint = self.interpolate_geodesic(centroid, i, minassigned)
                if fiberdirection is not None and cfpoint is not None:
                    self.find_startpoints(centroid, fiberdirection, self.facetnormals[i], cfpoint)
                    self.calc_geodesics(loc)
            loc = self.startpoints.shape[0]

    def assign_vertices(self):
        """
        Assign all vertices a uv coordinate based on the closest geodesic to the vertex
        Cleanup methods are employed here such as self.fill_missing_geodesics and self.fill_low_density_geodesics

        :return: Sets the self.geoparameterization
        """
        mask = np.ones((self.geoparameterization.shape[0]), dtype=bool)
        mask[np.unique(self.surface_vertexids)] = False
        for i in range(0, self.vertices.shape[0]):
            vidx = np.where(np.linalg.norm(self.vertices - self.vertices[i], axis=1) == np.min(np.linalg.norm(self.vertices - self.vertices[i], axis=1)))
            neighbors = np.unique(np.where((self.vertexids == vidx))[0])

            for neighbor in neighbors:
                try:
                    closest_geodesic, element, distance = self.find_close_geodesic([neighbor], self.vertices[i])
                except IndexError:
                    continue

                testval = np.dot(calcunitvector(np.cross(closest_geodesic[7]*closest_geodesic[3], distance[0])), self.facetnormals[element])[0]
                fpoint_t = np.array([closest_geodesic[6] + distance[1] + closest_geodesic[5][0], testval*np.linalg.norm(distance[0]) + closest_geodesic[5][1]])
                if np.isnan(fpoint_t[1]):
                    fpoint_t[1] = closest_geodesic[5][1]

                fiberrec = np.copy(self.geoparameterization)
                fiberrec[i] = fpoint_t

                if self.check_negative_area(fiberrec):
                    pass
                else:
                    self.geoparameterization[i] = fpoint_t
        assert not self.check_negative_area(self.geoparameterization)
        # assert np.where((np.isnan(self.geoparameterization).all(axis=1) & np.array(~mask)))[0].size == 0

    def create_parameterization(self):
        """
        Create the uv parameterization based on the computed geodesic paths
        If vertices can't be assigned or are missed cleanup methods are employed to attempt to solve coordinates for all
        vertices. This is difficult to make robust as geometric complexities can vary quite largely.
        """
        mask = np.ones((self.geoparameterization.shape[0]), dtype=bool)
        mask[np.unique(self.surface_vertexids)] = False

        missed_elements = np.setdiff1d(range(0, self.vertexids.shape[0]), self.georecord.keys())
        if missed_elements.size > 0:
            for minassigned in [1, 0]:
                # Fill in any elements that had been missed by the initial spawing of geodesics
                missing = self.fill_missing_geodesics(missed_elements, minassigned)
                if missing == 0:
                    break

        # Attempt to assign uv coordinates to all vertices
        self.assign_vertices()

        if np.where((np.isnan(self.geoparameterization).all(axis=1) & np.array(~mask)))[0].size > 0:
            # Couldn't assign all points so lets add some more geodesics in elements where the geodesic density is low
            self.fill_low_density_geodesics(0)
            self.assign_vertices()

        # Basic interpolation of uv coordinates based on nearby element fiber directions
        self.interpolate(np.where((np.isnan(self.geoparameterization).all(axis=1) & np.array(~mask)))[0], mask)

        # Last ditch attempt of merely averaging nearby uv coordinates
        self.average_fpoint(np.where((np.isnan(self.geoparameterization).all(axis=1) & np.array(~mask)))[0], mask)

        # Guarantee that we don't have any flipped triangles or that we missed any vertices
        # If vertices are missed consider trying other angle_error values
        assert not self.check_negative_area(self.geoparameterization)
        assert np.where((np.isnan(self.geoparameterization).all(axis=1) & np.array(~mask)))[0].size == 0

    def interpolate(self, leftover_idxs, mask):
        """
        Basic interpolation of uv coordinates based on nearby element fiber directions leftover from geodesic paths

        :param leftover_idxs: Indices of missed vertices during assignment
        :param mask: Surface vertex index mask
        :return: Any vertex indices that were missed by this cleanup method
        """
        timeout = 0
        while leftover_idxs.shape[0] > 0 and timeout < 50:
            for i in range(0, leftover_idxs.shape[0]):
                unassigned_facets = np.unique(np.where((self.vertexids == leftover_idxs[i]))[0])
                done = False
                for j in unassigned_facets:
                    elementvertices = self.vertices[self.vertexids[j]]
                    assigned_vertsids = np.where((~np.isnan(self.geoparameterization[self.vertexids[j]]).all(axis=1)))[0]
                    if assigned_vertsids.shape[0] > 0:
                        if not np.isnan(self.fiberdirections[j]).any():
                            fiberdirection = self.fiberdirections[j]
                            for k in assigned_vertsids:
                                assigned_fpoint = self.geoparameterization[self.vertexids[j]][k]
                                fdistance, closest_point = GEO.calcclosestpoint(fiberdirection, elementvertices[k],
                                                                                np.array([self.vertices[leftover_idxs[i]]]),
                                                                                self.facetnormals[j])
                                unassigned_fpoint = assigned_fpoint + fdistance

                                fiberrec = np.copy(self.geoparameterization)
                                fiberrec[leftover_idxs[i]] = unassigned_fpoint

                                if self.check_negative_area(fiberrec):
                                    pass
                                else:
                                    self.geoparameterization[leftover_idxs[i]] = unassigned_fpoint
                                    done = True
                                    break
                    if done:
                        break
            leftover_idxs = np.where((np.isnan(self.geoparameterization).all(axis=1) & np.array(~mask)))[0]
            timeout += 1
        return leftover_idxs

    def average_fpoint(self, leftover_idxs, mask):
        """
        Last ditch effort to assign the missed vertices. Simply average any nearby uv coordinates and make sure the new
        coordinate won't flip the triangle.

        :param leftover_idxs: Indices of missed vertices during assignment
        :param mask: Surface vertex index mask
        :return: Any vertex indices that were missed by this cleanup method
        """
        for i in range(0, leftover_idxs.shape[0]):
            neighbors = np.unique(np.where((self.vertexids == leftover_idxs[i]))[0])
            neighbor_fpoint = self.geoparameterization[self.vertexids[neighbors]]
            neighbor_fpoint = neighbor_fpoint.reshape(neighbor_fpoint.shape[0]*neighbor_fpoint.shape[1], 2)

            count = ~np.isnan(neighbor_fpoint).all(axis=1)

            average_fpoint = np.sum(neighbor_fpoint[count], axis=0)/neighbor_fpoint[count].shape[0]

            fiberrec = np.copy(self.geoparameterization)
            fiberrec[leftover_idxs[i]] = average_fpoint

            if self.check_negative_area(fiberrec):
                pass
            else:
                self.geoparameterization[leftover_idxs[i]] = average_fpoint

            del fiberrec
        leftover_idxs = np.where((np.isnan(self.geoparameterization).all(axis=1) & np.array(~mask)))[0]
        return leftover_idxs

    def layup(self, angle, orientation_locations=None, precision=1e-4, maxsteps=10000, lr=1e-3, decay=0.7, eps=1e-8, mu=0.8, plotting=False, save=False):
        """
        Once the parameterization has been computed we can calculate any fiber orientation without needing to compute a new geodesic
        mapping. Simply rotate the geoparameterization by angle and then minimize the strain energy.

        :param angle: Desired fiber orientation
        :param orientation_locations: Optional locations to calculate the fiber orientations at. Default is element centroids.
        :param precision: Termination threshold for the strain energy optimization
        :param maxsteps: Maximum number of optimizaiton iterations
        :param lr: Optimization rate (similar to learning rate in machine learning optimizers)
        :param decay: Decay rate
        :param eps: Error precision value
        :param mu: Momentum value
        :param plotting: Do we want to plot the results using matplotlib?
        :param save: Do we want to save the fiber orintations at orientation_locations to a .npy file
        :return: texcoords2inplane - Transformation matrix between 3D space and 2D space, once created any orientation at any
        point on the surface can be evaluated.
        """
        orientations = None
        # Remove the angle_error used to create the geodesic parameterization from these calculations
        angle -= self.error
        # Simple rotation matrix by angle
        rmatrix = np.array([[np.cos(np.deg2rad(angle)), -np.sin(np.deg2rad(angle))],
                            [np.sin(np.deg2rad(angle)), np.cos(np.deg2rad(angle))]])
        # Rotate the geodesic parameterization by angle
        parameterization = np.matmul(self.geoparameterization, rmatrix)

        # Optimize the rotate parameterization using an RMSprop optimizer. This will minimize the strain energy between the
        # 3D mesh surface and the uv mesh.
        optimizedparamterization, loss = self.fiberoptimize(parameterization, precision=precision, maxsteps=maxsteps, lr=lr, decay=decay, eps=eps, mu=mu)
        # Compute the transformation between the 3D model and the optimized uv mesh
        texcoords2inplane = self.calctransform(optimizedparamterization)

        if orientation_locations is not None:
            orientations = self.calcorientations_abaqus(orientation_locations, self.vertices, self.vertexids, self.inplanemat,
                                                        texcoords2inplane, self.obj.implpart.surfaces[0].boxes,
                                                        self.obj.implpart.surfaces[0].boxpolys,
                                                        self.obj.implpart.surfaces[0].boxcoords)

        if plotting:
            if orientations is None:
                # Calculate orientations at the centroid of each element
                orientation_locations = self.vertices[self.vertexids].sum(axis=1) / 3
                orientations = self.calcorientations_abaqus(orientation_locations, self.vertices, self.vertexids, self.inplanemat,
                                                            texcoords2inplane, self.obj.implpart.surfaces[0].boxes,
                                                            self.obj.implpart.surfaces[0].boxpolys,
                                                            self.obj.implpart.surfaces[0].boxcoords)

            if save:
                np.save("orientation_%s.npy" % angle, orientations)

            import matplotlib.pyplot as plt
            from mpl_toolkits.mplot3d import axes3d

            self.plot_geodesics()

            fig = plt.figure()
            plt.plot(range(len(loss)), loss)
            plt.title("Loss Function")
            plt.xlabel("Iteration")
            plt.ylabel("Energy (J/length)")

            fig = plt.figure()
            plt.scatter(parameterization[:, 0], parameterization[:, 1], c="b")
            plt.scatter(optimizedparamterization[:, 0], optimizedparamterization[:, 1], c="orange")
            plt.title("Parameterizations")
            plt.ylabel("V")
            plt.xlabel("U")
            plt.legend(["Original", "Optimized"])

            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(self.vertices[:, 0], self.vertices[:, 1], self.vertices[:, 2])
            ax.quiver(orientation_locations[:, 0], orientation_locations[:, 1], orientation_locations[:, 2], orientations[:, 0], orientations[:, 1], orientations[:, 2], arrow_length_ratio=0, length=1.0)
            plt.show()

        return texcoords2inplane

    def fiberoptimize(self, seed, precision=None, maxsteps=None, lr=None, decay=None, eps=None, mu=None):
        """
        Minimum strain energy optimization of a seed parameterization. Uses an RMSprop optimization algorithm.

        :param seed: Initial uv parameterization
        :param precision: Termination threshold for the strain energy optimization
        :param maxsteps: Maximum number of optimizaiton iterations
        :param lr: Optimization rate (similar to learning rate in machine learning optimizers)
        :param decay: Decay rate
        :param eps: Error precision value
        :param mu: Momentum value
        :return:
        """
        # Define the strain energy function
        def f(x, *args):
            return OP.computeglobalstrain(self.normalized_2d, x, self.vertexids, self.stiffness_tensor)

        # Define the strain energy gradient function
        def gradf(x, *args):
            oc = np.argmin(np.linalg.norm(self.geoparameterization, axis=1))
            return OP.computeglobalstrain_grad(self.normalized_2d, x, self.vertexids, self.stiffness_tensor, oc)

        print("Optimizing...")
        initenergy = OP.computeglobalstrain(self.normalized_2d, seed.flatten(), self.vertexids, self.stiffness_tensor)
        print("Initial strain energy: %s J/(model length)" % initenergy)

        start_time = time.time()

        # print("Optimizing with rmsprop...")
        optimizedparameterization, loss = OP.rmsprop_momentum(f, gradf, seed, precision=precision, maxsteps=maxsteps, lr=lr, decay=decay, eps=eps, mu=mu)

        stop_time = time.time()
        elapsed = stop_time - start_time
        print("Time to optimize: %f seconds" % elapsed)

        return optimizedparameterization, loss

    def calctransform(self, parameterization):
        """
        Calculate the transformation matrix between the 3D model and the uv parameterization using spatialNDE

        :param parameterization: uv parameterization
        :return: Transformation matrix
        """
        self.obj.implpart.surfaces[0].intrinsicparameterization = polygonalsurface_texcoordparameterization.new(self.obj.implpart.surfaces[0], parameterization, self.obj.implpart.surfaces[0].vertexidx, None)
        self.obj.implpart.surfaces[0].intrinsicparameterization.buildprojinfo(self.obj.implpart.surfaces[0])

        texcoords2inplane = self.obj.implpart.surfaces[0].intrinsicparameterization.texcoords2inplane

        return texcoords2inplane

    def plot_geodesics(self):
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import axes3d

        mask = np.ones((self.geoparameterization.shape[0]), dtype=bool)
        mask[np.unique(self.surface_vertexids)] = False
        leftover_idxs = np.where((np.isnan(self.geoparameterization).all(axis=1) & np.array(~mask)))[0]

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(self.initpoint[0], self.initpoint[1], self.initpoint[2])
        ax.scatter(self.vertices[:, 0], self.vertices[:, 1], self.vertices[:, 2], alpha=0.1)
        ax.scatter(self.startpoints[:, 0], self.startpoints[:, 1], self.startpoints[:, 2])
        # ax.scatter(self.vertices[leftover_idxs][:, 0], self.vertices[leftover_idxs][:, 1], self.vertices[leftover_idxs][:, 2], c="r")
        # ax.quiver(self.startpoints[:, 0], self.startpoints[:, 1], self.startpoints[:, 2], self.sfiberdirections[:, 0], self.sfiberdirections[:, 1], self.sfiberdirections[:, 2])
        for i in self.geoints:
            ax.plot(i[:, 0], i[:, 1], i[:, 2])

        fig = plt.figure()
        for i in range(0, self.vertexids.shape[0]):
            plt.plot(self.geoparameterization[self.vertexids[i]][:, 0], self.geoparameterization[self.vertexids[i]][:, 1])
        plt.scatter(self.startuv[:, 0], self.startuv[:, 1])

    def calcorientations_abaqus(self, modellocs, vertices, vertexids, inplanemat, texcoords2inplane, boxes, boxpolys, boxcoords):
        """
        Function optimized for Abaqus to determine the orientations of given locations.

        :param modellocs: Points of interest on the surface or close to the surface
        :param vertices: Vertices of the model
        :param vertexids: Indices of the vertices of the model for each element
        :param inplanemat: Orthogonal matrix define by spatialnde
        :param texcoords2inplane: Transformation matrix between 3D space and uv space
        :param boxes: spatialnde polygon box definition
        :param boxpolys: spatialnde polygon definition relative to boxes
        :param boxcoords: spatialnde box vertex coordinates
        :return: Orientations at modellocs
        """
        orientations = np.empty((modellocs.shape[0], 3))
        orientations[:] = np.nan

        rangex = np.array([boxcoords[:, 0], boxcoords[:, 3]]).T
        rangey = np.array([boxcoords[:, 1], boxcoords[:, 4]]).T
        rangez = np.array([boxcoords[:, 2], boxcoords[:, 5]]).T

        def in_box(p):
            containers = np.where(np.logical_and(np.logical_and(
                np.logical_and(p[0] > rangex[:, 0], p[0] < rangex[:, 1]),
                np.logical_and(p[1] > rangey[:, 0], p[1] < rangey[:, 1])),
                np.logical_and(p[2] > rangez[:, 0], p[2] < rangez[:, 1])))[0]
            bindx = boxes[containers][:, -1][boxes[containers][:, -1] != -1][0]

            polyslist = np.array([boxpolys[bindx]])
            count = 1
            while True:
                if boxpolys[bindx+count] == -1:
                    break
                polyslist = np.append(polyslist, boxpolys[bindx+count])
                count += 1

            return bindx, polyslist

        for i in range(0, modellocs.shape[0]):
            idx, element = None, None

            point = modellocs[i]

            try:
                box_idx, polys = in_box(point)
            except IndexError:
                idx = np.where(np.linalg.norm(vertices - point, axis=1) == np.min(np.linalg.norm(vertices - point, axis=1)))
                vert = vertices[idx][0]
                box_idx, polys = in_box(vert)

            for j in polys:
                check, projpnt = GEO.check_proj_inplane_pnt(point, vertices[vertexids[j]])
                if self.point_in_polygon_3d(vertices[vertexids][j], projpnt, inplanemat[j]):
                    element = j
                    break

            if element is None:
                if idx is None:
                    idx = np.where(np.linalg.norm(vertices - point, axis=1) == np.min(np.linalg.norm(vertices - point, axis=1)))
                vertneighbors = np.unique(np.where((vertexids == idx))[0])
                element = vertneighbors[0]

            if element is not None:
                red_texcoords2inplane = texcoords2inplane[element][:2, :2]
                texutexbasis = np.array([1.0, 0.0])
                texu2dbasis = np.dot(red_texcoords2inplane, texutexbasis)
                u3D = np.dot(inplanemat[element].T, texu2dbasis)
                orientations[i] = calcunitvector(u3D)
            else:
                print("Failed to find point on surface: %s" % point)
        return orientations

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
