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


import pickle
import numpy as np
from autofiber.generator import AutoFiber as AF

# Define the material properties of the composite fiber material. Since all of the models used here are using mm as
# as distance units then our material properties must be in GPa in order to get strain energy in J/mm. See
# optimization.py line 99.
materialproperties = (
    # Young's Modulus [E1, E2, E3]
    [1.415e11/1.e9, 8.5e9/1.e9, 8.5e9/1.e9],
    # Poisson's ratio [nu12, nu13, nu23]
    [0.33, 0.33, 0.33],
    # Shear Modulus [G12, G13, G23]
    [5.02e9/1.e9, 5.02e9/1.e9, 2.35e9/1.e9])

# angles = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0]
angles = [0.0]

for angle in angles:
    print("Angle: %s" % angle)
    # Flat plate - X3D - Anisotropic material properties
    # Path to X3D CAD model
    # test = AF('demos/FlatPlate.x3d',
    #           # Point close to the center of the model
    #           np.array([1.0, -1.0, 0.0]),
    #           # Vector pointing in fiber direction
    #           np.array([-np.cos(np.deg2rad(angle)), np.sin(np.deg2rad(angle)), 0]),
    #           # Vector pointing in the primary normal to the surface
    #           np.array([0, 0, 1]),
    #           # Material properties of the composite fiber material
    #           materialproperties=materialproperties,
    #           # Geodesic spacing interval
    #           fiberint=0.01)
    # Compute an orientation given a fiber angle
    # texcoords2inplane = test.layup(0.0, plotting=True)

    # Cylinder - X3D - Anisotropic material properties
    # test = AF('demos/Cylinder.x3d',
    #           np.array([0.0, 0.0, 1.0]),
    #           np.array([np.sin(-1 * np.deg2rad(angle)), np.cos(np.deg2rad(angle)), 0.0]),
    #           np.array([0, 0, 1.0]),
    #           materialproperties=materialproperties,
    #           fiberint=0.05)
    # texcoords2inplane = test.layup(0.0, plotting=True)

    # Saddle model - X3D - Anisotropic material properties
    test = AF('demos/SmallSaddle32.x3d',
              np.array([0.0, 1.0, 0.0]),
              np.array([-np.cos(np.deg2rad(angle)), 0, -np.sin(np.deg2rad(angle))]),
              np.array([0, 1, 0]),
              materialproperties=materialproperties,
              fiberint=0.01,
              # Further angle to be applied to the given fiber direction
              #  This is typically used as a tolerance in the fiber angle but
              #  can be used to rotate the base parameterization such as in the case
              #  of the saddle part we can actually define a parameterization easily at 45.0 degrees.
              #  Recommended to be any number except exactly zero
              angle_error=45.0)
    texcoords2inplane = test.layup(0.0, plotting=True)

    # Curved Composite Mold / De-La-Mo - DMObject - Anisotropic material properties
    # test = AF(pickle.load(open('demos/DMObject.pkl', 'rb')),     # mm
    #           np.array([0.0, 0.0, 3.0]),
    #           np.array([np.cos(np.deg2rad(angle)), np.sin(np.deg2rad(angle)), 0.0]),
    #           np.array([0, 0, 1]),
    #           materialproperties=materialproperties,
    #           fiberint=0.25)
    # Can define mesh coordinates at which the orientation can be computed
    # meshcoords = np.load("demos/curved_abaqus_mesh_coords.npy")
    # texcoords2inplane = test.layup(0.0, orientation_locations=meshcoords, plotting=True)
