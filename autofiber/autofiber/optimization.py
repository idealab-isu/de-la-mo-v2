# - * - coding: utf - 8 - * -
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


def calcunitvector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)


def calc2d(obj, points):
    """
    Calculate a 2D representation of a 3D model

    :param obj: spatialnde object
    :param points: 3D model points to be converted to 2D
    :return: points in 2D space
    """
    coord_sys = obj.implpart.surfaces[0].inplanemats
    coord_sys = np.transpose(coord_sys, axes=(0, 2, 1))
    points_2d = np.matmul(points, coord_sys)
    return points_2d


def minor(arr, i, j):
    """
    https://stackoverflow.com/questions/3858213/numpy-routine-for-computing-matrix-minors
    Calculate the minor of a matrix with ith row, jth column removed

    :param arr: Matrix of interest
    :param i: row to remove
    :param j: column to remove
    :return: minor of arr with ith row removed and jth column removed
    """
    # ith row, jth column removed
    return ((-1) ** (i + j)) * arr[:, np.array(list(range(i))+list(range(i+1, arr.shape[1])))[:, np.newaxis],
                               np.array(list(range(j))+list(range(j+1, arr.shape[2])))]


def build_checkerboard(w, h):
    """
    https://stackoverflow.com/questions/2169478/how-to-make-a-checkerboard-in-numpy
    Build a checkerboard array

    :param w: width of checkerboard
    :param h: height of checkerboard
    :return: checkerboard array of width w and height h
    """
    re = np.r_[w * [1, -1]]  # even-numbered rows
    ro = np.r_[w * [-1, 1]]  # odd-numbered rows
    return np.row_stack(h * (re, ro))[:w, :h]


def computeglobalstrain(normalized_2d, fiberpoints, vertexids, stiffness_tensor):
    r"""
    Compute the strain energy between a 2d representation of a surface and a uv (geodesic) parameterization

    :math:`p_{3D}` - 2D mapping of 3D model surface based on element normal and inplanemat

    :math:`p_{UV}` - Parameterization of 3D surface based on geodesic lines

    Compute the areas of each triangular element in UV space:

    :math:`A_{uv} = \frac{1}{2}det(p_{uv})`

    Compute the deformation gradient between each element in UV space and the corresponding element in 2D space:

    :math:`F = p_{3D} * p_{UV}^{-1}`

    Utilizing the Lagrangian finite strain tensor:

    :math:`\epsilon = \frac{1}{2}(C - I)`

    where :math:`C = F^{t}F` the right Cauchy–Green deformation tensor and :math:`I` is the identity matrix.

    Since :math:`\overrightarrow{\epsilon} = [\epsilon_{11}, \epsilon_{22}, 2\gamma_{12}]` we have to:

    :math:`\overrightarrow{\epsilon} = \frac{\epsilon}{[1.0, 1.0, 0.5]}^T`

    Therefore the total strain energy density of the surface with units (J/length):

    :math:`E_{total} = \sum_{k=0}^N\frac{1}{2}\overrightarrow{\epsilon_k}^2SA_{uv_k}`

    where :math:`S` is the stiffness tensor of the defined material and N is the total number of mesh elements.

    :param normalized_2d: 2D representation of a 3D model
    :param fiberpoints: uv parameterization
    :param vertexids: Vertex indices of each element in the 3D model
    :param stiffness_tensor: Stiffness tensor of the given material
    :return: The computed total strain energy between normalized_2d and fiberpoints
    """
    element_vertices_uv = fiberpoints.reshape(int(fiberpoints.shape[0]/2), 2)[vertexids]

    centroid_2d = np.sum(normalized_2d, axis=1) / 3
    centroid_uv = np.sum(element_vertices_uv, axis=1) / 3

    rel_uv = np.subtract(element_vertices_uv, centroid_uv[:, np.newaxis])
    rel_2d = np.subtract(normalized_2d, centroid_2d[:, np.newaxis])

    rel_uvw = np.pad(rel_uv, [(0, 0), (0, 0), (0, 1)], "constant", constant_values=1).transpose(0, 2, 1)
    rel_3d = np.pad(rel_2d, [(0, 0), (0, 0), (0, 1)], "constant", constant_values=1).transpose(0, 2, 1)

    # Compute the area of each element
    areas = 0.5 * np.linalg.det(rel_uvw)

    # Compute the deformation matrix for each element
    F = np.matmul(rel_3d, np.linalg.inv(rel_uvw))[:, :2, :2]

    # We can exclude the rotation of F by multiplying by it's transpose and gaining the strain
    # https://en.wikipedia.org/wiki/Finite_strain_theory
    strain = 0.5 * (np.matmul(F.transpose(0, 2, 1), F) - np.identity(F.shape[1]))

    m = np.array([1.0, 1.0, 0.5])[np.newaxis].T
    strain_vector = np.divide(np.array([[strain[:, 0, 0]], [strain[:, 1, 1]], [strain[:, 0, 1]]]).transpose((2, 0, 1)), m).squeeze()

    # http://homepages.engineering.auckland.ac.nz/~pkel015/SolidMechanicsBooks/Part_I/BookSM_Part_I/08_Energy/08_Energy_02_Elastic_Strain_Energy.pdf
    # J == Pa * m^3 -> J/m = Pa * m^2           ex. GPa * mm^3 = 1 J so to get J/mm -> GPa * mm^2
    strain_energy_density = 0.5*np.multiply(np.einsum("ei,ei->e", strain_vector, np.matmul(strain_vector, stiffness_tensor)), areas)

    total_strain_energy = np.sum(strain_energy_density)
    return total_strain_energy


# Create a helper matrix now to use later
duvw_duij_t = np.zeros((6, 3, 3))
for j in range(0, 3):
    for i in range(0, 2):
        duvw_duij_t[j*2+i, i, j] = 1


def computeglobalstrain_grad(normalized_2d, fiberpoints, vertexids, stiffness_tensor, oc):
    r"""
    Compute the gradient of the strain energy function defined above with respect to the movement of each point in the
    uv parameterization.

    :math:`p_{3D}` - 2D mapping of 3D model surface based on element normal and inplanemat

    :math:`p_{UV}` - Parameterization of 3D surface based on geodesic lines

    Compute the areas of each triangular element in UV space:

    :math:`A_{uv} = \frac{1}{2}det(p_{uv})`

    In order to calculate the derivative of the area of each element with respect to each nodal displacement we will
    calculate the (i,j)-minor of :math:`A_{uv}` by the determinate of the matrix created by removing the ith row and jth
    column in :math:`p_{uv}` for each element in :math:`p_{uv}`:

    :math:`M_{ij} = det(minor(p_{uv_{ij}}))`

    Then using :math:`M_{ij}` we can compute the cofactor matrix:

    :math:`Cof = ((-1)^{i+j}M_{ij})_{1 \leq i,j \leq n^*}`

    Therefore, the adjugate matrix of :math:`A_{uv}` is:

    :math:`adj(A_{uv}) = Cof^T`

    The derivative of the area of each element with respect to each nodal displacement can be calculated using
    Jacobi's formula as follows:

    :math:`\frac{dA_{uv}}{dp_{uv}} = -\frac{1}{2}trace(adj(A_{uv})*\frac{dp_{uv}}{dp_{uv_{ij}}})`

    where :math:`\frac{dp_{uv}}{dp_{uv_{ij}}}` is the derivative of each nodal displacement with respect to moving all
    the other nodes for each mesh element.

    Compute the deformation gradient between each element in UV space and the corresponding element in 2D space:

    :math:`F = p_{3D} * p_{UV}^{-1}`

    Utilizing the Lagrangian finite strain tensor:

    :math:`\epsilon = \frac{1}{2}(C - I)`

    where :math:`C = F^{t}F` the right Cauchy–Green deformation tensor and :math:`I` is the identity matrix.

    Since :math:`\overrightarrow{\epsilon} = [\epsilon_{11}, \epsilon_{22}, \gamma_{12}/2]` we have to:

    :math:`\overrightarrow{\epsilon} = \frac{\epsilon}{[1.0, 1.0, 0.5]}^T`

    The derivative of the deformation tensor with respect to each nodal displacement is as follows:

    :math:`\frac{dF}{dp_{uv}} = p_{uv}^{-1}\frac{dp_{uv}}{dp_{uv_{ij}}}p_{uv}^{-1}p_{3D}`

    The derivative of strain with respect to each nodal displacement:

    :math:`\frac{d\epsilon}{dp_{uv}} = \frac{1}{2}(\frac{dF}{dp_{uv}}^TF + F^T\frac{dF}{dp_{uv}})`

    Then to account for :math:`2\gamma_{12}`:

    :math:`\frac{d\overrightarrow{\epsilon}}{dp_{uv}} = \frac{\frac{d\epsilon}{dp_{uv}}}{[1.0, 1.0, 0.5]}`

    Finally, we can compute the derivative of the strain energy density with respect to each nodal displacement with the
    following application of the chain rule:

    :math:`\frac{dE}{dp_{uv}} = \overrightarrow{\epsilon}S\frac{d\overrightarrow{\epsilon}}{dp_{uv}}A_{uv} + \frac{1}{2}\overrightarrow{\epsilon}S\overrightarrow{\epsilon}\frac{dA_{uv}}{dp_{uv}}`

    :param normalized_2d: 2D representation of a 3D model
    :param fiberpoints: uv parameterization
    :param vertexids: Vertex indices of each element in the 3D model
    :param stiffness_tensor: Stiffness tensor of the given material
    :param oc: A vertex index which we want to constrain by fixing it's location
    :return: The gradient of strain energy with respect to movement of each point in the uv parameterization
    """
    element_vertices_uv = fiberpoints.reshape(int(fiberpoints.shape[0]/2), 2)[vertexids]

    centroid_2d = np.sum(normalized_2d, axis=1) / 3
    centroid_uv = np.sum(element_vertices_uv, axis=1) / 3

    rel_uv = np.subtract(element_vertices_uv, centroid_uv[:, np.newaxis])
    rel_2d = np.subtract(normalized_2d, centroid_2d[:, np.newaxis])

    rel_uvw = np.pad(rel_uv, [(0, 0), (0, 0), (0, 1)], "constant", constant_values=1).transpose(0, 2, 1)
    rel_3d = np.pad(rel_2d, [(0, 0), (0, 0), (0, 1)], "constant", constant_values=1).transpose(0, 2, 1)

    # Compute the area of each element
    areas = 0.5 * np.linalg.det(rel_uvw)

    # Compute the minor of each UV element by calculating the determinate of the minor resulting from deleting the ith
    # row and jth column
    minor_mat = np.zeros(rel_uvw.shape)
    for i in range(0, 3):
        for j in range(0, 3):
            minor_mat[:, i, j] = np.linalg.det(minor(rel_uvw, i, j))

    # Compute the adjugate matrix
    # https://en.wikipedia.org/wiki/Adjugate_matrix
    adj_mat = np.multiply(minor_mat, build_checkerboard(minor_mat.shape[1], minor_mat.shape[2])).transpose(0, 2, 1)

    # Compute the derivative of the area with respect to moving each uv point using Jacobi's formula
    # https://en.wikipedia.org/wiki/Jacobi%27s_formula
    # area = 0.5*det(UV) -> dareas_duv = trace(adj(UV)*dUV_{dUV_ij})
    dareas_duv = -0.5*np.trace(np.matmul(adj_mat[:, np.newaxis, :, :], duvw_duij_t), axis1=2, axis2=3)

    # Compute the deformation matrix
    F = np.matmul(rel_3d, np.linalg.inv(rel_uvw))[:, :2, :2]

    # We can exclude the rotation of F by multiplying by it's transpose and gaining the strain
    # https://en.wikipedia.org/wiki/Finite_strain_theory
    strain = 0.5 * (np.matmul(F.transpose(0, 2, 1), F) - np.identity(F.shape[1]))

    m = np.array([1.0, 1.0, 0.5])[np.newaxis].T
    strain_vector = np.divide(np.array([[strain[:, 0, 0]], [strain[:, 1, 1]], [strain[:, 0, 1]]]).transpose((2, 0, 1)), m).squeeze()

    # compute the derivative of the deformation matrix with respect to each uv coordinate
    # 3D*F = UV -> UV^-1*dUV_{dUV_ij}*UV^-1*3D = dF_duv
    dF_duv = np.matmul(rel_3d[:, np.newaxis, :, :], np.matmul(np.matmul(np.linalg.inv(rel_uvw)[:, np.newaxis, :, :], duvw_duij_t), np.linalg.inv(rel_uvw)[:, np.newaxis, :, :]))[:, :, :2, :2]

    # Compute the derivative of the strain with respect to each uv coordinate
    dstrainvector_duv = np.zeros((strain_vector.shape[0], strain_vector.shape[1], 6))
    for i in range(0, 6):
        dstrain_du = 0.5 * (np.matmul(dF_duv[:, i, :, :].transpose(0, 2, 1), F) + np.matmul(F.transpose(0, 2, 1), dF_duv[:, i, :, :]))
        dstrainvector_duv[:, :, i] = np.divide(np.array([[dstrain_du[:, 0, 0]], [dstrain_du[:, 1, 1]], [dstrain_du[:, 0, 1]]]).transpose((2, 0, 1)), m).squeeze()

    # Using the chain rule we can compute the derivative of the strain energy with respect to each uv coordinate
    # dE_du = strain*C*dstrain_duv*areas + 0.5*strain*C*strain*dareas_duv
    dE_du = (np.einsum("ei,e->ei", np.einsum("ei,eij->ej", np.matmul(strain_vector, stiffness_tensor), dstrainvector_duv), areas) + 0.5*np.einsum("e,ej->ej", np.einsum("ei,ei->e", np.matmul(strain_vector, stiffness_tensor), strain_vector), dareas_duv)).reshape(dstrainvector_duv.shape[0], 3, 2)

    point_strain_grad = np.zeros((int(fiberpoints.shape[0]/2), 2))
    for i in range(0, vertexids.shape[0]):
        ele_vertices = vertexids[i]
        ele_strain_grad = dE_du[i]

        point_strain_grad[ele_vertices] = point_strain_grad[ele_vertices] + ele_strain_grad

    point_strain_grad[oc][0] = 0.0
    point_strain_grad[oc][1] = 0.0

    return -1*point_strain_grad.flatten()


# https://medium.com/100-days-of-algorithms/day-69-rmsprop-7a88d475003b
def rmsprop_momentum(F, dF, x_0, precision=None, maxsteps=None, lr=None, decay=None, eps=None, mu=None):
    x = x_0.flatten()
    loss = []
    dx_mean_sqr = np.zeros(x.shape, dtype=float)
    momentum = np.zeros(x.shape, dtype=float)

    if F(x) < precision:
        return x.reshape(-1, 2), loss

    for _ in range(maxsteps):
        b0 = F(x)
        dx = dF(x)
        dx_mean_sqr = decay * dx_mean_sqr + (1 - decay) * dx ** 2
        momentum = mu * momentum + lr * dx / (np.sqrt(dx_mean_sqr) + eps)
        x -= momentum

        loss.append(F(x))

        sys.stdout.write("Residual: %s      \r" % abs(F(x) - b0))
        sys.stdout.flush()

        if abs(F(x) - b0) < precision:
            break
        if F(x) < 0:
            raise ValueError("Negative strain energy detected. Bad parameterization?")
    print("")
    return x.reshape(-1, 2), loss
