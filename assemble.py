from numpy import eye, zeros, ones, multiply, s_, linspace, ndindex, cos, pi
from numpy.linalg import solve as np_solve
from numpy.linalg import cond
from numpy import sum as np_sum
from numpy import atleast_2d as atleast2d
from numpy import vectorize
from utils.shape_functions import shape_functions_dict
from utils.element_def import element
from utils.constants import quad_pts, POLYNOMIAL, penalty_factor
from utils.boundary_conditions import has_essential_boundary_condition
from utils.boundary_conditions import has_neumann_boundary_condition
from utils.boundary_conditions import set_essential_boundary_conditions
from quadrature import coord_transform
import matplotlib.pyplot as plt


def construct_stiffness(
    elements,
    polynomial_order,
    alpha,
    beta,
    gamma,
    k,
    b,
    c
):
    #
    #
    quad_data = quad_pts[3]
    #
    # Specify the basis functions
    shape_functions = shape_functions_dict[polynomial_order]

    #
    # First construct the naive stiffness matrix
    # which does not account for boundary conditions
    # or reduction necessary to prevent singularity
    num_elements = len(elements)
    num_nodes = num_elements * polynomial_order + 1
    K_naive = zeros((num_nodes, num_nodes))
    #
    # Identify the locations of the 
    # submatrices in the master matrix
    local_stiffness_slices = compute_local_stiffness_slices(
        num_elements,
        polynomial_order
    )
    #
    # Iterate over the elements
    for i in range(num_elements):
        #
        # Construct local stiffness matrices
        K_el = construct_local_stiffness(
            elements[i],
            shape_functions,
            quad_data,
            k,
            b,
            c
        )
        local_slice = local_stiffness_slices[i]
        K_naive[local_slice] += K_el

    #
    # Grab values for penalty enforcement
    K0 = atleast2d(K_naive[:, 0]).transpose()
    KN = atleast2d(K_naive[:, -1]).transpose()

    #
    # Reduce matrix size to account for known data
    K = K_naive
    essential_bcs = has_essential_boundary_condition(alpha)
    neumann_bcs = has_neumann_boundary_condition(beta)
    if essential_bcs[0]:
        K = K[1:,1:]
    elif not neumann_bcs[0]:
        K[0, 0] -= k(0) * beta[0] / alpha[0]

    if essential_bcs[1]:
        K = K[:-1,:-1]
    elif not neumann_bcs[1]:
        K[-1, -1] += k(1) * beta[1] / alpha[1]

    #
    # Incorporate boundary conditions

    return (K, K0, KN)


def construct_load(
    elements,
    polynomial_order,
    alpha,
    beta,
    gamma,
    K0,
    KN,
    f,
    k
):
    #
    #
    quad_data = quad_pts[3]
    #
    # Specify the basis functions
    shape_functions = shape_functions_dict[polynomial_order]
    #
    # First construct the naive load vector
    # which does not account for boundary conditions
    # or reduction necessary to prevent singularity
    num_elements = len(elements)
    num_nodes = num_elements * polynomial_order + 1
    F_naive = zeros((num_nodes, 1))
    #
    # First, identify the locations of the
    # submatrices in the master matrix
    local_load_slices = compute_local_load_slices(
        num_elements,
        polynomial_order
    )

    #
    # Iterate over the elements
    for i in range(len(elements)):
        #
        # Construct local stiffness matrices
        F_el = construct_local_load(
            elements[i],
            shape_functions,
            quad_data,
            f
        )

        local_slice = local_load_slices[i]
        F_naive[local_slice] += F_el

    #
    # Reduce vector size to account for known data
    F = F_naive
    essential_bcs = has_essential_boundary_condition(alpha)
    neumann_bcs = has_neumann_boundary_condition(beta)
    
    REMOVE_FIRST = False
    REMOVE_LAST = False
    if essential_bcs[0]:
        F = F - K0 * gamma[0]
        REMOVE_FIRST = True
    elif neumann_bcs[0]:
        F[0, 0] -= k(0) * gamma[0] / alpha[0]
    else:
        # print("Natural BCs at 0")
        F[0, 0] -= k(0) * gamma[0] / alpha[0]

    if essential_bcs[1]:
        F = F - KN * gamma[1]
        REMOVE_LAST = True
    elif neumann_bcs[1]:
        F[-1, 0] += k(1) * gamma[1] / alpha[1]
    else:
        # print("Natural BCs at 1")
        F[-1, 0] += k(1) * gamma[1] / alpha[1]

    if REMOVE_FIRST: F = F[1:, :]
    if REMOVE_LAST: F = F[:-1, :]

    return F


def compute_local_stiffness_slices(
    num_elements,
    polynomial_order
):
    return [
        s_[
            (i * polynomial_order):((i+1) * polynomial_order + 1),
            (i * polynomial_order):((i+1) * polynomial_order + 1)
        ] for i in range(num_elements)
    ]


def compute_local_load_slices(
    num_elements,
    polynomial_order
):
    return [
        s_[
            (i * polynomial_order):((i+1) * polynomial_order + 1), 0:1
        ] for i in range(num_elements)
    ]


def construct_local_load(
    element,
    shape_functions,
    quad_data,
    f
):
    #
    # Set quad data
    x_quad = quad_data.x
    w_quad = quad_data.w

    #
    # Init empty matrix
    num_shape_functions = len(shape_functions)
    f_el = zeros((num_shape_functions, 1))

    #
    #
    x0 = element.x[0]
    xl = element.x[-1]

    #
    # Transform the quadrature points
    x_quad_transform = coord_transform(
        x_quad,
        x0,
        xl
    )

    #
    # Evaluate the functions at the quadrature points
    f_quad = f(x_quad_transform)
    psi_quad = [p.psi(p, x_quad) for p in shape_functions]

    for i in range(num_shape_functions):
        #
        # Perform quadrature
        f_el[i, 0] = ((xl - x0) / 2.) * np_sum(
            w_quad * (
                f_quad * psi_quad[i]
            )
        )
    #
    return f_el


def construct_local_stiffness(
    element,
    shape_functions,
    quad_data,
    k,
    b,
    c
):
    #
    # Set quad data
    x_quad = quad_data.x
    w_quad = quad_data.w

    #
    # Init empty matrix
    num_shape_functions = len(shape_functions)
    k_el = zeros((num_shape_functions, num_shape_functions))

    #
    #
    x0 = element.x[0]
    xl = element.x[-1]

    #
    # Transform the quadrature points
    x_quad_transform = coord_transform(
        x_quad,
        x0,
        xl
    )

    #
    # Evaluate the functions at the quadrature points
    k_quad = k(x_quad_transform)
    b_quad = b(x_quad_transform)
    # print("x_trans: ", x_quad_transform)
    # print("b_quad: ", b_quad)
    c_quad = c(x_quad_transform)
    psi_quad = [p.psi(p, x_quad) for p in shape_functions]
    # print("PQUAD: ", psi_quad)
    dpsi_quad = [p.dpsi(p, x_quad) for p in shape_functions]

    for index in ndindex(num_shape_functions, num_shape_functions):
        #
        # Get indices in matrix
        i = index[0]
        j = index[1]

        #
        # Perform quadrature
        jac = (xl - x0) / 2.
        k_el[i, j] = jac * np_sum(
            w_quad * (
                b_quad * psi_quad[i] * psi_quad[j] +
                c_quad * psi_quad[i] * dpsi_quad[j] * (1 / jac) +
                k_quad * dpsi_quad[i] * dpsi_quad[j] * (1 / jac)**2
            )
        )
        # print(k_el)

    return k_el

