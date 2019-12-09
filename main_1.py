from utils.element_def import element
from utils.constants import POLYNOMIAL
from assemble import construct_load, construct_stiffness
from utils.boundary_conditions import set_essential_boundary_conditions
from numpy import vectorize
from numpy import linspace
from numpy import savetxt
from numpy import exp as np_exp
from numpy.linalg import solve as np_solve
from math import exp
import matplotlib.pyplot as plt

nelems = 16
polynomial_order = POLYNOMIAL.LINEAR
x = linspace(0, 1, nelems + 1)
elements = [
    element(x[i], x[i + 1], polynomial_order) for i in range(nelems)
]
alpha = [0., 0.]
beta = [1., 1.]
gamma = [1., 4.]
k = lambda x: 2. + x
k = vectorize(k)
b = lambda x: 0.
b = vectorize(b)
c = lambda x: 0.
c = vectorize(c)
f = lambda x: -6 - 4*x
f = vectorize(f)

K, K0, KN = construct_stiffness(
    elements,
    polynomial_order,
    alpha,
    beta,
    gamma,
    k,
    b,
    c
)

savetxt("stiffness_quad.csv", K, delimiter=" & ", fmt="%.2f")

F = construct_load(
    elements,
    polynomial_order,
    alpha,
    beta,
    gamma,
    K0,
    KN,
    f,
    k
)

savetxt("load_quad.csv", F, delimiter=" \\\\ ", fmt="%.2f")

x = linspace(0, 1, 100)
t = lambda x: (1 + x)**2
t = vectorize(t)
y = t(x)

# print("COND: ", cond(K))
xx = linspace(0, 1, nelems * polynomial_order + 1)
uu = np_solve(K, F).flatten()
u = set_essential_boundary_conditions(alpha, gamma, uu)

plt.cla()
plt.plot(x, y, 'k-', label='True Sol.')
plt.plot(xx, u, 'ro', label='Numerical Sol.')
plt.legend()
plt.show()