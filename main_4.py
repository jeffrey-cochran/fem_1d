from utils.element_def import element
from utils.constants import POLYNOMIAL
from assemble import construct_load, construct_stiffness
from utils.boundary_conditions import set_essential_boundary_conditions
from numpy import vectorize
from numpy import linspace
from numpy import savetxt
from numpy import exp as np_exp
from numpy import zeros
from numpy import divide
from scipy.interpolate import interp1d
from scipy.integrate import quad
from numpy import cos, pi, sin, sqrt, log
from numpy.linalg import solve as np_solve
from math import exp
import matplotlib.pyplot as plt

elems = [4, 8, 16, 32]
err = []
for elem in elems:
    nelems = elem
    polynomial_order = POLYNOMIAL.LINEAR
    x = linspace(0, 1, nelems + 1)
    elements = [
        element(x[i], x[i + 1], polynomial_order) for i in range(nelems)
    ]
    alpha = [0., 1.]
    beta = [1., 0.]
    gamma = [1., -pi / 2.]
    k = lambda x: 1.
    k = vectorize(k)
    b = lambda x: 0.
    b = vectorize(b)
    c = lambda x: 0.
    c = vectorize(c)
    f = lambda x: 0.25 * pi**2 * cos(pi * x / 2.)
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
    t = lambda x: cos(pi * x / 2.)
    tprime = lambda x: -pi / 2. * sin(pi * x / 2.)
    tprime = vectorize(tprime)
    t = vectorize(t)
    y = t(x)

    xx = linspace(0, 1, nelems * polynomial_order + 1)
    uu = np_solve(K, F).flatten()
    u = set_essential_boundary_conditions(alpha, gamma, uu)
    uprime = zeros(xx.shape)
    uprime[:-1] = (u[1:] - u[:-1]) / (xx[1:] - xx[:-1])
    uprime[-1] = tprime(1)
    u = interp1d(xx,u)
    uprime = interp1d(xx, uprime)

    def integrand(zzz):
        try:
            return sqrt((u(zzz) - t(zzz))**2 + (uprime(zzz) - tprime(zzz))**2)
        except:
            print(zzz)
            return 0

    print(x[0], x[-1])
    cur_err = quad(integrand, xx[0], xx[-1])
    err.append(cur_err[0])

print(elems)
print(err)

aa = (log(err[1:]) - log(err[:-1]))
bb = (log(elems[1:]) - log(elems[:-1]))
print(divide(aa, bb))
plt.cla()
plt.plot(log(elems), log(err), 'ro', label='log(H1 error).')
plt.xlabel("log(num elems)")
plt.ylabel("log(H1 err)")
plt.legend()
plt.show()