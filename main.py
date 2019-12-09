from utils.element_def import element
from utils.constants import POLYNOMIAL
from assemble import construct_load, construct_stiffness
from numpy import vectorize
from numpy import linspace
nelems = 2
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
b = lambda x: 1.
b = vectorize(b)
c = lambda x: x
c = vectorize(c)
# f = lambda x: 0.25 * pi**2 * cos(pi * x / 2.)
f = lambda x: 3.*x**2 - 5.
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


F = construct_load(
    elements,
    polynomial_order,
    alpha,
    beta,
    gamma,
    K0,
    KN,
    f
)

x = linspace(0, 1, 100)
# t = lambda x: -(5 * x**2 / 4. - 25 * x / 24.) if (x <= 2/3.) else -0.5*(x**2 - 5 * x / 6. - 1 / 6.)
# t = lambda x: cos(pi * x / 2.)
t = lambda x: (1 + x)**2
t = vectorize(t)
y = t(x)

print("COND: ", cond(K))
xx = linspace(0, 1, nelems * polynomial_order + 1)
uu = np_solve(K, F).flatten()
u = set_essential_boundary_conditions(alpha, gamma, uu)

plt.plot(x, y, 'k-', xx, u, 'ro')
plt.show()