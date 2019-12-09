from numpy import exp, linspace
from numpy.linalg import solve as np_solve
from assemble import construct_stiffness
from assemble import construct_load
from utils.constants import POLYNOMIAL

def solve(
    nelem,
    h,
    deg,
    alpha,
    beta,
    gamma,
    k=None,
    b=None,
    c=None,
    f=None
):
    #
    # Cast polynomial degree to enums
    deg = POLYNOMIAL(deg)
    #
    # Check inputs make sense
    assert( nelem > 0 )
    assert( nelem == len(h) )
    assert( deg in [POLYNOMIAL.LINEAR, POLYNOMIAL.QUADRATIC] )
    assert( len(alpha) == 2 )
    assert( len(beta) == 2 )
    assert( len(gamma) == 2 )

    #
    # Make sure functions are callable on ndarrays
    #
    # NOTE: functions passed as None will be treated
    # as a constant equal to zero. 
    x_test = linspace(0,1,3)
    callable_functions = [k, b, c, f]
    for fun in callable_functions:
        if fun is not None:
            try:
                fun(x_test)
            except:
                try: 
                    fun = lambda x: float(fun)
                except:
                    raise(Exception("One of the functions k, b, c, f is neither callable on ndarrays nor a number. Aborting!"))
        else:
            fun = lambda x: 0.
    
    #
    # Assemble the system
    K = construct_stiffness()
    F = construct_load()

    #
    # Solve the system
    u = np_solve(K, F)

    #
    # Solve auxiliary equations
    if False:
        print("Hello")
    
    #
    return u

# nelem = 8
# x = linspace(0, 1, nelem + 1)
# h = x[1:] - x[:-1]
# deg = 2
# alpha = [1., 0.]
# beta = [-2., 1.]
# gamma = [-1., exp(1.)]
# k = lambda x: 1 + x
# b = lambda x: x**2
# c = lambda x: -x**2 + x
# f = lambda x: -2 * exp(x)

nelem = 8
x = linspace(0, 1, nelem + 1)
h = x[1:] - x[:-1]
deg = 2
alpha = [0., 0.]
beta = [1., 1.]
gamma = [0., 0.]
k = lambda x: 1.

u = solve(
    nelem,
    h,
    deg,
    alpha,
    beta,
    gamma,
    k=k
)

print(u)