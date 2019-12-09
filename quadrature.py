from math import sqrt
from numpy import linspace
from utils.constants import quad_pts

#
# General integration
def gauss_int(f=None, num_pts=None, a=None, b=None):
    
    #
    # Get points and weights
    Q = quad_pts[num_pts]
    x = Q.x
    w = Q.w
    
    #
    # Approximate the integral
    val = 0.
    c1 = (b - a) / 2.
    c0 = (b + a) / 2.

    #
    # Evaluate at each quadrature points and sum
    for i in range( len( x ) ):
        wi = w[i]
        xi = x[i]

        val = f(c1 * xi + c0) * wi + val
    #
    return val * c1


def coord_transform(
    xi,
    a,
    b
):
    return a + 0.5 * (b - a) * (1 + xi)
