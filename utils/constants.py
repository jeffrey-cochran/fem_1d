from collections import namedtuple
from numpy import sqrt, asarray
from enum import IntEnum

#
# Pair quadrature points and weights in a data structure
QuadData = namedtuple('QuadData', 'x w')


#
# Create enums for polynomial order
class POLYNOMIAL(IntEnum):
    LINEAR = 1
    QUADRATIC = 2
    CUBIC = 3


quad_pts = {
    1: QuadData(
        asarray([0.0]),
        asarray([2.])
    ),
    2: QuadData(
        asarray([-1 / sqrt(3), 1 / sqrt(3)]),
        asarray([1., 1.])
    ),
    3: QuadData(
        asarray([-sqrt(3 / 5.), 0.0, sqrt(3 / 5.)]),
        asarray([5/9., 8/9., 5/9.])
    )
}

#
# Penalty factor for essential conditions
penalty_factor = 1e10