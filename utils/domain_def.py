from .element_def import element
from numpy import linspace
from .constants import POLYNOMIAL

class domain(object):

    def __init__(
        self,
        x0,
        xl,
        number_elements,
        polynomial_order
    ):
        #
        # Make sure domain is sensible
        assert( xl > x0)

        #
        # Make sure num elmeents is sensible
        assert( number_elements > 0 )

        #
        # Make sure polynomial order makes sense
        assert( polynomial_order in [ POLYNOMIAL.LINEAR, POLYNOMIAL.QUADRATIC, POLYNOMIAL.CUBIC ] )

        self.element_boundaries = linspace(x0, xl, number_elements + 1)
        self.elements = [ element( self.element_boundaries[i], self.element_boundaries[i+1], polynomial_order) for i in range(number_elements) ]

        #
        return
