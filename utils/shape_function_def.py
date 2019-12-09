from numpy import zeros, asarray, sqrt
from numpy.linalg import solve
from .constants import POLYNOMIAL
import scipy

class shape_function_polynomial(object):

    def __init__(
        self,
        order=POLYNOMIAL.CUBIC,
        basis_number=0
    ):
        #
        # Make sure basis function exists for polynomial order
        assert(basis_number >= 0 and basis_number <= int(order) + 1)
        
        #
        # Set evaluation function based on
        # poly order and basis number
        self.basis_number = basis_number
        if order == POLYNOMIAL.LINEAR:
            self.psi = self.linear_psi()
            self.dpsi = self.linear_dpsi()
        elif order == POLYNOMIAL.QUADRATIC:
            self.psi = self.quadratic_psi()
            self.dpsi = self.quadratic_dpsi()
        else:
            self.psi = self.cubic_psi()
            self.dpsi = self.cubic_dpsi()
        
        #
        return

    def psi(self, x): 
        return

    def dpsi(self, x): 
        return

    def linear_psi(
        self
    ):
        #
        #
        if self.basis_number == 0:
            return lambda self, x: .5 * (1 - x)

        #
        return lambda self, x: .5 * (1 + x)

    def quadratic_psi(
        self
    ):
        #
        #
        if self.basis_number == 0:
            return lambda self, x: x * (x - 1.) * .5
        elif self.basis_number == 1:
            return lambda self, x: 1 - x**2

        #
        return lambda self, x: x * (x + 1.) * .5

    def cubic_psi(
        self
    ):
        #
        #
        if self.basis_number == 0:
            return lambda self, x: 9/16. * (1/9.  - x**2) * (x - 1)
        elif self.basis_number == 1:
            return lambda self, x: 27/16. * (1 - x**2) * (1/3. - x)
        elif self.basis_number == 2:
            return lambda self, x: 27/16. * (1 - x**2) * (1/3. + x)

        #
        return lambda self, x: -9/16. * (1/9. - x**2) * (1 + x)

    def linear_dpsi(
        self
    ):
        #
        # 
        if self.basis_number ==0:
            dpsi = lambda self, x: -.5
        else:
            dpsi = lambda self, x: .5
        
        #
        return dpsi

    def quadratic_dpsi(
        self
    ):
        #
        #
        if self.basis_number == 0:
            return lambda self, x: x - .5
        elif self.basis_number == 1:
            return lambda self, x: -2 * x

        #
        return lambda self, x: x + .5
    
    def cubic_dpsi(
        self
    ):
        #
        #
        if self.basis_number == 0:
            return lambda self, x: -9/16. * (3 * x**2 - 2 * x - 1/9.)
        elif self.basis_number == 1:
            return lambda self, x: 27/16. * (3 * x**2 - 2/3. * x - 1)
        elif self.basis_number == 2:
            return lambda self, x: 27/16. * (-3 * x**2 - 2/3. * x + 1)

        #
        return lambda self, x: -9/16. * (-3 * x**2 - 2 * x + 1/9.)
        