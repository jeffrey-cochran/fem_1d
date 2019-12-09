from .shape_function_def import shape_function_polynomial as sfp
from .constants import POLYNOMIAL

linear_shape_functions = [
    sfp(order=POLYNOMIAL.LINEAR, basis_number=i) for i in range(2)
]

quadratic_shape_functions = [
    sfp(order=POLYNOMIAL.QUADRATIC, basis_number=i) for i in range(3)
]

cubic_shape_functions = [
    sfp(order=POLYNOMIAL.CUBIC, basis_number=i) for i in range(4)
]

shape_functions_dict = {
    POLYNOMIAL.LINEAR: linear_shape_functions,
    POLYNOMIAL.QUADRATIC: quadratic_shape_functions,
    POLYNOMIAL.CUBIC: cubic_shape_functions
}