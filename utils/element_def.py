from numpy import linspace


class element(object):

    def __init__(
        self,
        xl,
        xr,
        polynomial_order
    ):
        self.x = linspace(xl, xr, int(polynomial_order) + 1)
        return
