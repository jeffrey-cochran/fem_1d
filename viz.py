# from func_lib import linear, quadratic
from numpy import linspace, zeros
from basis_function_def import basis_function_polynomial

def visualize_solution(
    ax,
    x,
    w,
    poly_order,
    num_elements
):
    master_xplot = zeros(100*num_elements+1)
    master_fplot = zeros(100*num_elements+1)
    for i in range(0, len(x), poly_order):
        if i == len(x) - 1:
            break

        xe = x[i:i+poly_order+1]
        we = w[i:i+poly_order+1].tolist()

        x_plot = linspace(xe[0], xe[-1], 101)
        f_plot = zeros(x_plot.shape)

        for j in range(len(we)):
            cur_f = basis_function_polynomial(
                left_boundary=xe[0],
                spacing=(xe[-1]-xe[0]),
                order=poly_order,
                basis_number=j+1
            )
            f_plot = we[j] * cur_f(x_plot) + f_plot
            
        #
        master_xplot[i*100:(i+1)*100+1]=x_plot
        master_fplot[i*100:(i+1)*100+1]=f_plot
    #
    return (master_xplot, master_fplot)