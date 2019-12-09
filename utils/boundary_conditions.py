from numpy import zeros

def has_essential_boundary_condition(
    alpha
):
    return [a == 0 for a in alpha]


def has_neumann_boundary_condition(
    beta
):
    return [a == 0 for a in beta]

def set_essential_boundary_conditions(
    alpha,
    gamma,
    uu
):
    n = len(uu)

    has_essential_bcs = has_essential_boundary_condition(alpha)

    if not (has_essential_bcs[0] or has_essential_bcs[1]):
        return uu

    if has_essential_bcs[0]:
        n += 1

    if has_essential_bcs[1]:
        n += 1

    u = zeros((n,))

    if has_essential_bcs[0] and has_essential_bcs[1]:
        u[1:-1] = uu
        u[0] = gamma[0]
        u[-1] = gamma[1]
    elif has_essential_bcs[0]:
        u[1:] = uu
        u[0] = gamma[0]
    else:
        u[:-1] = uu
        u[-1] = gamma[1]

    print("u: ", u.shape)
    return u