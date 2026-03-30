import sympy as sp

def hamiltionian_partials(H, x, lam, u):
    ''' Takes in Sympy expressions for H, x, lam, 
    returns the state and costate eqns as well as Dh/Du
    as vectors of sympy expressions'''

    # make them each Sympy matrices
    x = xp.Matrix(x)
    lam = sp.Matrix(lam)

    # sympy can only do one derivative at a time, stack them into vectors
    x_dot = sp.Matrix([sp.diff(H, l) for l in lam])     
    lam_dot = -sp.Matrix([sp.diff(H, xi) for xi in x])

    # Can't remember if we'll use Dh/Du = 0
    zero = sp.Matrix([sp.diff(H, ui) for ui in u])

    return x_dot, lam_dot, zero 