# TODO: add control utilities

import numpy as np
import sympy as sy
from naming import *                            # For the model definitions
from model import sy_params_dynamic, sy_f_simple
# from model import sy_f_full as f              # Defining them both here for the time being

### FOR SOLVE_BVP ###
def define_variables(no_dead=True):
    """Define sympy variables for use in other functions."""
    u1, u2 = sy.symbols("u1, u2")
    S, I, T1, D = sy.symbols("S, I, T1, D")
    u = sy.Matrix([u1, u2])

    if no_dead:
        x = sy.Matrix([S, I, T1])
        lam = sy.symbols("l1, l2, l3")
    else:
        x = sy.Matrix([S, I, T1, D])
        lam = sy.symbols("l1, l2, l3, l4")
    return x, u, sy.Matrix([lam])

def L(x, u, no_dead=True):
    """Inside of the cost integral, x.TQx + u.TRu."""
    R = sy.Identity(2)

    if no_dead:
        Q = sy.Matrix([[0, 0, 0],
                       [0, 3, 0], 
                       [0, 0, 0]])
        
    else:
        Q=sy.Matrix([[0, 0, 0, 0],
                    [0, 0, 0, 0],
                    [0, 0, 0, 0],
                    [0, 0, 0, 1]])
    return x.T@Q@x + u.T@R@u

def H(f, L, lam):
    """Define the hamiltonian for a system defined by dx = f, J[u] = int_0^inf L"""
    return lam @ f - L

def hamiltionian_partials(H, x, lam, u):
    ''' Takes in Sympy expressions for H, x, lam, 
    returns the state and costate eqns as well as Dh/Du
    as vectors of sympy expressions'''

    # make them each Sympy matrices
    x = sp.Matrix(x)
    lam = sp.Matrix(lam)

    # sympy can only do one derivative at a time, stack them into vectors
    x_dot = sp.Matrix([sp.diff(H, l) for l in lam])     
    lam_dot = -sp.Matrix([sp.diff(H, xi) for xi in x])

    # Can't remember if we'll use Dh/Du = 0
    zero = sp.Matrix([sp.diff(H, ui) for ui in u])

    return x_dot, lam_dot, zero 



### SYMPY MODEL DEFINITION
# define model
def sy_f_full(seasonal=True, no_dead=True):
    """This is the method to generate a full formulation of our model.

    - non-seasonal: autonomous, so we can use infinite horizon LQR
    - seasonal: non-autonomous, so infinite horizon LQR is no good
    
    Substitutes the dynamic variables into sy_f_simple. This can be used to actually
    evaluate when parameter values are known, but is much more painful to visualize than
    sy_f_simple."""
    # get dynamic parameter names and values
    names = sy_params_dynamic_names()
    tosub = sy_params_dynamic(seasonal=seasonal)

    # get simple f and sub in params
    if no_dead:
        fsimp = sy_f_simple_no_dead()
    else:
        fsimp = sy_f_simple()
    f = fsimp.subs(dict(zip(names, tosub)))

    return f

def sy_f_simple_no_dead():
    """Uses the sympy definitions below to construct a sympy expression for f(x)
    where x' = f(x) is our ODE. **This expression should only be used to construct
    a dynamic expression if any of the parameters are dynamic, and to double check the equations
    look right**. Can be used to construct both seasonal and autonomous versions of the ODE

    None of the sympy variables are functions in the output - and some of the parameters
    depend on t, u1, u2, and the state, so if you want to take derivatives or evaluate, sub in the dynamic
    expressions first.
    """
    # Removed the D class from the model
    S, I, T1, _ = sy_vars_model()
    alpha, beta, theta, xi, K = sy_params_dynamic_names()
    gamma, zeta, eta = sy_params_static()
    f = sy.Matrix([                 # f = [S', I', T1', D']^T
        [- alpha * I * S + gamma * T1 + eta * xi * I - zeta * S * (1 - T1 / K)],
        [alpha * I * S - zeta * I * (1 - T1 / K) - eta * xi * I - beta * I],
        [zeta * (I + S) * (1 - T1 / K) - gamma * T1]
    ])
    
    return f