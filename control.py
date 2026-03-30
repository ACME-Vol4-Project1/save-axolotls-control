# TODO: add control utilities

import numpy as np
import sympy as sy

### FOR SOLVE_BVP ###
def define_variables():
    """Define sympy variables for use in other functions."""
    u1, u2 = sy.symbols("u1, u2")
    S, I, T1, D = sy.symbols("S, I, T1, D")
    u = sy.Matrix([u1, u2])
    x = sy.Matrix([S, I, T1, D])
    return u, x

def cost_int(x, u, Q=sy.Matrix([[0, 0, 0, 0],
                               [0, 0, 0, 0],
                               [0, 0, 0, 0],
                               [0, 0, 0, 1]]), R=sy.Identity(2)):
    """Inside of the cost integral, x.TQx + u.TRu."""
    return x.T@Q@x + u.T@R@u

def hamiltonian(f, L, lam):
    """Define the hamiltonian for a system defined by dx = f, J[u] = int_0^inf L"""
    return lam @ f - L