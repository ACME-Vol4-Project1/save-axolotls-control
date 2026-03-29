import numpy as np
import sympy as sy
from matplotlib import pyplot as plt
from scipy.integrate import solve_ivp

# define the model
def model_equation(t, y, theta_func, control1, control2, alpha_func, gamma, zeta, beta_func):
    """ODE equation for our Chytrid fungus model. The y contains the values of
    each compartment at time t, where the compartments are
        - S: susceptible
        - I: infected
        - T1: treatment 1 (frog sauna)
        - T2: treatment 2 (antifungal bath)
        - D: deceased

        Parameters:
        - theta: season. 0 mod 2pi indicates peak summer, pi mod 2pi indicates peak winter
        - control1: control parameters u1 and u2 (functions of t)
        - control2: functions of u1 and u2
        - other parameters: express movement between components
    """
    S, I, T1, D = y
    u1, u2 = control1
    K_func, eta_func, xi_func = control2
    K = K_func(u1)
    eta = eta_func(u2)
    xi = xi_func(u2, I)

    # apply season
    theta = theta_func(t)
    alpha = alpha_func(theta)
    beta = beta_func(theta)

    # precompute T1 transition term
    enter_habitat = zeta*(1 - T1/K) if K > 0 else 0

    # define transition equations
    dS = -alpha*I*S + gamma*T1 + eta*xi*I - S * enter_habitat
    dI = alpha*I*S - I * enter_habitat - eta*xi*I - beta*I
    dT1 = (I + S) * enter_habitat - gamma*T1
    dD = beta*I

    return np.array([dS, dI, dT1, dD])

def model_parameters():
    alpha, gamma, zeta, beta = .00001, 0.1, 0.15, 0.01 # central time average parameter estimates
    control2 = K, eta, xi # other parameters

    return 

# SYMPY CODE - DEFINES ODES + MODEL
# updated to model 1.2

# progress towards lambdifying sympy models
# # lambdify sympy models
# class symodel():
#     def __init__(self):
#         self.f = lambdify_f()

# def lambdify_f():
#     """Lambdifies sympy expression for the model"""
#     t = sy_vars_temporal()
#     S, I, T1, D = sy_vars_model()
#     u1, u2 = sy_vars_control()
#     alpha, beta, theta, xi, K = sy_params_dynamic()
#     gamma, zeta, eta = sy_params_static()

#     # syf = sympy_f()

#     # generate a function that will evaluate the dynamic params

#     # generate a function that we can pass in everything to raw
#     flong = sy.lambdify([
#         t,                              # t
#         [S, I, T1, D],                  # y
#         [u1, u2],                       # u
#         [alpha, beta, theta, xi, K],    # dynamic params
#         [gamma, zeta, eta]              # static params
#     ],
#     syf)

#     # then generate a wrapper function that will compute the dynamic variables first



# define model
def sy_f_full():
    """Substitutes the dynamic variables into sy_f_simple. This can be used to actually
    evaluate when parameter values are known, but is much more painful to visualize."""
    # get dynamic parameter names and values
    names = sy_params_dynamic_names()
    tosub = sy_params_dynamic()

    # get simple f and sub in params
    fsimp = sy_f_simple()
    f = fsimp.subs(dict(zip(names, tosub)))

    return f

def sy_f_simple():
    """Uses the sympy definitions below to construct a sympy expression for f(x)
    where x' = f(x) is our ODE. **This expression should only be used to construct
    a dynamic expression if any of the parameters are dynamic, and to double check the equations
    look right**.

    None of the sympy variables are functions in the output - and some of the parameters
    depend on t, u1, u2, and the state, so if you want to take derivatives or evaluate, sub in the dynamic
    expressions first.
    """
    S, I, T1, D = sy_vars_model()
    alpha, beta, theta, xi, K = sy_params_dynamic_names()
    gamma, zeta, eta = sy_params_static()
    f = sy.Matrix([                 # f = [S', I', T1', D']^T
        [- alpha * I * S + gamma * T1 + eta * xi * I - zeta * S * (1 - T1 / K)],
        [alpha * I * S - zeta * I * (1 - T1 / K) - eta * xi * I - beta * I],
        [zeta * (I + S) * (1 - T1 / K) - gamma * T1],
        [beta * I]
    ])
    return f

def sy_params_dynamic():
    """Constructs the dynamic variables in terms of t, u1, and u2.
    CHECK sy_params_dynamic_names and ensure these two match!!"""
    t = sy_vars_temporal()
    S, I, T1, D = sy_vars_model()
    u1, u2 = sy_vars_control()

    # construct theta: the time of year.
    # 0 mod 2pi = summer solstice.
    # pi mod 2pi = winter solstice.
    # assume that t=0 is in the summer and years are 365 days
    theta = 2 * sy.pi * t / 365

    # helper function: generate sinusoidal periodic function
    sinu = lambda base, omega, phi: base * (1 + omega + sy.cos(theta + phi))

    # construct alpha: periodic sinusoidal
    alpha = sinu(*sy.symbols("α_0,ω_α,φ_α"))

    # construct beta: periodic sinusoidal
    beta = sinu(*sy.symbols("β_0,ω_β,φ_β"))

    # construct xi: depends on u2 and I
    xi = u2 / ((I + 50) * 100)

    # construct K: depends on u1
    K = 365 * u1 / 2

    return [alpha, beta, theta, xi, K]

# define parameter names
def sy_params_static():
    return sy.symbols("γ,ζ,η")

def sy_params_dynamic_names():
    return sy.symbols("α,β,θ,ξ,K")

# define variables
def sy_vars_control():
    return sy.symbols("u1,u2")

def sy_vars_model():
    return sy.symbols("S,I,T1,D")

def sy_vars_temporal():
    return sy.symbols("t")
