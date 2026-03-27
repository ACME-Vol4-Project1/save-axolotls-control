import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import solve_ivp

# define the model
def model_equation(t, y, control1, control2, alpha, gamma, zeta, beta):
    """ODE equation for our Chytrid fungus model. The y contains the values of
    each compartment at time t, where the compartments are
        - S: susceptible
        - I: infected
        - T1: treatment 1 (frog sauna)
        - T2: treatment 2 (antifungal bath)
        - D: deceased

        Parameters:
        - control1: control parameters u1 and u2 (functions of t)
        - control2: functions of u1 and u2
        - other parameters: express movement between components
    """
    S, I, T1, T2, D = y
    u1, u2 = control1
    K_func, eta_func, xi_func = control2
    K = K_func(u1)
    eta = eta_func(u2)
    xi = xi_func(u2, I)

    # precompute T1 transition term
    enter_habitat = zeta*(1 - T2/K) if K > 0 else 0

    # define transition equations
    dS = -alpha*I*S + gamma*T1 + eta*T2 - S * enter_habitat
    dI = alpha*I*S - I * enter_habitat - xi*I - beta*I
    dT1 = (I + S) * enter_habitat - gamma*T1
    dT2 = xi*I - eta*T2 - (1 - eta)*T2
    dD = beta*I + (1 - eta)*T2

    return np.array([dS, dI, dT1, dT2, dD])