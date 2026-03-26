import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import solve_ivp

# define the model
def model_equation(t, y, control1, control2, alpha, gamma, zeta, beta, delta, mu):
    """ODE equation for our Chytrid fungus model. The y contains the values of
    each compartment at time t, where the compartments are
        - S: susceptible
        - I: infected
        - T1: treatment 1 (frog sauna)
        - T2: treatment 2 (antifungal bath)
        - R: deceased but can still transmit the fungus

        Parameters:
        - control1: control parameters u1 and u2 (functions of t)
        - control2: functions of u1 and u2
        - other parameters: express movement between components
    """
    S, I, T1, T2, R = y
    u1, u2 = control1
    K, eta, xi = control2

    # define transition equations
    dS = -alpha*I*S + gamma*T1 + eta(u2)*T2 - zeta*S*(1 - (S + I)/K(u1))
    dI = alpha*I*S - zeta*I*(1 - (S + I)/K(u1)) - xi(u2)*I - beta*I
    dT1 = zeta*(I + S)*(1 - (S + I)/K(u1)) - gamma*T1
    dT2 = xi(u2)*I - eta(u2)*T2 - (1 - eta(u2))*T2
    dR = beta*I + (1 - eta(u2))*T2

    return [dS, dI, dT1, dT2, dR]