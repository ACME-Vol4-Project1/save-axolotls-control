import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import solve_ivp

# define the model
def model_equation(t, y, params=[.09, .005, .001, .004, .001, .005, .001]):
    """ODE equation for our Chytrid fungus model. The y contains the values of
    each compartment at time t, where the compartments are
        - S: susceptible
        - I: infected
        - T1: treatment 1 (frog sauna)
        - T2: treatment 2 (antifungal bath)
        - R1: deceased but can still transmit the fungus
        - R2: host and fungus deceased
    """
    S, I, T1, T2, R1, R2 = y
    alpha, gamma, eta, zeta, beta, delta, xi = params

    # # introduce periodic seasonal transmission
    # sigma = 6    # amplitude of seasonal variation in transmission
    # alpha = alpha1*(1 + sigma*np.cos(2*np.pi*t))
    
    # define transition equations
    dS = -alpha*I*S - gamma*I*R1 + eta*T2 - zeta*S*(sum(y) - (T1 + T2 + R1 + R2))
    dI = alpha*I*S + gamma*I*R1 - zeta*I*(sum(y) - (T1 + T2 + R1 + R2)) - xi*I - beta*I
    dT1 = zeta*(I + S)*(sum(y) - (T1 + T2 + R1 + R2)) 
    dT2 = xi*I - eta*T2
    dR1 = beta*I - delta*R1                            
    dR2 = delta*R1                                      

    return [dS, dI, dT1, dT2, dR1, dR2]