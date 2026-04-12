import sympy as sy

# define parameter names
def sy_params_static():
    """Standardized parameter names: γ,ζ,η"""
    return sy.symbols("γ,ζ,η")

def sy_params_dynamic_names():
    """Standardized parameter names: α,β,θ,ξ,K"""
    return sy.symbols("α,β,θ,ξ,K")

# define variables
def sy_vars_control():
    """Standardized parameter names: u1, u2"""
    return sy.symbols("u1,u2")

def sy_vars_model():
    """Standardized parameter names: S,I,T1,D"""
    return sy.symbols("S,I,T1,D")

def sy_vars_temporal():
    """Standardized parameter names: t"""
    return sy.symbols("t")


### modifications for reintroduction model (same temporal and control variables)
def sy_params_static2():
    """Define parameters that don't depend on control or season
        - d2: natural death rate of mature frogs
        - gamma: rate of frogs leaving saunas
        - zeta: rate of frogs entering saunas
        - phi: rate of juvenile frogs surviving to maturity
        - a, b: scaling factor for juvenile frogs death rate and treatment rate
    """
    return sy.symbols("d2, gamma, zeta, phi, a, b")

def sy_params_dynamic_names2():
    """Parameters that depend on season or control.
        - eps: spontaneous infection rate (from water or other species)
        - d1: birth rate
        - a1: infection rate of juvenile frogs
        - a2: infection rate of mature frogs
        - beta: death rate from infection
        - theta: helps define seasonality
        - nu: rate of introduction of vaccinated frogs
        - K: carrying capacity of frog saunas
    """
    return sy.symbols("eps, d1, a1, a2, beta, theta, nu, K")

def sy_vars_model2():
    """Model compartment names
        - S1: juvenile frogs
        - S2: mature frogs (susceptible)
        - V: vaccinated
        - I: infected
        - T: frog sauna treatment
        - D: dead
    """
    return sy.symbols("S1, S2, V, I, T, D")