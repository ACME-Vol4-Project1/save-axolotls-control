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
