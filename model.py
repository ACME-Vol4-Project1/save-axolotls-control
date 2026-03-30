import numpy as np
import sympy as sy
from naming import *
from params import default_params
from copy import deepcopy

### Model Class
# initialize a model with the desired attributes, then we can have class functions with the jacobian, other stuff as functions...
class Model():
    def __init__(self, seasonal=True, params=None):
        """An instance of this class provides numeric functions relevant to solving and
        optimizing the model, constructed from the symbolic model representation and with 
        default or overridden parameters.

        Arguments
            seasonal=True: specify if you want the seasonal variant of our model
            params=None or dictionary: use if you want to override any or all of the default parameters.

        Returns
            Model object, with the following functions:
                Model.f(t, x, u): compute x' = f(t, x, u) for use in solve_ivp and other places
                Model.f_x(t, x, u): compute the jacobian D_xf(t, x, u) for use in LQR
                Model.f_u(t, x, u): compute the jacobian D_uf(t, x, u) for use in LQR
        """
        # get the right params dict, overriding defaults with input dict
        temp_params = default_params(seasonal=seasonal)
        if params is not None:
            temp_params.update(params)
        params = temp_params
        self.params = params

        # construct sympy expression + fill parameters
        sy_f = sy_f_full(seasonal=seasonal)
        self.sy_f = sy_f.subs(params)

        # now that parameters are filled in, we need to define the desired model functions
        # get standard names
        t = sy_vars_temporal()
        x = sy_vars_model()
        u = sy_vars_control()
        args = [t, x, u]

        # compute f and the jacobian of f, separating cases where
        # 1. u1 is nonzero (the sauna transfer term should be nonzero unless at equilibrium)
        # 2. T1 and u1 are both zero (no saunas utilized, so the sauna transfer term should be zero)
        # 3. u1 is zero but T1 is nonzero (raise exception)

        # this situation may produce unstable results, to analyze later.

        # Cases 1: u1 is nonzero
        sy_f_x = sy_f.jacobian(x)
        sy_f_u = sy_f.jacobian(u)

        self._f_case1 = sy.lambdify(
            args,
            sy_f
        )

        self._f_x_case1 = sy.lambdify(
            args,
            sy_f_x,
        )

        self._f_u_case1 = sy.lambdify(
            args,
            sy_f_u,
        )

        # Cases 2: T1 and u1 are both zero
        params_case2 = deepcopy(params)
        params_case2.update({"ζ": 0})
        sy_f_case2 = sy_f.subs(params_case2)

        sy_f_x_case2 = sy_f_case2.jacobian(x)
        sy_f_u_case2 = sy_f_case2.jacobian(u)

        self._f_case2 = sy.lambdify(
            args,
            sy_f_case2
        )

        self._f_x_case2 = sy.lambdify(
            args,
            sy_f_x_case2,
        )

        self._f_u_case2 = sy.lambdify(
            args,
            sy_f_u_case2,
        )

    def f(self, t, x, u):
        """Compute numerical jacobian of f with respect to x"""
        if not u[0] == 0:
            return self._f_case1(t, x, u).ravel() # ravel so it's not a 2d array
        
        elif u[0] == 0 and x[2] == 0:
            return self._f_case2(t, x, u).ravel() # ravel so it's not a 2d array
        
        else:
            raise ValueError("f is nan when $u_1 = 0$ and $T_1 \\neq 0$")

    def f_x(self, t, x, u):
        """Compute numerical jacobian of f with respect to x"""
        if not u[0] == 0:
            return self._f_x_case1(t, x, u)
        
        elif u[0] == 0 and x[2] == 0:
            return self._f_x_case2(t, x, u)
        
        else:
            raise ValueError("f is nan when $u_1 = 0$ and $T_1 \\neq 0$")

    def f_u(self, t, x, u):
        """Compute numerical jacobian of f with respect to u"""
        if not u[0] == 0:
            return self._f_u_case1(t, x, u)
        
        elif u[0] == 0 and x[2] == 0:
            return self._f_u_case2(t, x, u)
        
        else:
            raise ValueError("f is nan when $u_1 = 0$ and $T_1 \\neq 0$")



### NUMERICAL MODEL IMPLEMENTATION
# define the model
# def model_equation(t, y, theta_func, control1, control2, alpha_func, gamma, zeta, beta_func):
#     """ODE equation for our Chytrid fungus model. The y contains the values of
#     each compartment at time t, where the compartments are
#         - S: susceptible
#         - I: infected
#         - T1: treatment 1 (frog sauna)
#         - T2: treatment 2 (antifungal bath)
#         - D: deceased

#         Parameters:
#         - theta: season. 0 mod 2pi indicates peak summer, pi mod 2pi indicates peak winter
#         - control1: control parameters u1 and u2 (functions of t)
#         - control2: functions of u1 and u2
#         - other parameters: express movement between components
#     """
#     S, I, T1, D = y
#     u1, u2 = control1
#     K_func, eta_func, xi_func = control2
#     K = K_func(u1)
#     eta = eta_func(u2)
#     xi = xi_func(u2, I)

#     # apply season
#     theta = theta_func(t)
#     alpha = alpha_func(theta)
#     beta = beta_func(theta)

#     # precompute T1 transition term
#     enter_habitat = zeta*(1 - T1/K) if K > 0 else 0

#     # define transition equations
#     dS = -alpha*I*S + gamma*T1 + eta*xi*I - S * enter_habitat
#     dI = alpha*I*S - I * enter_habitat - eta*xi*I - beta*I
#     dT1 = (I + S) * enter_habitat - gamma*T1
#     dD = beta*I

#     return np.array([dS, dI, dT1, dD])


### SYMPY MODEL DEFINITION
# define model
def sy_f_full(seasonal=True):
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
    fsimp = sy_f_simple()
    f = fsimp.subs(dict(zip(names, tosub)))

    return f

def sy_f_simple():
    """Uses the sympy definitions below to construct a sympy expression for f(x)
    where x' = f(x) is our ODE. **This expression should only be used to construct
    a dynamic expression if any of the parameters are dynamic, and to double check the equations
    look right**. Can be used to construct both seasonal and autonomous versions of the ODE

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

def sy_params_dynamic(seasonal=True):
    """Constructs the dynamic variables in terms of t, state, u1, and u2.
    CHECK sy_params_dynamic_names and ensure these two match!!"""
    t = sy_vars_temporal()
    S, I, T1, D = sy_vars_model()
    u1, u2 = sy_vars_control()

    ### construct non-seasonal parameters
    # construct xi: depends on u2 and I
    xi = u2 / ((I + 50) * 100)

    # construct K: depends on u1
    K = 365 * u1 / 2

    ### construct seasonal parameters
    if seasonal:
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
    
    else:
        # theta is a dummy variable at this point that we don't need if 
        # we're not building a seasonal model, so we can just replace it 
        # with its own name
        # honestly, also the same with alpha and beta seeing as they don't
        # depend on the controls.
        alpha, beta, theta = sy_params_dynamic_names()[:3]

    return [alpha, beta, theta, xi, K]
