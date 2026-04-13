import numpy as np
from naming import sy_params_static, sy_params_dynamic_names , sy_vars_control
from naming import sy_params_static2, sy_params_dynamic_names2, sy_vars_model2

### Default Parameters
def default_params(seasonal=True, const_u1=None):
    # setup parameter dictionary for sympy substitution
    params = dict()

    # parameters that don't differ between seasonal and not
    statics = dict(zip(
        sy_params_static(), # γ,ζ,η
        [0.1, 0.15, 1.]
        ))
    params.update(statics)

    # if u_1 is being held constant:
    if const_u1 is not None:
        if const_u1 == 0:
            params["ζ"] = 0 # get rid of carrying capacity term
        else:
            params[sy_vars_control()[0]] = const_u1 # fix u_1

    # if model is autonomous: add alpha and beta
    if not seasonal: 
        dynamics = dict(zip(
            sy_params_dynamic_names()[:2], # alpha, beta
            [0.000035, 0.003]
        ))
        params.update(dynamics)

    # if the model is seasonal: add seasonal alpha and beta
    if seasonal:
        dynamics = dict()
        alphas = dict(zip(
            "α_0,ω_α,φ_α".split(','),
            # [0.000025, 0.5, 1.25 * np.pi]
            [0.000025, 0.5, 1.25 * np.pi] # optimized in test_model.ipynb
        ))
        dynamics.update(alphas)

        betas = dict(zip(
            "β_0,ω_β,φ_β".split(','),
            # [0.0035, 0.86, np.pi]
            [0.0020, 0.86, np.pi] # optimized in test_model.ipynb
        ))
        dynamics.update(betas)

        params.update(dynamics)

    # return parameter dict
    return params

def default_params2(seasonal=True, const_u1=None):
    """Get dictionary of default parameters for substitution in sympy expression"""
    params = dict()

    # static parameters
    statics = dict(zip(
        sy_params_static2(),   # d2, gamma, zeta, phi, a, b
        [0.1, 0.1, 0.15, 1., .25, .25]   # TODO: get acutal params for d2, phi, a, b
    ))
    params.update(statics)

    # check for fixed value of u1
    if const_u1 is not None:
        if const_u1 == 0:
            params["gamma"] = 0
        else:
            params[sy_vars_control()[0]] = const_u1   # fixed u1

    # seasonal parameters
    if seasonal:
        dynamics = dict()

        alphas = dict(zip(
            "a1_0, w_a1, p_a1, a2_0, w_a2, p_a2".split(','),
            [0.000025, 0.5, 1.25*np.pi, 0.000025, 0.5, 1.25*np.pi]    # TODO: get actual parameters a1
        ))
        dynamics.update(alphas)
        
        betas = dict(zip(
            "b_0, w_b, p_b".split(','),
            [0.0020, 0.86, np.pi]    # TODO: get actual parameters
        ))
        dynamics.update(betas)

        deltas = dict(zip(
            "d_0, w_d, p_d".split(','),
            [0.0020, 0.86, np.pi]    # TODO: get actual parameters
        ))
        dynamics.update(deltas)

        epsilons = dict(zip(
            "e_0, w_e, p_e".split(','),
            [0.0020, 0.86, np.pi]    # TODO: get actual parameters
        ))
        dynamics.update(epsilons)

    else:
        dynamics = dict(zip(
            sy_params_dynamic_names2()[:5],
            [0.0001, 100, 0.00035, 0.000035, 0.003]
        ))

    params.update(dynamics)
    return params