import numpy as np
from naming import sy_params_static, sy_params_dynamic_names 

### Default Parameters
def default_params(seasonal=True):
    # setup parameter dictionary for sympy substitution
    params = dict()

    # parameters that don't differ between seasonal and not
    statics = dict(zip(
        sy_params_static(), # γ,ζ,η
        [0.1, 0.15, 1.]
        ))
    params.update(statics)

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