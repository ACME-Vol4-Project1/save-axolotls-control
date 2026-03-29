# Constraints 1.1 - Cost Functionals and Constraints

## Functionals
1. Quadratic, infinite horizon

Satisfies LQR control infinite horizon assumptions if the ODE is autonomous, so we can use iterated LQR to optimize this cost functional with the autonomous model - but unsure if this works for the seasonal model

$$J[x, u_1, u_2] = \int_0^\infty c D^2 + u_1^2 + u_2^2 dt$$

2. Quadratic, finite horizon

Satisfies LQR control finite horizon assumptions (we think). We may be able to use finite horizon LQR to solve this one, possibly, even with the non-autonomous model.
$$J[x, u_1, u_2] = \int_0^{t_f} c D^2 + u_1^2 + u_2^2 dt$$

3. Linear (doesn't appear to satisfy LQR control assumptions, but we could try to solve with `solve_bvp`)
$$J_{l}[x, u_1, u_2] = \int_0^\infty c D + u_1 + u_2 dt$$

### Future options but not yet:
- Consider endpoint costs - we don't want the frogs to live for most of the interval and then die at the end...?
- Possible Budget Constraint? Like the following
$$\int_0^\infty u_1 + u_2 dt \leq B$$
- Should we phrase functionals in terms of more living rather than fewer dead frogs? Not as easy for LQR as we have things currently, but we could look at ways around it, or do it with `solve_bvp`

## Approaches 

### Solve via linearization + iterated LQR
- Jeremy and Josie are exploring linearization options + will discuss with Dr. Jarvis Monday
- Henry developed sympy version of $f$ so we can linearize it at each step with sympy. For iterated LQR we probably will need to stick to the autonomous model, so use the autonomous kwarg = True.

### Solve using `solve_bvp`
- Chammelia, Emeline, Henry going to explore this for at least one of the above functionals.

### Change Log
- Constraints 1.1: Cleaned up and focused current efforts on the quadratic LQR-appropriate functionals after team discussion on 28 March.
- Constraints 1.0: Were written on board or discussed in lab on 26 March