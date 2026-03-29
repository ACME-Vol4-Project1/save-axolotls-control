# Constraints 1.1 - Cost Functionals and Constraints

## 1. Fully Quadratic

$$J[x, u_1, u_2] = \int_0^{t_f} c D^2 + u_1^2 + u_2^2 dt$$

This is convex and quadratic, good candidate for LQR, squared money has pros and cons for interpretations as does squared dead frogs. But a good start.

### Planned Approaches
Satisfies LQR control infinite horizon assumptions if the ODE is autonomous, so we can use iterated LQR to optimize this cost functional with the autonomous model - but unsure if this works for the seasonal model. So possible approaches include

- **Iterated LQR**: autonomous dynamics
    - Change to infinite horizon assumption
    - Linearize every small amount of time and solve for small intervals

- **`solve_bvp`**: autonomous or seasonal dynamics
    - Need to find the Hamiltonian first, then specify endpoint + other constraints
    - Perhaps iterated LQR can provide us with a good initial guess.

### Potential approaches
- **Finite Time Horizon LQR**: autonomous (and seasonal?) dynamics
    - *Problem: Nonlinear dynamics. Can we solve by linearizing across small time intervals and maybe working backward from the endpoint?*
    - This would make endpoint cost easy to model
    - Requires formulating differential ricatti equation and find endpoint conditions
    - Need to solve differential ricatti equation - possible with `solve_ivp`?

- **Iterated LQR**: seasonal dynamics
    - *Problem: Dynamics change over time. Could we freeze time values across each small interval to make each interval autonomous?*
    - Change to infinite horizon assumption
    - Linearize every small amount of time and solve for small intervals

## 2. Partially or non-quadratic

$$J_{l}[x, u_1, u_2] = \int_0^\infty c D + u_1 + u_2 dt$$

### Planned Approaches
Satisfies LQR control infinite horizon assumptions if the ODE is autonomous, so we can use iterated LQR to optimize this cost functional with the autonomous model - but unsure if this works for the seasonal model. So possible approaches include

- **`solve_bvp`**: autonomous or seasonal dynamics
    -  Need to find the Hamiltonian first, then specify endpoint + other constraints

## Future options but not yet:
- Consider endpoint costs - we don't want the frogs to live for most of the interval and then die at the end...?
- Possible Budget Constraint? Like the following
$$\int_0^\infty u_1 + u_2 dt \leq B$$
- Should we phrase functionals in terms of more living rather than fewer dead frogs? Not as easy for LQR as we have things currently, but we could look at ways around it, or do it with `solve_bvp`

## Change Log
- Constraints 1.1: Cleaned up and focused current efforts on the quadratic LQR-appropriate functionals after team discussion on 28 March.
- Constraints 1.0: Were written on board or discussed in lab on 26 MarchSatisfies LQR control finite horizon assumptions (we think). We may be able to use finite horizon LQR to solve this one, possibly, even with the non-autonomous model, or we could use solve_bvp.