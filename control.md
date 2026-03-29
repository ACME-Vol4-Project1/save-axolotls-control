## Constraints 1.0 - Possible Cost Functionals and Constraints

### Functionals
Quadratic (satisfies LQR control assumptions)
$$J_{q}[x, u_1, u_2] = \int_0^\infty c D^2 + u_1^2 + u_2^2 dt$$

Linear (also may satisfy LQR control assumptions)
$$J_{l}[x, u_1, u_2] = \int_0^\infty c D + u_1 + u_2 dt$$

### Constraints

Possible Budget Constraint?
$$\int_0^\infty u_1 + u_2 dt \leq B$$

### Parameters

- $c$: represents relative cost of a frog life
- $B$: represents budget

### Change Log

- Possible changes
    - Can linear be written as LQR?
    - Can we implement budget constraints?
    - Add cost functional favoring more living frogs rather than fewer dead - could be especially important if we add frog growth term
- Model 1.0: Were written on board or discussed in lab on 26 March