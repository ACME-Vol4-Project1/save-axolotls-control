# Introduction

## Model 1.1 - Initial Model

We are modeling Green + Golden Bell Frogs in Australia and Chitrid fungus infections. We are trying to use optimal control to prioritize treatment strategies for optimal effect in preserving frog populations. Our controllable parameters are the construction of hotspot shelters ('frog saunas') and the active treatment of infected frogs using antifungal baths.

### Variables

#### State Variables
Frog Compartments (units of `frogs`)
- $S$: susceptible
- $I$: infected
- $T_1$: treatment 1 (frog sauna)
- $T_2$: treatment 2 (antifungal bath)
- $R$: deceased

#### Control Variables
Resource Allocation (units of `$ / time`)
- $u_1$: frog sauna resource spend rate
- $u_2$: antifungal bath resource spend rate

### Equations

$$\begin{align*}
    \dot{S} &= - \alpha I S + \gamma T_1 + \eta(u_2) T_2 - \zeta S \left(1 - \frac{S + I}{K}\right)\\
    \dot{I} &= \alpha I S - \zeta I \left(1 - \frac{S + I}{K(u_1)}\right) - \xi I - \beta I\\
    \dot{T_1} &= \zeta (I + S) \left(1 - \frac{S + I}{K(u_1)}\right) - \gamma T_1\\
    \dot{T_2} &= \xi I - \eta(u_2) T_2 - (1 - \eta(u_2)) T_2\\
    &= \xi(u_2) I - T_2\\
    \dot{R} &= \beta I + (1 - \eta(u_2)) T_2\\
\end{align*}$$

### Model Parameters

Controllable Parameters
- $K$: carrying capacity of frog saunas, a function of $T_1$
- $\xi$: rate of transfer into antifungal treatment

Constant Parameters
- $\alpha$: infection rate from living hosts
- $\beta$: death rate from disease
- $\gamma$: infection rate from deceased hosts
- $\delta$: death rate of fungus in dead hosts
- $\zeta$: rate of transfer to heat treatment
- $\eta$: rate at which individuals leave treatment



### Change Log

- Possible changes
    - Add reproduction and natural death to healthy frogs
- Model 1.2 (proposed):
    - Possibly remove $T_2$
    - Possibly add a variable for number of frog saunas
- Model 1.1: We decided on in lab on 26 March
- Model 1.0: Emeline's initial push


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


## Background Information

### Frog Treatments background

Hotspot shelters:
* __https://www.nature.com/articles/s41586-024-07582-y__  Study showed that sunlight-heated artificial refugia attract endangered frogs and enable body temperatures high enough to clear infections, should have details on type of shelter and materials

Antifungal bath:
* __https://conbio.onlinelibrary.wiley.com/doi/10.1111/csp2.12762__ study showed treated frogs were 4x more likely to survive first year, reduced Bd infection intensity (It seems like the age of the frog affects their susceptibility; it’s worse in tadpoles. Maybe we should incorporate that?)
__https://www.int-res.com/journals/dao/articles/dao02813__  treated frogs w/ antifungal drug itraconazole, treated frogs had a higher survival rate after 5 wks. Treatment reduced growth rate of frogs though (weird)