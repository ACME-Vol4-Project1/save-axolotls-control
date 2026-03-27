# Introduction

## Model 1.1 - Initial Model

We are modeling Green + Golden Bell Frogs in Australia and Chitrid fungus infections. We are trying to use optimal control to prioritize treatment strategies for optimal effect in preserving frog populations. Our controllable parameters are the construction of hotspot shelters ('frog saunas') and the active treatment of infected frogs using antifungal baths.

### Units

- `f`: quantity of frogs
- `$`: cost
- `t`: time

### Variables

#### State Variables
Frog Compartments (units of `f`)
- $S$: susceptible
- $I$: infected
- $T_1$: treatment 1 (frog sauna)
- $T_2$: treatment 2 (antifungal bath)
- $D$: deceased

#### Control Variables
Resource Allocation (units of `$ / t`)
- $u_1$: frog sauna resource spend rate
- $u_2$: antifungal bath resource spend rate

### Equations

$$\begin{align*}
    \dot{S} &= - \alpha I S + \gamma T_1 + \eta T_2 - \zeta S \left(1 - \frac{S + I}{K}\right)\\
    \dot{I} &= \alpha I S - \zeta I \left(1 - \frac{S + I}{K}\right) - \xi I - \beta I\\
    \dot{T_1} &= \zeta (I + S) \left(1 - \frac{S + I}{K}\right) - \gamma T_1\\
    \dot{T_2} &= \xi I - \eta T_2 - (1 - \eta) T_2\\
    &= \xi I - T_2\\
    \dot{D} &= \beta I + (1 - \eta) T_2\\
\end{align*}$$

### Model Parameters

Controllable Parameters

name | units | description
--- | --- | --- 
$K$ | `f` | carrying capacity of frog saunas, a function of $u_1$
$\xi$ | `1 / t` | rate of transfer into antifungal treatment, a function of $u_2$

Constant Parameters

name | units | description
--- | --- | --- 
$\alpha$ | `1 / f / t` | infection rate from living hosts
$\beta$ | `1 / t` | death rate from disease
$\gamma$ | `1 / f / t` | rate that frogs leave sauna
$\zeta$ | `1 / t` | rate of transfer to heat treatment
$\eta$ | `1` (unitless) | success rate of antifungal bath treatment

### Model Parameter Estimation

name | value | reasoning
--- | --- | --- 
$K$ | $\left(\frac{365}{2} u_2\right) \frac{\text{frogs} \cdot \text{days}}{\$}$ | Assuming we are averaging cost over a period of $~7$ years. According to Claude's summary of the sauna paper, 20 frogs were studied in the hot sauna in the sauna study. A Sauna costs $~\$150 \text{ USD}$ to build (rounding up from Claude estimate) and $20 \text{ USD}$ to maintain (rounding up again) every year with a lifetime with maintanance of $~7.5$ years, totaling a cost of $\$40$ USD per year or $\$40/365$ per day if we are considering a seven year period. The number of effective saunas we model should be proportional to the control expenditure rate divided by this cost rate. So the total carrying capacity should be $20 \times u_2 \times (365 / 40 / \$) = (365/2/\$) u_2$.
$\xi$ | $\frac{1}{I \cdot \$100} u_2$ | Assuming we can catch and treat a frog at the economic rate suggested by Claude. The scaling by $I$ holds the per-frog cost constant no matter the number of frogs. The cost of treatment was about $\$100$ USD per frog. For less economics-driven numbers, see Claude `literature_params_1-0-1.pdf`, these ($0.005 - 0.020 - 0.040 / \text{day}$) numbers were estimated by Claude from field studies. lowest corresponds to monthly surveys, the the highest to active trapping, in SPRING AND SUMMER. In winter, apparently, the frogs are much harder to catch, so we might want to include that in the model at some point.
$\alpha$ | $10^{-3} / \text{frog} / \text{day}$ | see Claude `literature_params` parameter estimation documents - I chose something near the middle of the general range it gave. This according to Claude doesn't have a direct estimation in the literature for the Green and Golden Bell Frogs, but the frogs get sick from much less fungus than most other frogs, so this infection parameter is on the higher range for all frog species.
$\beta$ | $0.025 - .002 / \text{day}$ | Winter - Summer estimates from Campbell et al. 2019 as extracted by Claude. Apparently the temperature (winter v summer) has an extremely high impact on the frogs' mortality rate from the fungus. Also, it seems that the Green and Golden Bell frogs have a very low disease-causing threshold of fungus, so they correspondingly have a higher mortality than many other species.  See `literature_params_1-0-1.pdf`.
$\gamma$ | $0.1 / \text{day}$ | Apparently it takes frogs 6-11 days in sauna for treatment. Let's pretend they hang out for 10 days before leaving. This could be tuned better.
$\zeta$ | $0.15 / \text{day}$ | base rate of entry to frog saunas - scaled down by carrying capacity. This is a spitball estimate picking the middle of a range from Claude.
$\eta$ | $1$ | Assuming that frog bath actually cures all the frogs without presenting any risks. This parameter isn't necessarily realistic and could be changed.

### Change Log

- Possible changes (prioritized by Henry based on some literature review)
    1. Consider adding a seasonal component
    2. Consider removing $T_2$ in favor of 'instant' treatment
    3. Add reproduction and natural death to healthy frogs
    4. Possibly add recovered frogs + immunity from saunas
    5. Possibly add a variable for number of frog saunas
- Model 1.1: We decided on in lab on 26 March, also pivoted to Golden and Green Bell frogs specifically.
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