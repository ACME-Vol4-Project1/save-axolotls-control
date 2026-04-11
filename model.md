# Introduction

## Model 1.3.2 - Refined Version of Initial Model

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
- $T_1$: active frog sauna users
- $D$: deceased

Variable $T_2$ has been deprecated.

#### Control Variables
Resource Allocation (units of `$ / t`)
- $u_1$: frog sauna resource spend rate
- $u_2$: antifungal bath resource spend rate

### Equations

$$\begin{align*}
    \dot{S} &= - \alpha I S + \gamma T_1 + \eta \xi I - \zeta S \left(1 - \frac{T_1}{K}\right)\\
    \dot{I} &= \alpha I S - \zeta I \left(1 - \frac{T_1}{K}\right) - \eta \xi I - \beta I\\
    \dot{T_1} &= \zeta (I + S) \left(1 - \frac{T_1}{K}\right) - \gamma T_1\\
    \dot{D} &= \beta I\\
\end{align*}$$

Variable $T_2$ has been deprecated.

### Model Parameters

Controllable Parameters

name | units | description
--- | --- | --- 
$K$ | `f` | carrying capacity of frog saunas, a function of $u_1$. If $K$ is set to zero, we assume that no frogs are entering saunas through the carrying capacity term to avoid dividing by zero.
$\xi$ | `1 / t` | rate of transfer into antifungal treatment, a function of $u_2$

Constant Parameters

name | units | description
--- | --- | --- 
$\alpha$ | `1 / f / t` | infection rate from living hosts
$\beta$ | `1 / t` | death rate from disease
$\gamma$ | `1 / f / t` | rate that frogs leave sauna
$\zeta$ | `1 / t` | rate of transfer to heat treatment
$\eta$ | `1` (unitless) | success rate of antifungal bath treatment
$\theta$ | `1` (unitless) | radial measure of season

### Model Parameter Estimation

name | value | reasoning
--- | --- | --- 
$\theta$ | Set based on the time of year, and $\cos(\theta) \in [-1, 1]$ | $\cos(\theta) = 1$ corresponds to summer solstice, $\cos(\theta) = -1$ corresponds to winter solstice
$K$ | $\left(\frac{365}{2} u_1\right) \frac{\text{frogs} \cdot \text{days}}{\$}$ | Assuming we are averaging cost over a period of $~7$ years. According to Claude's summary of the sauna paper, 20 frogs were studied in the hot sauna in the sauna study. A Sauna costs $~\$150 \text{ USD}$ to build (rounding up from Claude estimate) and $20 \text{ USD}$ to maintain (rounding up again) every year with a lifetime with maintanance of $~7.5$ years, totaling a cost of $\$40$ USD per year or $\$40/365$ per day if we are considering a seven year period. The number of effective saunas we model should be proportional to the control expenditure rate divided by this cost rate. So the total carrying capacity should be $20 \times u_2 \times (365 / 40 / \$) = (365/2/\$) u_1$.
$\xi$ | $\frac{1}{(I + 50 \text{ frogs}) \cdot \$100} u_2$ | Assuming we can catch and treat a frog at the economic rate suggested by Claude. The scaling by $I + 50 \text{ frogs}$ holds the per-frog cost constant (in large number of infected frogs regime, when $I >> 50$) but decreases effectiveness when there are very few sick frogs (is more realistic - one sick frog is hard to find - and helps avoid $-\infty$ explosion of $I$). The cost of treatment was about $\$100$ USD per frog. For less economics-driven numbers, see Claude `literature_params_1-0-1.pdf`, these ($0.005 - 0.020 - 0.040 / \text{day}$) numbers were estimated by Claude from field studies. lowest corresponds to monthly surveys, the the highest to active trapping, in SPRING AND SUMMER. In winter, apparently, the frogs are much harder to catch, so we might want to include that in the model at some point.
$\alpha$ | $10^{-5} (1 + 0.5 \cos(\theta + 1.25 \pi)) / \text{frog} / \text{day}$ | see Claude `literature_params` parameter estimation documents - as for the seasonal component, see the season document. Transmission is supposed to peak 1-2 months before mortality. As for the mean value, I chose something near the middle of the general range it gave. This according to Claude doesn't have a direct estimation in the literature for the Green and Golden Bell Frogs, but the frogs get sick from much less fungus than most other frogs, so this infection parameter is on the higher range for all frog species.
$\beta$ | $\frac{(0.025 + .002)}{2} + \frac{(0.025 - .002)(\cos(\theta + \pi))}{2}/ \text{day}$ | $0.025 - .002$ are Winter - Summer estimates from Campbell et al. 2019 as extracted by Claude. We can interpolate with a cosine peaking in winter. Apparently the temperature (winter v summer) has an extremely high impact on the frogs' mortality rate from the fungus. Also, it seems that the Green and Golden Bell frogs have a very low disease-causing threshold of fungus, so they correspondingly have a higher mortality than many other species.  See `literature_params_1-0-1.pdf`.
$\gamma$ | $0.1 / \text{day}$ | Apparently it takes frogs 6-11 days in sauna for treatment. Let's pretend they hang out for 10 days before leaving. This could be tuned better.
$\zeta$ | $0.15 / \text{day}$ | base rate of entry to frog saunas - scaled down by carrying capacity. This is a spitball estimate picking the middle of a range from Claude.
$\eta$ | $1$ | Assuming that frog bath actually cures all the frogs without presenting any risks. This parameter isn't necessarily realistic and could be changed.

### Change Log
- Possible changes (prioritized based on some literature review)
    1. Add frog vaccination + vaccinated compartment. Possibly add recovered frogs + immunity from saunas as well. Probably also need baby frogs adn a spontaneous infection rate
    2. Maybe consider frog age stages + differential vulnerability
    3. Add reproduction and maybe natural death to healthy frogs
- Model 1.3.2: Add option `const_u1` for $u_1$ to be a constant tunable parameter rather than a control variable. May simplify the system for control problems. Should still impact the cost functional as if it were a function, but will require other methods to optimize.
- Model 1.3.1: Migrate to a sympy-inheriting ODE definition. Model reworking may be necessary if $u_1 = 0, T_1 \neq 0$ case becomes relevanat in optimal control solutions.
- Model 1.3: Added a seasonal component (merged 28 March).
- Model 1.2: Removed the $T_2$ variable to make bath treatment instantaneous. Also changed the parameterization of $\xi$ to give it a saturated and unsaturated regime (easier to catch sick frogs when there are a lot of them, to a point. Harder to catch frogs when there aren't). Also fixed erroneous estimate for $\alpha$ - it was orders of magnitude off by mistake.
- Model 1.1.1: Updated model carrying capacity terms to saturate in $T_1$ instead of $S + I$. Need to check parameters to make sure that is still valid. Also fixed incorrect control terms and fixed a dividing by zero issue when carrying capacity is zero.
- Model 1.1: We decided on in lab on 26 March, also pivoted to Golden and Green Bell frogs specifically.
- Model 1.0: Emeline's initial push

## Background Information

### Frog Treatments background

Hotspot shelters:
* __https://www.nature.com/articles/s41586-024-07582-y__  Study showed that sunlight-heated artificial refugia attract endangered frogs and enable body temperatures high enough to clear infections, should have details on type of shelter and materials

Antifungal bath:
* __https://conbio.onlinelibrary.wiley.com/doi/10.1111/csp2.12762__ study showed treated frogs were 4x more likely to survive first year, reduced Bd infection intensity (It seems like the age of the frog affects their susceptibility; it’s worse in tadpoles. Maybe we should incorporate that?)
__https://www.int-res.com/journals/dao/articles/dao02813__  treated frogs w/ antifungal drug itraconazole, treated frogs had a higher survival rate after 5 wks. Treatment reduced growth rate of frogs though (weird)


## Model 1.4.1: Reintroduction Model

Parameters needed:
- Natural birth and death rates (birth will be seasonal, death will not be)
- Distance between summer solstice and peak time
- Difference in infection rate for juvenile vs adult frogs
- Vulnerability measure for young frogs? 
- What age frogs can use the Saunas? tadpoles no but little frogs yes? What percentage? 
- What is the spontaneous infection rate (environment/other species) ? 
- Reasonable bounds on u_2, the rate of adding vaccinated frogs
- $\epsilon , \delta_b, \delta_d $ should depend on t for seasonality, 



### Equations

$$\begin{align*}
    \dot{S_1} &= \delta_1 (\delta_2 + V + I + T) - \varphi S_1 - \alpha_1 S_1 I - \epsilon S_1 - a \delta_2 S_1 - b \zeta \delta_1 (1 - \frac{T}{K - \zeta V}) \\
    \dot{S_2} &= \varphi S_1 - \alpha_2 S_2 I - \epsilon S_2 - \delta_2 S_2 - \zeta S_2 (1 - \frac{T}{K - \zeta V}) + \gamma (\frac{S_2 I}{V + S_2 + I}) T \\
    \dot{I} &= \alpha_2 S_2 I - \beta I + \alpha_1 S_1 I + \epsilon (S_1 + S_2)\\
    \dot{V} &= \nu V - \delta_2 V - \zeta V \left(1 - \frac{T}{K - \zeta V}\right) + \gamma (\frac{V}{V + S_2 + I})T\\
    \dot{T} &= \zeta (I + b S_1 + S_2 + V) \left(1 - \frac{T}{K - \zeta V}\right) - \gamma T\\
    \dot{D} &= a \delta_2 S_1 + \delta_2 S_2 + \beta I\\
\end{align*}$$


### Model Parameters

Controllable Parameters (units are not correct yet)

name | units | description
--- | --- | --- 
$K$ | `f` | carrying capacity of frog saunas, a function of $u_1$. If $K$ is set to zero, we assume that no frogs are entering saunas through the carrying capacity term to avoid dividing by zero.


Constant Parameters

name | units | description
--- | --- | --- 
$\alpha_1$ | `1 / f / t` | infection rate for baby frogs
$\alpha_2$ | `1 / f/ t` | infection rate for adult frogs
$\delta1$ | unknown | birth rate
$\delta2$ | unknown | death rate
$\varphi$ | unknown | rate the frogs grow up (mature)
$a$ | unknown | scale of frog death for juveniles
$b$ | unknown | fraction of juveniles that can go to sauna T
$\beta$ | `1 / t` | death rate from disease
$\gamma$ | `1 / f / t` | rate that frogs leave sauna
$\zeta$ | `1 / t` | rate of transfer to sauna T
$\theta$ | `1` (unitless) | radial measure of season