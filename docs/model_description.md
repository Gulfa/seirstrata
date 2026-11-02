# Model description

`seirstrata` contains one model in two versions: a deterministic one (ODEs) and a stochastic one (discrete-time binomial draws). Both versions share the same parameters and outputs.

## Structure

The population is divided into **groups** $g = 1,\dots,n$ and **immunity strata** $k = 1,\dots,K$ (`strata`, from no immunity upwards). In the R interface, groups are the combinations of age group $a$ and region $r$. The model code only sees $n = n_\text{age} \cdot n_\text{region}$ groups, numbered with age varying fastest: $g = (r-1)\,n_\text{age} + a$.

S, E, I and H are kept per group and stratum. R is one compartment per group:

```
            λ_g·sus_gk             σ                γ·p_hosp_gk                (1 - p_death_hosp_gk) / dur_hosp_g
  S_gk ───────────────► E_gk ──────────► I_gk ───────────────────► H_gk ──────────────────────────────────► R_g
   │  ▲                                    │                          │                                      │
   │  │ ω_g,k+1·S_g,k+1                    │                          └─ p_death_hosp_gk / dur_hosp_g ──► D_hosp_gk
   │  │ (waning from above)                ├─ γ·(1-p_hosp_gk)·p_death_comm_gk ──► D_comm_gk
   │  │                                    └─ γ·(1-p_hosp_gk)·(1-p_death_comm_gk) ──► R_g
   └──┴─ ω_gk·S_gk (waning to k-1)                                                                           │
                                                                       ω^R_g·R_g ──► S_g,R_to_stratum ◄──────┘
```

- **S:** susceptible, on a ladder of immunity strata. Stratum 1 has no immunity. Each stratum has its own susceptibility, infectiousness and severity.
- **E:** exposed (latent, not infectious). Mean duration is the latent period, $1/\sigma$.
- **I:** infectious. Mean duration is the infectious period, $1/\gamma$. When people leave I, they are hospitalised with probability `p_hosp`. Otherwise they die in the community with probability `p_death_comm` or recover.
- **H:** in hospital, not infectious. Mean length of stay is `dur_hosp`. On leaving, patients die with probability `p_death_hosp` or recover.
- **R:** recovered and fully protected. Everyone who recovers goes to R, whatever their stratum. R wanes at rate `waning_rate_R` into the ladder at `R_to_stratum` (default: the top stratum).
- **D_hosp, D_comm:** cumulative deaths in hospital and in the community.

An infection keeps the stratum of the person infected: E, I and H in stratum $k$ use that stratum's infectiousness and severity, so breakthrough infections are milder if the stratum says so.

**Waning down the ladder:** people in $S_{gk}$ move to $S_{g,k-1}$ at rate `waning_rate[k]` ($\omega_{gk}$). With the same rate on every step, the time to reach stratum 1 from stratum $k$ is gamma (Erlang) distributed, as with a chain of waning compartments.

**Initial state:** a fraction `initial_R` of each group starts in R. The rest are spread over the strata by `initial_strata`. `I0` and `E0` are placed in stratum 1, taken from its share.

Vaccination status is not tracked. Doses from `vacc_schedule` move people up the ladder (see [Vaccination over time](#vaccination-over-time)), and someone whose protection has waned below the vaccine stratum can be vaccinated again.

## Force of infection

$$
\lambda_g(t) = \beta\, b(t) \sum_{h} C_{gh}\, \min\big(m_g(t), m_h(t)\big)\; \frac{\sum_k \text{trans}_{hk} I_{hk}}{N_h}
$$

- $N_h$ is the living population of group $h$ ($S+E+I+H$ over all strata, plus $R$).
- $C_{gh}$ is the number of contacts per day a person in group $g$ has with people in group $h$.
- $\text{sus}_{gk}$ and $\text{trans}_{hk}$ are the relative susceptibility and infectiousness: the group's `susceptibility` / `transmissibility` times the stratum's `rel_susceptibility` / `rel_transmissibility`.
- $b(t)$ is the global transmission multiplier `beta_t` (see [Transmission multiplier over time](#transmission-multiplier-over-time)).
- $m_g(t)$ is the contact multiplier from `beta_schedule`: the fraction of group $g$'s contacts that still happen. Contacts between $g$ and $h$ are scaled by $\min(m_g, m_h)$, so the more restricted group decides:
  - A multiplier applied to everyone scales all transmission by exactly that factor.
  - "65+ at 0.5" halves every contact involving someone aged 65+, including their contacts with other age groups.
  - "Region A at 0.25" cuts all contacts involving region A residents to 25%, including their contacts with other regions.
  - Contact rates stay reciprocal, because the factor is the same from both sides.
  - Where two restricted groups meet, the stricter multiplier applies; the two are not multiplied together.

New infections in group $g$, stratum $k$, occur at rate $\lambda_g\,\text{sus}_{gk}\,S_{gk}$.

### Transmission multiplier over time

`beta_t` multiplies transmission for everyone. It is meant for seasonality, holidays and similar effects, which are built outside the package (in `sesong`) and passed in as numbers:

```r
days <- 0:209
seasonal <- 1 + 0.25 * cos(2 * pi * (days - 100) / 365)
holiday <- ifelse(days %in% 90:105, 0.85, 1)
build_model(..., beta_t = seasonal * holiday)          # one knot per day
build_model(..., beta_t = data.frame(time = c(0, 89, 90, 105, 106),
                                     mult = c(1, 1, 0.85, 0.85, 1)))
```

- **Knots, linear in between.** The values are knots: a daily vector has one at each day from 0, a data.frame one at each `time` (starting at 0). Between knots the multiplier is linear, so a sharp step is two knots one day apart. It is linear rather than piecewise constant because a step every day forces the ODE solver to restart at each day boundary: with a daily multiplier this is about 3 times faster.
- **After the last knot the last value holds.** For a daily vector that is usually a mistake, such as a projection longer than the vector built for the fit, so `run_model()` warns. A data.frame is assumed to hold its last value on purpose.
- **It combines with the targeted contact multipliers** from `beta_schedule`: transmission between $g$ and $h$ is scaled by $b(t) \min(m_g, m_h)$. Scenarios that set `beta_schedule` therefore keep the fitted seasonality.
- **`R0` does not include `beta_t`.** It is the reproduction number at $b = 1$, so R0 values are comparable between seasons whatever the multiplier is on the start day. `R_eff(model, time)` includes $b(\text{time})$.

### Combining age and region mixing

$$
C_{(a,r),(b,s)} = C^\text{age}_{ab} \cdot M^\text{region}_{rs}
\qquad\Longleftrightarrow\qquad
C = M^\text{region} \otimes C^\text{age}
$$

- `C_age[a, b]` is contacts per day an age-$a$ person has with age-$b$ people.
- `M_region[r, s]` is the fraction of a region-$r$ resident's contacts made with residents of region $s$. Rows sum to 1.

This assumes that everyone has the same age-mixing pattern and that the region matrix only redistributes where contacts happen. A 1×1 matrix removes the corresponding dimension.

### Mixing that changes over time

`mixing_schedule` switches the matrices at given times, so the *pattern* of contacts can change, not just its level:

```r
build_model(..., C_age = C_term, M_region = M_normal,
  mixing_schedule = list(
    list(time = 40, C_age = C_holiday),                  # schools close
    list(time = 40, M_region = M_restricted),            # ... and travel is restricted
    list(time = 70, C_age = C_term, M_region = M_normal) # both back to normal
  ))
```

Each entry replaces `C_age`, `M_region` or both from its time onwards; whichever is not given carries over. `C_age` and `M_region` themselves apply from time 0. Internally this is one `n × n × n_time` array, interpolated piecewise-constant, so the cost is a little memory and nothing at run time.

`beta_schedule` and `mixing_schedule` are independent and compose: the multipliers scale whichever matrix is in force. Use `mixing_schedule` when the *shape* changes (school holidays, closed borders) and `beta_schedule` when the *level* does (distancing, a group reducing its contacts).

**R0 is calibrated from the mixing at time 0.** `R_eff(model, time)` uses the matrix in force at `time`, so it reflects later changes.

## Vaccination over time

`vacc_schedule` gives doses per day from given times, for everyone or for chosen age groups and regions:

```r
build_model(...,
  strata = c("none", "partial", "vaccinated"),
  vacc_stratum = "vaccinated",
  vacc_schedule = data.frame(time  = c(40, 70),
                             doses = c(5000, 0),     # 5,000 a day from day 40 to day 70
                             age   = "65+"))
```

- Rows follow the same rules as `beta_schedule`: from `time` onwards, the groups a row matches (`NA`, or a missing `age`/`region` column, matches all) receive `doses` per day between them, shared in proportion to their population, until a later row matching them overrides it. Before the first row, no doses are given. `doses = 0` ends a campaign.
- Doses go to people in **S strata below `vacc_stratum`** and to people in **R**, in proportion to their size. Those in S move to `vacc_stratum`. Doses to R are counted but change nothing, as giving a dose to someone recently infected adds little. People in E, I or H, and in S at or above `vacc_stratum`, are not vaccinated.
- Vaccination status is not tracked, so people in R stay eligible. If the S strata below `vacc_stratum` run out, the remaining doses keep going to R and `cum_doses` is not capped by the population. With real dose data this is harmless; with made-up schedules, keep doses below the eligible population.
- **Protection is immediate**: a vaccinated person has the vaccine stratum's parameters from the moment of the dose. To approximate a delay until protection, start the schedule correspondingly later.
- In group $g$, the doses $d_g(t)$ are given at per-capita rate $\min\big(d_g(t) / (\sum_{k < k_v} S_{gk} + R_g),\ 10\big)$ per day among the eligible. This delivers exactly $d_g$ a day while eligible people remain.
- Doses given are recorded as `cum_doses` and `inc_doses`. They have no stratum (`"all"` in `get_output(by = "stratum")`). Initial immunity from `initial_strata` is not counted.
- `R_eff(model, time)` uses the initial strata only, not the doses given since.

## Stratum effects

Each stratum's parameters are the group's parameters times the stratum multipliers. Multipliers can be the same for all ages (a vector by stratum) or differ by age (an age × stratum matrix).

| Parameter in stratum $k$ | = group parameter × |
|---|---|
| `sus` | `rel_susceptibility[k]` |
| `trans` | `rel_transmissibility[k]` |
| `p_hosp` | `rel_hosp[k]` |
| `p_death_hosp`, `p_death_comm` | `rel_death[k]` |

`rel_hosp` and `rel_death` act **given infection**. The total protection against hospitalisation in stratum $k$ is $1 - \text{rel\_susceptibility}_k \cdot \text{rel\_hosp}_k$. To use a vaccine effectiveness (VE), set the multiplier to $1 - \text{VE}$.

## Calibrating β from R0

When `R0` is given, β is chosen so that the spectral radius of the next-generation matrix is `R0`. That radius is computed for a population entirely in stratum 1 (no immunity), with all contact multipliers set to 1:

$$
K_{gh} = \beta\, \text{sus}_{g1}\, C_{gh}\, \text{trans}_{h1}\, \frac{N_g}{N_h} \cdot \frac{1}{\gamma}
$$

`R_eff(model, time)` uses the same construction over all strata (types $(g, k)$), with the initial susceptible populations in each stratum (R excluded) and the multipliers at `time`. `print(model)` shows both values.

**Known limitation:** the next-generation R0 describes generation-to-generation growth once infections follow the dominant eigenvector's age and region distribution. If the epidemic is seeded differently (for example in one region or one age group), early growth can be faster or slower than R0 suggests until that distribution is reached. The same applies to the calibrated β.

## Stochastic version

The stochastic model uses tau-leaping with a fixed time step `dt`. Each step:

| Transition | Draw |
|---|---|
| leaving S | $n_S = \text{Bin}(S, 1 - e^{-(\lambda\,\text{sus} + \omega + v)\,dt})$, with $\omega$ the waning rate and $v$ the vaccination rate (0 at or above the vaccine stratum) |
| S → E | $\text{Bin}(n_S, \lambda\,\text{sus} / (\lambda\,\text{sus} + \omega + v))$ |
| S → stratum below | $\text{Bin}(n_S - n_{SE}, \omega / (\omega + v))$; the rest is vaccinated |
| leaving R | $n_R = \text{Bin}(R, 1 - e^{-(\omega^R + v)\,dt})$ |
| R → ladder | $\text{Bin}(n_R, \omega^R / (\omega^R + v))$; the rest are doses given to R, who stay in R |
| E → I | $\text{Bin}(E, 1 - e^{-\sigma dt})$ |
| leaving I | $n_I = \text{Bin}(I, 1 - e^{-\gamma dt})$ |
| I → H | $\text{Bin}(n_I, p_\text{hosp})$ |
| I → D_comm | $\text{Bin}(n_I - n_{IH}, p_\text{death,comm})$ |
| leaving H | $n_H = \text{Bin}(H, 1 - e^{-dt/\text{dur\_hosp}})$ |
| H → D_hosp | $\text{Bin}(n_H, p_\text{death,hosp})$ |

The chained binomials guarantee that competing outcomes never exceed the number of people leaving a compartment, and counts stay whole numbers.

**Time-step bias.** In discrete time, a person stays in each stage for $dt/(1-e^{-\text{rate}\cdot dt}) \approx 1/\text{rate} + dt/2$ on average. The infectious period is therefore slightly too long, and epidemics are slightly larger than in the ODE model. For R0 = 2 the mean final size is:

| `dt` | final size vs. ODE |
|---|---|
| 1 | +6.1% |
| 0.5 | +3.2% |
| 0.25 | +1.7% |
| 0.1 (default) | +0.7% |
| 0.05 | +0.4% |

## Plot styling

Plots use ColorBrewer **Dark2** for groups (age, region, scenario, MCMC chain). With more than eight groups the palette is interpolated rather than recycled, so no two groups share a colour. Spaghetti lines fade as they get more numerous but never below 30% opacity; `alpha` overrides that on `plot()`. Maps keep viridis, which suits a continuous scale.

## Outputs

All outputs are recorded daily and stratified by age, region and immunity stratum (R, doses and weekly admissions have no stratum). [`get_output()`](../R/output.R) sums over any strata you don't ask for.

| Variable | Meaning |
|---|---|
| `S`, `E`, `I`, `H`, `R` | number in each compartment (`H` = hospital occupancy) |
| `D_hosp`, `D_comm`, `D` | cumulative deaths (hospital, community, both) |
| `cum_infections`, `cum_hosp`, `cum_doses` | cumulative infections, hospital admissions, doses given by `vacc_schedule` |
| `inc_infections`, `inc_hosp`, `inc_deaths_hosp`, `inc_deaths_comm`, `inc_deaths`, `inc_doses` | new infections / admissions / deaths / doses during day $(t-1, t]$ |
| `inc_hosp_week` | admissions during $(t-7, t]$; meaningful at $t = 7, 14, \dots$ |

Initial infections (`I0`, `E0`) are not counted in `cum_infections` or `inc_infections`.

## Fitting to hospital admissions

**The likelihood is part of the model.** The last block of both odin files compares modelled admissions with observed admissions through a negative binomial distribution, with size `nb_size` (variance $\mu + \mu^2/\text{nb\_size}$). dust2 evaluates it, as an "unfilter" for the ODE model and as a particle filter for the stochastic model.

**Daily or weekly.** The models count admissions both per day (`inc_hosp`, reset every day) and per week (`inc_hosp_week`, by group, reset at days 7, 14, …). The parameter `obs_weekly` selects which is compared with the data. Weekly data are admissions during $(t-7, t]$ at $t = 7, 14, \dots$, so **model day 0 must be the first day of a reporting week**, which is a Monday for ISO weeks.

**Several series.** Data can be split by age and/or region. Each series is the weighted sum of the model groups that `obs_map[series, group]` maps to it. With `model_data(..., age_groups = )`, several model age groups can count towards one coarser data age group.

**From fitted values to model parameters.** `model_likelihood(model, data, pars, process, initial)` returns a monty model. Each fitted parameter is either a scalar model parameter used directly (e.g. `beta`, `nb_size`), or goes through a `process` function that returns model parameters to replace, at the level of `model$pars`:

```r
lik <- model_likelihood(m, data, pars = c("scale", "amp", "severity", "nb_size"),
  process = function(x) list(
    beta = m$beta * x$scale,
    beta_t_value = my_seasonality(x$amp),       # built outside the package
    p_hosp = m$pars$p_hosp * x$severity),
  initial = c(scale = 1, amp = 0.2, severity = 1, nb_size = 20))
prior <- monty::monty_dsl({ scale ~ Uniform(0.5, 2); ... })
posterior <- lik$likelihood + prior
```

`process` must keep parameter sizes fixed (e.g. the number of `beta_t` knots), because array sizes are fixed when the model is compiled for a run. It works below `build_model()`, so nothing derived is recalculated: to fit R0 rather than `beta`, convert in `process`.

**Projections.** `likelihood_pars(lik, x)` gives the model parameters that the fit sets at a value `x`, and `run_model(model, pars = ...)` runs with them. `model` can be a scenario variant of the fitted model.

Priors, samplers and the choice of free parameters belong to the packages that use `seirstrata`.
