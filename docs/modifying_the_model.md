# Modifying the model

The model equations are in two [odin2](https://mrc-ide.github.io/odin2/) files:

| File | Version |
|---|---|
| `inst/odin/seir_ode.R` | deterministic: `deriv()` equations |
| `inst/odin/seir_stoch.R` | stochastic: `update()` equations with `Binomial()` draws |

odin2 translates these into C++, which dust2 compiles and runs. The generated code lives in `R/dust.R`, `R/cpp11.R`, `src/` and `inst/dust/`. **Never edit these by hand**; regenerate them instead (step 3 below).

The R layer takes care of the rest:

| File | Role |
|---|---|
| `R/build_model.R` | turns user arguments into the parameter list the odin models expect (`model$pars`) |
| `R/params.R` | expands parameters to age × region, builds the contact multiplier and vaccination schedules |
| `R/mixing.R` | combines the age and region mixing matrices |
| `R/r0.R` | next-generation matrix, R0 and R_eff |
| `R/run.R` | creates and runs the dust2 system, reshapes output to `[age, region, stratum, sim, time]` |
| `R/output.R`, `R/summary_plot.R` | `get_output()`, `summary()`, `plot()` |

## Development setup

Development needs `odin2` (to regenerate the C++ code) and `devtools`, on top of the packages in `DESCRIPTION`:

```r
install.packages(c("odin2", "dust2", "monty"),
                 repos = c("https://mrc-ide.r-universe.dev",
                           "https://cloud.r-project.org"))
install.packages("devtools")
```

The usual development commands, run from the package root:

```r
odin2::odin_package(".")   # regenerate C++ from inst/odin/*.R
devtools::document()       # regenerate NAMESPACE and man/
devtools::load_all()       # load the package for interactive work
devtools::test()           # run the tests
devtools::check()          # full R CMD check
devtools::install()        # install it
```

Note that `load_all()` compiles without optimisation, so the stochastic model and particle MCMC are several times slower than in an installed build. Install the package before timing anything.

## Checklist for a change

1. **Edit both odin files the same way.** They must keep the same parameters and the same states, in the same order. `run_model()` relies on this to treat the two versions identically. In `seir_stoch.R`:
   - turn every rate into a per-step probability, `1 - exp(-rate * dt)`;
   - draw counts with `Binomial()`;
   - split competing outcomes with chained binomials (see how leaving `I` is handled).
2. **Keep states as `[n_group, n_imm]` arrays** (group × immunity stratum) or `[n_group]` arrays (no stratum, like `R`). Then `run_model()`, `get_output()` and `plot()` pick up new states automatically: every state from `dust_unpack_state()` is reshaped to `[age, region, stratum, sim, time]`, with stratum `"all"` for `[n_group]` states.
3. **Regenerate the code:** `odin2::odin_package(".")`, then reinstall or use `devtools::load_all()`.
4. **Pass any new parameter from R:**
   - add an argument to `build_model()`;
   - expand it with `g(<arg>, "<name>", ...)` if it varies by group;
   - add it to the `pars` list;
   - document it with a roxygen `@param` line, then run `devtools::document()`.

   Any argument of `build_model()` automatically works with `update_model()`.
5. **Update the next-generation matrix** in `R/r0.R` if the change affects who can infect whom. Examples: a new infectious compartment, or a change to the force of infection. Otherwise R0 calibration silently becomes wrong.
6. **Update `docs/model_description.md`** and add a test in `tests/testthat/`. At minimum:
   - population is still conserved (`test-model.R`);
   - deterministic and stochastic still agree.

## Worked example: adding an asymptomatic compartment

Suppose a fraction `p_asymp` of infections are asymptomatic. They are infectious with relative infectiousness `rel_inf_asymp`, and are never hospitalised and never die.

In **`seir_ode.R`**:

```r
p_asymp       <- parameter(0)
rel_inf_asymp <- parameter(1)

# E splits into I and A
n_EA[, ] <- n_EI_total[i, j] * p_asymp
n_EI[, ] <- n_EI_total[i, j] - n_EA[i, j]
n_EI_total[, ] <- sigma * E[i, j]
n_A_out[, ] <- gamma * A[i, j]
dim(n_EA, n_EI_total, n_A_out) <- c(n_group, 2)

initial(A[, ]) <- 0
dim(A) <- c(n_group, 2)
deriv(A[, ]) <- n_EA[i, j] - n_A_out[i, j]
deriv(R[, ]) <- ... + n_A_out[i, j]          # add to the existing R equation

# A counts in N and contributes to transmission
N[] <- ... + sum(A[i, ])
prev[] <- (trans[i, 1] * (I[i, 1] + rel_inf_asymp * A[i, 1]) +
           trans[i, 2] * (I[i, 2] + rel_inf_asymp * A[i, 2])) / N[i]
```

In **`seir_stoch.R`**, make the same changes with binomial draws:

```r
n_EI_total[, ] <- Binomial(E[i, j], p_EI)
n_EA[, ] <- Binomial(n_EI_total[i, j], p_asymp)
n_EI[, ] <- n_EI_total[i, j] - n_EA[i, j]
n_A_out[, ] <- Binomial(A[i, j], p_I_out)
update(A[, ]) <- A[i, j] + n_EA[i, j] - n_A_out[i, j]
```

Then:
- add `p_asymp` and `rel_inf_asymp` to `build_model()` and its `pars` list;
- in `R/r0.R`, multiply `trans` by `(1 - p_asymp) + p_asymp * rel_inf_asymp` inside the next-generation matrix (asymptomatic and symptomatic periods are equally long here);
- regenerate, document and test.

## Changing the observation model

The comparison with data is the last block of both odin files: `admissions <- data()` and `admissions[] ~ NegativeBinomial(...)`, with `adm_group` choosing daily (`inc_hosp`) or weekly (`inc_hosp_week`) admissions per group. To fit to a different quantity, for example deaths:

1. Add a weekly accumulator for it (like `inc_hosp_week`, with `zero_every = 7`) and use it in `adm_group` in both files.
2. Regenerate the code.
3. Rename the data column where `R/likelihood.R` checks and builds `admissions`.

To fit two data streams at once, add a second `data()` and `~` line, and extend `model_data()` to build both.

**odin2 pitfall:** a literal whole number inside a `~` comparison, such as `NegativeBinomial(size = 10, ...)`, fails to compile ("deduced conflicting types ‘double’ and ‘int’"). Make the value a `parameter()` instead, as is done for `nb_size`.

## Useful odin2 syntax

| Syntax | Meaning |
|---|---|
| `x <- parameter()`, `parameter(0.5)` | required parameter / one with a default |
| `n <- parameter(type = "integer", constant = TRUE)` | a dimension |
| `dim(a, b) <- c(n_group, 2)` | declare array dimensions (several at once) |
| `x[, ] <- y[i, j] * z[i]` | array expressions; `i`, `j` index the dimensions |
| `sum(x[i, ])` | partial sums |
| `initial(x[, ], zero_every = 1) <- 0` | reset every day (daily incidence) |
| `interpolate(t, y, "constant")` | piecewise-constant time-varying input |
| `x <- data()`, `x ~ Poisson(mu)` | data and its likelihood (used by dust2 filters) |

See the [odin2 documentation](https://mrc-ide.github.io/odin2/) and the [dust2 documentation](https://mrc-ide.github.io/dust2/) for more.
