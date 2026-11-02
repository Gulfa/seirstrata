# seirstrata

An age- and region-stratified SEIR transmission model with hospitalisation and
deaths, written in [odin2](https://mrc-ide.github.io/odin2/) and run with
[dust2](https://mrc-ide.github.io/dust2/). It has a deterministic (ODE)
version and a stochastic (tau-leaping) version of the same model.

Derived from [seirvax](https://github.com/Gulfa/modelling_prep_ex) (commit
`64face5`). seirvax stays the preparedness-exercise package; `seirstrata` is
the production base model for the seasonal respiratory work in `sesong`.

## Scope

`seirstrata` contains the model and nothing specific to a pathogen, country
or fitting strategy:

- the transmission model, stratified by age, region and immunity stratum (a
  ladder of protection levels in S, waning down the ladder, and one fully
  protected R per group that wanes back into it);
- a global daily transmission multiplier (planned), so seasonality and
  holidays can be built outside the package;
- the observation model (the likelihood used for fitting), and an exported
  likelihood constructor (planned);
- running, output, summaries and plots.

Priors, samplers, choice of free parameters, seasonality and holiday
functions, and country data live in the packages that use it.

## Status

Work in progress. Steps 1.1 (copy from seirvax and trim) and 1.2 (immunity
strata, waning, single R) of the plan are done.
Fitting, the bundled fictional-country data and the line-list simulator from
seirvax have been removed.

## Development

```r
odin2::odin_package(".")   # regenerate C++ from inst/odin/*.R
devtools::document()       # regenerate NAMESPACE and man/
devtools::test()
```

At FHI this runs in the `fhi` docker container, where this directory is
mounted at `/fhi/seirstrata`. See
[`docs/modifying_the_model.md`](docs/modifying_the_model.md) and
[`docs/model_description.md`](docs/model_description.md).
