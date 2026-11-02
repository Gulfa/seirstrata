# seirstrata

An age- and region-stratified SEIR transmission model with hospitalisation and
deaths, written in [odin2](https://mrc-ide.github.io/odin2/) and run with
[dust2](https://mrc-ide.github.io/dust2/). It has a deterministic (ODE)
version and a stochastic (tau-leaping) version of the same model.

It started as a copy of the model in
[modelling_prep_ex](https://github.com/Gulfa/modelling_prep_ex), the
preparedness-exercise package, and is developed independently from it.
`seirstrata` is the base model for the seasonal respiratory work in `sesong`.

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

Work in progress: the immunity strata, waning and R are in place; the daily
transmission multiplier and the exported likelihood constructor are next.

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
