#' seirstrata: age- and region-stratified SEIR transmission model
#'
#' Build a model with [build_model()], run it with [run_model()], and look at
#' the results with [get_output()], `summary()` and `plot()`. The model
#' equations are described in
#' `docs/model_description.md`.
#'
#' @useDynLib seirstrata, .registration = TRUE
#' @import data.table
#' @importFrom ggplot2 .data
#' @keywords internal
"_PACKAGE"

# data.table non-standard evaluation
utils::globalVariables(c(".", "value", "variable", "time", "sim", "measure",
                         "median", "lower", "upper", ".SD", ".trace",
                         ".N", "draw", "predicted", "pred_lower", "pred_upper",
                         "admissions", "series", "step", "chain", "share",
                         "country", "region", "age_group", "population",
                         "profile", "area_km2", "density_per_km2",
                         ".value", "inc_infections", "inc_hosp",
                         "inc_deaths_comm", "date_onset", "date_hospital",
                         "date_death", "died", "id",
                         "hospitalised", "vaccinated", "age_mid", "scenario"))
