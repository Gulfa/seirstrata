output_dims <- c("age", "region", "stratum")

#' Extract results as a table
#'
#' Output variables (all per day, from time 0 to `t_max`):
#'
#' * Compartments (prevalence): `S`, `E`, `I`, `H` (hospital occupancy), `R`,
#'   and cumulative deaths `D_hosp`, `D_comm`, `D` (= both).
#' * Cumulative counts: `cum_infections`, `cum_hosp` (hospital admissions),
#'   `cum_doses` (doses given by `vacc_schedule`, including those given to
#'   people in `R`).
#'
#' Infections, admissions and deaths are recorded in the immunity stratum of
#' the person infected. `R` and the doses have no stratum: with
#' `by = "stratum"` they appear under stratum `"all"`.
#' * Daily incidence during day `(t - 1, t]`: `inc_infections`, `inc_hosp`
#'   (admissions), `inc_deaths_hosp`, `inc_deaths_comm`, `inc_deaths`,
#'   `inc_doses`.
#'
#' @param run Result of [run_model()].
#' @param vars Output variables to extract (default: all).
#' @param by Stratifications to keep, any of `"age"`, `"region"` and
#'   `"stratum"`. Everything else is summed over. The default,
#'   `character(0)`, gives population totals.
#'
#' @return A `data.table` with columns `time`, `sim`, the `by` columns
#'   (factors, in model order), `variable` and `value`. Deterministic runs
#'   have `sim = 1`. For a comparison there is also a `scenario` column.
#' @export
get_output <- function(run, vars = NULL, by = character(0)) {
  UseMethod("get_output")
}

#' @export
get_output.default <- function(run, vars = NULL, by = character(0)) {
  stop("'run' must come from run_model() or compare_scenarios()")
}

#' @rdname get_output
#' @export
get_output.seirstrata_run <- function(run, vars = NULL, by = character(0)) {
  if (is.null(vars)) {
    vars <- names(run$states)
  }
  unknown <- setdiff(vars, names(run$states))
  if (length(unknown) > 0) {
    stop(sprintf("Unknown variable(s): %s. Available: %s",
                 paste(unknown, collapse = ", "),
                 paste(names(run$states), collapse = ", ")))
  }
  by <- check_by(by)
  keep <- c(match(by, output_dims), 4L, 5L)
  drop <- setdiff(1:5, keep)

  res <- lapply(vars, function(v) {
    x <- aggregate_array(run$states[[v]], keep, drop)
    dt <- as.data.table(as.table(x), stringsAsFactors = FALSE)
    setnames(dt, "N", "value")
    dt[, variable := v]
    dt
  })
  res <- rbindlist(res)
  levels <- dimnames(run$states[[1]])
  levels$stratum <- c(run$model$dims$stratum, "all")
  for (b in by) {
    set(res, j = b, value = factor(res[[b]], levels[[b]]))
  }
  res[, time := as.numeric(time)]
  res[, sim := as.integer(sim)]
  setcolorder(res, c("time", "sim", by, "variable", "value"))
  setorderv(res, c("variable", by, "sim", "time"))
  res[]
}

check_by <- function(by) {
  by <- unique(as.character(by))
  bad <- setdiff(by, output_dims)
  if (length(bad) > 0) {
    stop(sprintf("'by' must be a subset of: %s",
                 paste(output_dims, collapse = ", ")))
  }
  output_dims[output_dims %in% by]
}

# Sum an array over the `drop` dimensions, keeping `keep` (in order).
aggregate_array <- function(x, keep, drop) {
  if (length(drop) == 0) {
    return(aperm(x, keep))
  }
  y <- colSums(aperm(x, c(drop, keep)), dims = length(drop))
  array(y, dim(x)[keep], dimnames(x)[keep])
}
