# Comparing scenarios: several models, runs or projections side by side.

#' Compare scenarios
#'
#' Collects several scenarios into one object that `summary()`, `plot()` and
#' [get_output()] treat as a single thing, with a `scenario` column.
#'
#' Scenarios can be given as models (which are run for you) or as runs from
#' [run_model()], e.g. projections over posterior draws made by a fitting
#' package. They must share the same age groups and regions.
#'
#' @param scenarios A named list of models or runs. Names become
#'   the scenario labels, in the order given.
#' @param t_max,type,n_sim Passed to [run_model()] for any entry that is a
#'   model rather than an already-run scenario.
#' @param vars Output variables to keep when running models, as in
#'   [run_model()]. `summary()` needs `cum_infections`, `cum_hosp`, `D` and
#'   `H`, which the default keeps.
#' @param ... Further arguments for [run_model()], e.g. `seed`, `dt`,
#'   `n_threads`.
#'
#' @return An object of class `seirstrata_comparison`.
#' @export
#' @examples
#' m <- build_model(N = 1e5, R0 = 2, I0 = 10, p_hosp = 0.02)
#' cmp <- compare_scenarios(list(
#'   baseline = m,
#'   distancing = update_model(m, beta_schedule = data.frame(time = 40, mult = 0.7))),
#'   t_max = 250)
#' summary(cmp)
#' plot(cmp, var = "H")
compare_scenarios <- function(scenarios, t_max = 365,
                              type = c("deterministic", "stochastic"),
                              n_sim = 100, vars = NULL, ...) {
  type <- match.arg(type)
  if (!is.list(scenarios) || length(scenarios) < 2) {
    stop("'scenarios' must be a list of at least two scenarios")
  }
  labels <- names(scenarios)
  if (is.null(labels) || any(!nzchar(labels)) || anyDuplicated(labels)) {
    stop("'scenarios' must have unique names, which label the scenarios")
  }
  runs <- lapply(scenarios, function(s) {
    if (inherits(s, "seirstrata_model")) {
      run_model(s, type, t_max = t_max, n_sim = n_sim, vars = vars, ...)
    } else if (inherits(s, "seirstrata_run")) {
      s
    } else {
      stop("Each scenario must be a model or a run from run_model()")
    }
  })
  dims <- lapply(runs, function(r) r$model$dims)
  if (!all(vapply(dims, identical, TRUE, dims[[1]]))) {
    stop("All scenarios must have the same age groups and regions")
  }
  structure(list(runs = runs, scenarios = labels), class = "seirstrata_comparison")
}

#' @export
print.seirstrata_comparison <- function(x, ...) {
  cat(sprintf("<seirstrata_comparison> %d scenarios: %s\n", length(x$runs),
              paste(x$scenarios, collapse = ", ")))
  for (i in seq_along(x$runs)) {
    r <- x$runs[[i]]
    cat(sprintf("  %-20s %s, day %d to %d, %d simulation(s)\n", x$scenarios[i],
                r$type, min(r$time), max(r$time), r$n_sim))
  }
  cat("  Use get_output(), summary() and plot() to compare them\n")
  invisible(x)
}

#' @rdname get_output
#' @export
get_output.seirstrata_comparison <- function(run, vars = NULL,
                                          by = character(0)) {
  d <- rbindlist(lapply(run$scenarios, function(s) {
    get_output(run$runs[[s]], vars, by)[, scenario := s]
  }))
  d[, scenario := factor(scenario, run$scenarios)]
  setcolorder(d, c("scenario", "time"))
  d[]
}

#' Summarise a comparison
#'
#' The outcomes of [summary.seirstrata_run()] for every scenario, in one table.
#'
#' @param object A comparison from [compare_scenarios()].
#' @param by Stratifications, as in [get_output()].
#' @param level Width of the interval for scenarios with several simulations.
#' @param ... Unused.
#'
#' @return A `data.table` with a `scenario` column.
#' @export
summary.seirstrata_comparison <- function(object, by = character(0),
                                       level = 0.95, ...) {
  d <- rbindlist(lapply(object$scenarios, function(s) {
    summary(object$runs[[s]], by = by, level = level)[, scenario := s]
  }), fill = TRUE)
  d[, scenario := factor(scenario, object$scenarios)]
  setcolorder(d, c("scenario", by, "measure"))
  setorderv(d, c("measure", by, "scenario"))
  d[]
}

#' Plot a comparison
#'
#' One variable over time, coloured by scenario. Scenarios with several
#' simulations (stochastic runs, or projections over posterior draws) are
#' shown as the median with a central interval, since overlaid spaghetti from
#' several scenarios is hard to read; `style = "spaghetti"` draws the
#' individual trajectories instead.
#'
#' @param x A comparison from [compare_scenarios()].
#' @param var Output variable, see [get_output()].
#' @param by Stratification to facet by, as in [get_output()].
#' @param style `"ribbon"` or `"spaghetti"`.
#' @param level Width of the interval.
#' @param n_traces Maximum number of trajectories per scenario for
#'   `style = "spaghetti"`.
#' @param alpha Opacity of the individual trajectories (`NULL` chooses by how
#'   many there are).
#' @param ... Unused.
#'
#' @return A ggplot object.
#' @export
plot.seirstrata_comparison <- function(x, var = "H", by = character(0),
                                    style = c("ribbon", "spaghetti"),
                                    level = 0.95, n_traces = 20,
                                    alpha = NULL, ...) {
  if (length(var) != 1) {
    stop("'var' must be a single variable")
  }
  style <- match.arg(style)
  by <- check_by(by)
  d <- get_output(x, var, by)

  if (style == "spaghetti") {
    d[, .trace := paste(scenario, sim)]
    p <- ggplot2::ggplot(d) +
      ggplot2::geom_line(ggplot2::aes(time, value, colour = scenario,
                                      group = .data$.trace),
                         data = d[sim <= n_traces],
                         alpha = alpha %||% trace_alpha(n_traces),
                         linewidth = 0.4)
  } else {
    a <- (1 - level) / 2
    q <- d[, .(median = stats::median(value),
               lower = stats::quantile(value, a, names = FALSE),
               upper = stats::quantile(value, 1 - a, names = FALSE)),
           by = c("scenario", "time", by)]
    p <- ggplot2::ggplot(q) +
      ggplot2::geom_ribbon(ggplot2::aes(time, ymin = lower, ymax = upper,
                                        fill = scenario), alpha = 0.2) +
      ggplot2::geom_line(ggplot2::aes(time, median, colour = scenario))
  }
  if (length(by) > 0) {
    p <- p + ggplot2::facet_wrap(by, scales = "free_y")
  }
  p + dark2_scales(x$scenarios) +
    ggplot2::labs(x = "Day", y = var, colour = NULL, fill = NULL) +
    ggplot2::theme_bw()
}
