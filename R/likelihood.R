# Data, parameters and likelihood for fitting.
#
# The comparison with data is part of the odin models (the end of
# inst/odin/*.R): modelled admissions are compared with observed admissions
# through a negative binomial distribution. dust2 evaluates it, as an
# "unfilter" for the deterministic model and a particle filter for the
# stochastic one. This file prepares the data, builds the parameter mapping
# and returns a monty likelihood. Priors, samplers and the choice of free
# parameters belong to the packages that use seirstrata.

#' Prepare admissions data for the likelihood
#'
#' Checks observed hospital admissions and converts them to the format the
#' model's comparison uses: one row per time, with one value per data series.
#'
#' @param data A data.frame with columns `time` and `admissions`, and
#'   optionally `age` and/or `region` to give separate series (e.g. admissions
#'   by age group). Without them the data are national totals. `NA`
#'   admissions are missing and do not contribute to the likelihood.
#' @param model A model from [build_model()].
#' @param period `"week"` or `"day"`. Weekly data are admissions during days
#'   `(time - 7, time]` and must be at `time` = 7, 14, ..., so model day 0 is
#'   the first day of a reporting week (a Monday for ISO weeks). Daily data
#'   are admissions during `(time - 1, time]`, at whole days >= 1.
#' @param age_groups Optional mapping from the model's age groups to coarser
#'   age groups in the data: a character vector named with the model's age
#'   labels, whose values are data age labels, e.g.
#'   `c("0-4" = "0-14", "5-14" = "0-14", "15-64" = "15-64")`. Without it, the
#'   data's `age` labels must be the model's.
#'
#' @return An object of class `seirstrata_data`.
#' @export
model_data <- function(data, model, period = c("week", "day"),
                       age_groups = NULL) {
  if (!inherits(model, "seirstrata_model")) {
    stop("'model' must be created with build_model()")
  }
  period <- match.arg(period)
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  if (!all(c("time", "admissions") %in% names(data))) {
    stop("'data' must have columns 'time' and 'admissions'")
  }
  extra <- setdiff(names(data), c("time", "admissions", "age", "region"))
  if (length(extra) > 0) {
    stop(sprintf("Unknown column(s) in 'data': %s", paste(extra, collapse = ", ")))
  }
  step <- if (period == "week") 7 else 1
  if (anyNA(data$time) || any(data$time < step) ||
      any(data$time %% step != 0)) {
    stop(if (period == "week")
      "Weekly 'data$time' must be 7, 14, ... (day 0 starts a reporting week)"
      else "Daily 'data$time' must be whole days >= 1 (day 0 is the model start)")
  }
  a <- data$admissions[!is.na(data$admissions)]
  if (any(a < 0) || any(a != round(a))) {
    stop("'data$admissions' must be non-negative whole numbers (or NA)")
  }

  dims <- model$dims
  group_age <- rep(dims$age, length(dims$region))
  group_region <- rep(dims$region, each = length(dims$age))
  if (!is.null(age_groups)) {
    if (!is.character(age_groups) || is.null(names(age_groups)) ||
        !setequal(names(age_groups), dims$age)) {
      stop("'age_groups' must be a character vector named with every model age group")
    }
    group_age <- unname(age_groups[group_age])
  }
  labels <- list(age = unique(group_age), region = dims$region)

  keys <- intersect(c("age", "region"), names(data))
  for (k in keys) {
    data[[k]] <- as.character(data[[k]])
    bad <- setdiff(stats::na.omit(data[[k]]), labels[[k]])
    if (length(bad) > 0) {
      stop(sprintf("Unknown %s in 'data': %s", k, paste(bad, collapse = ", ")))
    }
  }
  if (length(keys) == 0) {
    series <- data.frame(label = "all")
    data$.series <- 1L
  } else {
    series <- unique(data[keys])
    rownames(series) <- NULL
    key <- do.call(paste, c(data[keys], sep = "\r"))
    data$.series <- match(key, do.call(paste, c(series, sep = "\r")))
    series$label <- apply(series[keys], 1, function(x) {
      paste(ifelse(is.na(x), "all", x), collapse = " / ")
    })
  }
  if (anyDuplicated(data[c("time", ".series")])) {
    stop("'data' has more than one row for the same time and series")
  }

  n_obs <- nrow(series)
  obs_map <- matrix(0, n_obs, length(group_age))
  for (k in seq_len(n_obs)) {
    match <- rep(TRUE, length(group_age))
    if (!is.null(series$age) && !is.na(series$age[k])) {
      match <- match & group_age == series$age[k]
    }
    if (!is.null(series$region) && !is.na(series$region[k])) {
      match <- match & group_region == series$region[k]
    }
    obs_map[k, match] <- 1
  }

  times <- sort(unique(data$time))
  y <- matrix(NA_real_, length(times), n_obs)
  y[cbind(match(data$time, times), data$.series)] <- data$admissions
  data_dust <- data.frame(time = times)
  data_dust$admissions <- I(lapply(seq_along(times), function(i) y[i, ]))

  data$series <- series$label[data$.series]
  data$.series <- NULL
  structure(list(data_dust = data_dust, obs_map = obs_map, series = series,
                 data_long = data, period = period, dims = dims),
            class = "seirstrata_data")
}

#' @export
print.seirstrata_data <- function(x, ...) {
  cat(sprintf("<seirstrata_data> %s admissions, %d series, %d time point(s) (days %g to %g)\n",
              if (x$period == "week") "weekly" else "daily", nrow(x$series),
              nrow(x$data_dust), min(x$data_dust$time), max(x$data_dust$time)))
  cat(sprintf("  Series: %s\n", paste(x$series$label, collapse = ", ")))
  invisible(x)
}

#' Model parameters with changes
#'
#' The full list of parameters the compiled model uses (`model$pars`), with
#' some of them replaced. This is the level at which the likelihood's
#' `process` function works, and what [run_model()] takes as `pars`: changes
#' here skip [build_model()], so they are fast, but nothing derived from them
#' is recalculated (e.g. `beta` is not recalibrated from R0).
#'
#' @param model A model from [build_model()].
#' @param overrides A named list of parameters to replace. Each must exist in
#'   `model$pars` and keep its size: array sizes are fixed when the model is
#'   compiled for a run.
#' @return A named list of parameters.
#' @export
model_pars <- function(model, overrides = list()) {
  if (!inherits(model, "seirstrata_model")) {
    stop("'model' must be created with build_model()")
  }
  if (length(overrides) == 0) {
    return(model$pars)
  }
  if (!is.list(overrides) || is.null(names(overrides)) ||
      any(!nzchar(names(overrides)))) {
    stop("'overrides' must be a named list")
  }
  unknown <- setdiff(names(overrides), names(model$pars))
  if (length(unknown) > 0) {
    stop(sprintf("Unknown model parameter(s): %s", paste(unknown, collapse = ", ")))
  }
  pars <- model$pars
  for (nm in names(overrides)) {
    old <- pars[[nm]]
    new <- overrides[[nm]]
    same_size <- if (is.null(dim(old))) length(new) == length(old) else
      identical(as.integer(dim(new)), as.integer(dim(old)))
    if (!same_size) {
      stop(sprintf("'%s' must keep its size (%s)", nm,
                   paste(dim(old) %||% length(old), collapse = " x ")))
    }
    pars[[nm]] <- new
  }
  pars
}

#' Likelihood of the model given admissions data
#'
#' Builds a monty model (log-likelihood) from the model's comparison with
#' observed admissions. Add a prior (e.g. with [monty::monty_dsl()]) and use
#' monty to optimise or sample.
#'
#' The fitted parameters `pars` are mapped onto the model's parameters in one
#' of two ways:
#'
#' * Without `process`, each fitted parameter must be a scalar parameter of
#'   the model (e.g. `beta`, `nb_size`) and is used directly.
#' * With `process`, a function of a named list of the fitted values that
#'   returns a named list of model parameters to replace (see [model_pars()]),
#'   for example `function(x) list(beta = base_beta * x$scale)`. It runs once
#'   per likelihood evaluation, so keep it cheap; it must always return the
#'   same names, with the same sizes. `initial` is needed to find out which
#'   names it returns.
#'
#' Everything not set by the fitted parameters stays as in `model`.
#'
#' @param model A model from [build_model()].
#' @param data Admissions data: a data.frame (passed to [model_data()] with
#'   `period` and `age_groups`), or the result of [model_data()].
#' @param pars Names of the fitted parameters.
#' @param process Optional function mapping fitted values to model
#'   parameters, see Details.
#' @param initial Named starting values of `pars`. Required with `process`.
#' @param type `"deterministic"` (the ODE model) or `"stochastic"` (particle
#'   filter on the stochastic model).
#' @param nb_size Negative binomial size (dispersion) used unless it is
#'   fitted: larger means less noise.
#' @param period,age_groups Passed to [model_data()] when `data` is a
#'   data.frame.
#' @param n_particles,dt,n_threads Particle filter settings for
#'   `type = "stochastic"`.
#'
#' @return An object of class `seirstrata_likelihood`, with
#'   * `likelihood`: the monty model;
#'   * `packer`: the monty packer from fitted values to model parameters;
#'   * `filter`: the dust2 filter or unfilter;
#'   * `data`, `model`, `pars`, `type`.
#'
#'   Use [likelihood_pars()] to get the model parameters at a fitted value,
#'   e.g. to run a projection with [run_model()].
#' @export
model_likelihood <- function(model, data, pars, process = NULL, initial = NULL,
                             type = c("deterministic", "stochastic"),
                             nb_size = 10, period = c("week", "day"),
                             age_groups = NULL, n_particles = 200, dt = 0.25,
                             n_threads = 1) {
  if (!inherits(model, "seirstrata_model")) {
    stop("'model' must be created with build_model()")
  }
  type <- match.arg(type)
  if (!inherits(data, "seirstrata_data")) {
    data <- model_data(data, model, period = match.arg(period),
                       age_groups = age_groups)
  } else if (!identical(data$dims, model$dims)) {
    stop("'data' was prepared for a model with different groups")
  }
  if (!is.character(pars) || length(pars) == 0 || anyNA(pars) ||
      anyDuplicated(pars)) {
    stop("'pars' must be the unique names of the fitted parameters")
  }
  check_range(nb_size, "nb_size", min_open = TRUE)

  full <- model$pars
  full$n_obs <- nrow(data$obs_map)
  full$obs_map <- data$obs_map
  full$obs_weekly <- as.numeric(data$period == "week")
  full$nb_size <- nb_size

  if (is.null(process)) {
    bad <- setdiff(pars, names(full))
    if (length(bad) > 0) {
      stop(sprintf(paste("Without 'process', 'pars' must be model parameters;",
                         "not: %s"), paste(bad, collapse = ", ")))
    }
    not_scalar <- pars[lengths(full[pars]) != 1]
    if (length(not_scalar) > 0) {
      stop(sprintf("Fitted parameters must be scalars: %s",
                   paste(not_scalar, collapse = ", ")))
    }
    set_by_process <- character(0)
  } else {
    if (!is.function(process)) {
      stop("'process' must be a function")
    }
    if (is.null(initial) || !setequal(names(initial), pars)) {
      stop("'initial' must give a value for each of 'pars'")
    }
    probe <- process(as.list(initial[pars]))
    if (!is.list(probe) || is.null(names(probe))) {
      stop("'process' must return a named list")
    }
    # checks names and sizes
    model_pars(structure(list(pars = full), class = "seirstrata_model"), probe)
    set_by_process <- names(probe)
    clash <- intersect(set_by_process, pars)
    if (length(clash) > 0) {
      stop(sprintf(paste("'process' must not return the fitted parameters",
                         "themselves: %s"), paste(clash, collapse = ", ")))
    }
  }
  fixed <- full[setdiff(names(full), c(pars, set_by_process))]
  packer <- monty::monty_packer(pars, fixed = fixed, process = process)

  if (type == "deterministic") {
    filter <- dust2::dust_unfilter_create(seir_ode(), 0, data$data_dust)
  } else {
    filter <- dust2::dust_filter_create(seir_stoch(), 0, data$data_dust,
                                        n_particles = n_particles, dt = dt,
                                        n_threads = n_threads)
  }
  structure(list(likelihood = dust2::dust_likelihood_monty(filter, packer),
                 packer = packer, filter = filter, data = data, model = model,
                 pars = pars, type = type,
                 sets = setdiff(c(intersect(pars, names(full)), set_by_process),
                                c("n_obs", "obs_map", "obs_weekly", "nb_size"))),
            class = "seirstrata_likelihood")
}

#' Model parameters at a fitted value
#'
#' The model parameters that the fitted parameters set, at the value `x`:
#' pass them to [run_model()] as `pars`, e.g. to project over posterior draws.
#' Observation-model parameters (`nb_size`) are left out.
#'
#' @param likelihood Result of [model_likelihood()].
#' @param x Values of the fitted parameters, in the order of
#'   `likelihood$pars` (or named).
#' @return A named list of model parameters, for [model_pars()].
#' @export
likelihood_pars <- function(likelihood, x) {
  if (!is.null(names(x))) {
    x <- x[likelihood$pars]
  }
  p <- likelihood$packer$unpack(unname(x))
  p[likelihood$sets]
}

#' @export
print.seirstrata_likelihood <- function(x, ...) {
  cat(sprintf("<seirstrata_likelihood> %s, fitting: %s\n", x$type,
              paste(x$pars, collapse = ", ")))
  print(x$data)
  invisible(x)
}
