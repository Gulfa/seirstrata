# Expanding user-supplied parameters to the n_age x n_region group structure.

# A parameter can be given as
#   * a single value (same for every group),
#   * a vector by age (length n_age),
#   * a vector by region (length n_region, or named with the region labels),
#   * an n_age x n_region matrix.
# Returns an n_age x n_region matrix. as.vector() of the result is in backend
# group order (age varying fastest).
expand_param <- function(x, name, dims) {
  n_age <- length(dims$age)
  n_region <- length(dims$region)
  if (!is.numeric(x) || anyNA(x)) {
    stop(sprintf("'%s' must be numeric without missing values", name))
  }
  if (is.matrix(x)) {
    if (!identical(dim(x), c(n_age, n_region))) {
      stop(sprintf("'%s' must be a %d x %d (age x region) matrix", name,
                   n_age, n_region))
    }
    ret <- x
  } else if (length(x) == 1) {
    ret <- matrix(x, n_age, n_region)
  } else if (!is.null(names(x)) && setequal(names(x), dims$region) &&
             !setequal(names(x), dims$age)) {
    ret <- matrix(x[dims$region], n_age, n_region, byrow = TRUE)
  } else if (length(x) == n_age) {
    if (!is.null(names(x)) && setequal(names(x), dims$age)) {
      x <- x[dims$age]
    }
    ret <- matrix(x, n_age, n_region)
  } else if (length(x) == n_region) {
    ret <- matrix(x, n_age, n_region, byrow = TRUE)
  } else {
    stop(sprintf(paste(
      "'%s' must be a single value, a vector by age (length %d), a vector by",
      "region (length %d) or an age x region matrix"), name, n_age, n_region))
  }
  dimnames(ret) <- list(dims$age, dims$region)
  ret
}

check_range <- function(x, name, min = 0, max = Inf, min_open = FALSE) {
  bad <- if (min_open) x <= min else x < min
  if (any(bad) || any(x > max)) {
    range <- sprintf("%s%s, %s]", if (min_open) "(" else "[", min, max)
    stop(sprintf("'%s' must be in %s", name, range))
  }
  invisible(x)
}

# Population: a scalar (one group), a vector by age (one region), or an
# n_age x n_region matrix. Labels come from names / dimnames.
parse_population <- function(N) {
  if (is.data.frame(N)) {
    stop("'N' must be a number, a named vector by age or an age x region matrix")
  }
  if (is.matrix(N)) {
    age <- rownames(N)
    region <- colnames(N)
  } else {
    age <- names(N)
    region <- NULL
    N <- matrix(N, ncol = 1)
  }
  if (is.null(age)) {
    age <- if (nrow(N) == 1) "all" else paste0("age", seq_len(nrow(N)))
  }
  if (is.null(region)) {
    region <- if (ncol(N) == 1) "all" else paste0("region", seq_len(ncol(N)))
  }
  if (anyDuplicated(age) || anyDuplicated(region)) {
    stop("Age group and region labels must be unique")
  }
  if (!is.numeric(N) || anyNA(N) || any(N <= 0)) {
    stop("'N' must be positive in every age group and region")
  }
  if (any(abs(N - round(N)) > 1e-8)) {
    stop("'N' must contain whole numbers of people")
  }
  dimnames(N) <- list(age, region)
  list(N = round(N), dims = list(age = age, region = region))
}

# Build the per-group piecewise-constant contact multiplier.
#
# `schedule` is a data.frame with columns `time` and `mult`, and optionally
# `age` and/or `region` (NA or missing = all). Rows are applied in time order:
# from `time` onwards, the matching groups get multiplier `mult`, until a later
# row matching them overrides it. Before any row applies the multiplier is 1.
#
# Returns list(time = <n_t>, value = <n_group x n_t matrix>).
build_beta_schedule <- function(schedule, dims) {
  build_group_schedule(schedule, dims, "mult", "beta_schedule", default = 1,
                       row_value = function(value, match) value)
}

# Build the per-group piecewise-constant vaccination rate (doses per day).
#
# Same rules as build_beta_schedule(), with a column `doses`: from `time`
# onwards, the matching groups share `doses` per day in proportion to their
# population. Before any row applies, no doses are given.
build_vacc_schedule <- function(schedule, dims, N) {
  N <- as.vector(N)
  steps <- build_group_schedule(schedule, dims, "doses", "vacc_schedule",
                                default = 0, row_value = function(value, match) {
                                  value * N[match] / sum(N[match])
                                })
  ret <- ramp_steps(steps$time, steps$value)
  ret$steps <- steps
  ret
}

# Turn a piecewise-constant schedule (value[, k] from time[k] onwards) into
# knots for linear interpolation: each change becomes a linear ramp over one
# day centred on its time (narrower if changes are closer than a day). The
# area under a symmetric ramp equals that under the step, so totals (e.g.
# doses) are exact, while the input stays continuous: steps make the ODE
# solver stop at every change and slow it down several times. A final knot
# far in the future holds the last value.
ramp_steps <- function(time, value, width = 1) {
  n <- length(time)
  value <- as.matrix(value)
  if (n == 1) {
    return(list(time = c(0, 1e6), value = cbind(value, value)))
  }
  gaps <- diff(time)
  # half-width at each change: at most width / 2 and half the gaps around it
  # (change k is at time[k + 1]: gaps[k] is the gap before it, gaps[k + 1]
  # the gap after it)
  half <- pmin(width / 2, gaps / 2, c(gaps[-1], Inf) / 2)
  t_knot <- c(time[1], as.vector(rbind(time[-1] - half, time[-1] + half)),
              time[n] + 1e6)
  idx <- c(1, as.vector(rbind(seq_len(n - 1), seq_len(n - 1) + 1)), n)
  # a change at time 0 cannot ramp from before the start
  keep <- t_knot >= time[1]
  t_knot <- t_knot[keep]
  idx <- idx[keep]
  # drop knots at the same time (zero-width ramps), keeping the later one
  dup <- c(diff(t_knot) == 0, FALSE)
  list(time = t_knot[!dup], value = value[, idx[!dup], drop = FALSE])
}

# The common part of the schedules: a data.frame with columns `time` and
# `value_col`, and optionally `age` and/or `region` (NA or missing = all).
# Rows are applied in time order; from `time` onwards the groups a row matches
# take `row_value(value, match)` (one number per matched group), until a later
# row matching them overrides it.
build_group_schedule <- function(schedule, dims, value_col, name, default,
                                 row_value) {
  n_age <- length(dims$age)
  n_region <- length(dims$region)
  n_group <- n_age * n_region
  if (is.null(schedule)) {
    return(list(time = 0, value = matrix(default, n_group, 1)))
  }
  schedule <- as.data.frame(schedule, stringsAsFactors = FALSE)
  if (!all(c("time", value_col) %in% names(schedule))) {
    stop(sprintf("'%s' must have columns 'time' and '%s'", name, value_col))
  }
  extra <- setdiff(names(schedule), c("time", value_col, "age", "region"))
  if (length(extra) > 0) {
    stop(sprintf("Unknown column(s) in '%s': %s", name,
                 paste(extra, collapse = ", ")))
  }
  check_range(schedule$time, paste0(name, "$time"))
  check_range(schedule[[value_col]], paste0(name, "$", value_col))
  for (v in c("age", "region")) {
    if (is.null(schedule[[v]])) {
      schedule[[v]] <- NA_character_
    }
    schedule[[v]] <- as.character(schedule[[v]])
    bad <- setdiff(stats::na.omit(schedule[[v]]), dims[[v]])
    if (length(bad) > 0) {
      stop(sprintf("Unknown %s in '%s': %s (known: %s)", v, name,
                   paste(bad, collapse = ", "), paste(dims[[v]], collapse = ", ")))
    }
  }

  group_age <- rep(dims$age, n_region)
  group_region <- rep(dims$region, each = n_age)
  schedule <- schedule[order(schedule$time), , drop = FALSE]
  times <- sort(unique(c(0, schedule$time)))
  value <- matrix(NA_real_, n_group, length(times))
  current <- rep(default, n_group)
  for (k in seq_along(times)) {
    rows <- schedule[schedule$time == times[k], , drop = FALSE]
    for (r in seq_len(nrow(rows))) {
      match <- (is.na(rows$age[r]) | group_age == rows$age[r]) &
        (is.na(rows$region[r]) | group_region == rows$region[r])
      current[match] <- row_value(rows[[value_col]][r], match)
    }
    value[, k] <- current
  }
  list(time = times, value = value)
}

# Immunity strata --------------------------------------------------------------

check_strata <- function(strata) {
  strata <- as.character(strata)
  if (length(strata) == 0 || anyNA(strata) || any(!nzchar(strata)) ||
      anyDuplicated(strata)) {
    stop("'strata' must be unique, non-empty labels")
  }
  if ("all" %in% strata) {
    stop("'all' is reserved and cannot be a stratum label")
  }
  strata
}

# A stratum parameter can be given as
#   * a single value (same for every age group and stratum),
#   * a vector by stratum (length n_imm),
#   * an n_age x n_imm matrix.
# Returns an n_age x n_imm matrix.
expand_stratum_param <- function(x, name, dims) {
  n_age <- length(dims$age)
  n_imm <- length(dims$stratum)
  if (!is.numeric(x) || anyNA(x)) {
    stop(sprintf("'%s' must be numeric without missing values", name))
  }
  if (is.matrix(x)) {
    if (!identical(dim(x), c(n_age, n_imm))) {
      stop(sprintf("'%s' must be a %d x %d (age x stratum) matrix", name,
                   n_age, n_imm))
    }
    ret <- x
  } else if (length(x) == 1) {
    ret <- matrix(x, n_age, n_imm)
  } else if (length(x) == n_imm) {
    if (!is.null(names(x)) && setequal(names(x), dims$stratum)) {
      x <- x[dims$stratum]
    }
    ret <- matrix(x, n_age, n_imm, byrow = TRUE)
  } else {
    stop(sprintf(paste(
      "'%s' must be a single value, a vector by stratum (length %d) or an",
      "age x stratum matrix"), name, n_imm))
  }
  dimnames(ret) <- list(dims$age, dims$stratum)
  ret
}

# A stratum given as a label or a number -> its index
stratum_index <- function(x, name, strata) {
  if (length(x) != 1 || is.na(x)) {
    stop(sprintf("'%s' must be a single stratum", name))
  }
  k <- if (is.character(x)) match(x, strata) else x
  if (is.na(k) || k < 1 || k > length(strata) || k != round(k)) {
    stop(sprintf("'%s' must be one of: %s", name,
                 paste(strata, collapse = ", ")))
  }
  as.integer(k)
}

# Global transmission multiplier ------------------------------------------------

# beta_t as a numeric vector of daily values (knots at days 0, 1, ...) or a
# data.frame with columns time and mult (knots at `time`, starting at 0). The
# model interpolates linearly between knots; a final knot far in the future
# makes the last value hold. Returns list(time, value, end), where end is the
# last day of a daily vector (NA for a data.frame, whose last value is meant
# to hold).
build_beta_t <- function(beta_t) {
  if (is.null(beta_t)) {
    time <- 0
    value <- 1
    end <- NA_real_
  } else if (is.data.frame(beta_t)) {
    if (!all(c("time", "mult") %in% names(beta_t))) {
      stop("'beta_t' as a data.frame must have columns 'time' and 'mult'")
    }
    extra <- setdiff(names(beta_t), c("time", "mult"))
    if (length(extra) > 0) {
      stop(sprintf("Unknown column(s) in 'beta_t': %s",
                   paste(extra, collapse = ", ")))
    }
    time <- beta_t$time
    value <- beta_t$mult
    if (!is.numeric(time) || anyNA(time) ||
        is.unsorted(time, strictly = TRUE)) {
      stop("'beta_t$time' must be strictly increasing")
    }
    if (time[1] != 0) {
      stop("'beta_t$time' must start at 0")
    }
    end <- NA_real_
  } else {
    if (!is.numeric(beta_t) || length(beta_t) == 0) {
      stop("'beta_t' must be a numeric vector of daily values or a data.frame")
    }
    value <- as.vector(beta_t)
    time <- seq_along(value) - 1
    end <- max(time)
  }
  if (!is.numeric(value) || anyNA(value) || any(value < 0)) {
    stop("'beta_t' values must be non-negative numbers")
  }
  hold <- max(time) + 1e6
  list(time = as.numeric(c(time, hold)),
       value = as.numeric(c(value, value[length(value)])), end = end)
}

# beta_t at the given times, as the model interpolates it
beta_t_at <- function(model, time) {
  p <- model$pars
  stats::approx(p$beta_t_time, p$beta_t_value, xout = time, rule = 2)$y
}
