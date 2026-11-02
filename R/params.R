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
  build_group_schedule(schedule, dims, "doses", "vacc_schedule", default = 0,
                       row_value = function(value, match) {
                         value * N[match] / sum(N[match])
                       })
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
