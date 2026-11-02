#' Combine an age mixing matrix and a region mixing matrix
#'
#' The backend model works on `n_age * n_region` groups, numbered with age
#' varying fastest: group `(r - 1) * n_age + a` is age `a` in region `r`. The
#' contact rate between age `a` in region `r` and age `b` in region `s` is
#'
#' `C[(a, r), (b, s)] = C_age[a, b] * M_region[r, s]`
#'
#' i.e. everyone has the age-mixing pattern `C_age`, and `M_region` decides in
#' which regions those contacts take place.
#'
#' @param C_age `n_age x n_age` matrix; `C_age[a, b]` is the number of contacts
#'   per day a person of age `a` has with people of age `b`.
#' @param M_region `n_region x n_region` matrix; `M_region[r, s]` is the
#'   fraction of a region-`r` resident's contacts made with residents of region
#'   `s`. Rows must sum to 1. The identity matrix gives isolated regions.
#'
#' @return The `n_age * n_region` square contact matrix.
#' @export
combine_mixing <- function(C_age, M_region) {
  kronecker(M_region, C_age)
}

check_C_age <- function(C_age, n_age, age_labels) {
  if (is.null(C_age)) {
    if (n_age > 1) {
      stop("'C_age' is required when there is more than one age group")
    }
    return(matrix(1, 1, 1, dimnames = list(age_labels, age_labels)))
  }
  C_age <- as.matrix(C_age)
  if (!identical(dim(C_age), c(n_age, n_age))) {
    stop(sprintf("'C_age' must be a %d x %d matrix (one row/column per age group)",
                 n_age, n_age))
  }
  if (anyNA(C_age) || any(C_age < 0)) {
    stop("'C_age' must be non-negative")
  }
  dimnames(C_age) <- list(age_labels, age_labels)
  C_age
}

check_M_region <- function(M_region, n_region, region_labels) {
  if (is.null(M_region)) {
    if (n_region > 1) {
      stop("'M_region' is required when there is more than one region")
    }
    return(matrix(1, 1, 1, dimnames = list(region_labels, region_labels)))
  }
  M_region <- as.matrix(M_region)
  if (!identical(dim(M_region), c(n_region, n_region))) {
    stop(sprintf("'M_region' must be a %d x %d matrix (one row/column per region)",
                 n_region, n_region))
  }
  if (anyNA(M_region) || any(M_region < 0)) {
    stop("'M_region' must be non-negative")
  }
  if (any(abs(rowSums(M_region) - 1) > 1e-6)) {
    stop("Rows of 'M_region' must sum to 1 (fraction of contacts made in each region)")
  }
  dimnames(M_region) <- list(region_labels, region_labels)
  M_region
}

# Build the piecewise-constant contact matrix over time.
#
# `C_age` and `M_region` give the mixing from time 0. Each entry of
# `schedule` is a list with `time` and `C_age` and/or `M_region`; from that
# time on, the given matrix replaces the one in use (the other one carries
# over). Returns list(time = <n_t>, value = <n x n x n_t array>).
build_mixing_schedule <- function(C_age, M_region, schedule, dims) {
  combined <- function(C, M) combine_mixing(C, M)
  if (is.null(schedule)) {
    return(list(time = 0, value = array(combined(C_age, M_region),
                                        c(rep(length(dims$age) * length(dims$region), 2), 1))))
  }
  if (!is.list(schedule) || !all(vapply(schedule, is.list, TRUE))) {
    stop("'mixing_schedule' must be a list of lists, each with 'time' and ",
         "'C_age' and/or 'M_region'")
  }
  times <- vapply(schedule, function(e) {
    if (!is.numeric(e$time) || length(e$time) != 1 || is.na(e$time)) {
      stop("Each entry of 'mixing_schedule' needs a single 'time'")
    }
    e$time
  }, 1)
  check_range(times, "mixing_schedule time")
  for (e in schedule) {
    unknown <- setdiff(names(e), c("time", "C_age", "M_region"))
    if (length(unknown) > 0) {
      stop(sprintf("Unknown element(s) in 'mixing_schedule': %s",
                   paste(unknown, collapse = ", ")))
    }
    if (is.null(e$C_age) && is.null(e$M_region)) {
      stop("Each entry of 'mixing_schedule' must change 'C_age' or 'M_region'")
    }
  }

  schedule <- schedule[order(times)]
  times <- sort(times)
  if (times[1] > 0) {
    times <- c(0, times)
    schedule <- c(list(list(time = 0)), schedule)
  }
  n <- length(dims$age) * length(dims$region)
  value <- array(NA_real_, c(n, n, length(times)))
  current_C <- C_age
  current_M <- M_region
  for (k in seq_along(times)) {
    e <- schedule[[k]]
    if (!is.null(e$C_age)) {
      current_C <- check_C_age(e$C_age, length(dims$age), dims$age)
    }
    if (!is.null(e$M_region)) {
      current_M <- check_M_region(e$M_region, length(dims$region), dims$region)
    }
    value[, , k] <- combined(current_C, current_M)
  }
  list(time = times, value = value)
}
