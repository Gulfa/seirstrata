#' Run a model
#'
#' Runs the deterministic (ODE) or stochastic (tau-leaping) version of the
#' model and returns daily output from time 0 to `t_max`.
#'
#' @param model A model from [build_model()].
#' @param type `"deterministic"` or `"stochastic"`.
#' @param t_max Number of days to simulate.
#' @param n_sim Number of stochastic simulations (ignored for deterministic).
#' @param dt Time step (days) of the stochastic model; must divide one day
#'   evenly. The discrete-time approximation lengthens every stage by about
#'   `dt / 2` on average, which makes epidemics slightly larger than in the
#'   deterministic model (about +0.7% infections at `dt = 0.1` and +1.7% at
#'   `0.25` for R0 = 2). Smaller is more accurate but slower.
#' @param seed Random seed for the stochastic model (`NULL` = random).
#' @param n_threads Number of threads used to run stochastic simulations in
#'   parallel.
#' @param vars Output variables to keep (default: all, see [get_output()]).
#'   Keeping fewer saves memory for large stochastic runs; `summary()` needs
#'   `cum_infections`, `cum_hosp`, `D` and `H`.
#' @param state Optional state to start from, as returned by dust2 for this
#'   model, instead of the model's initial conditions, e.g. to continue
#'   from a particle filter's state.
#' @param t_start First day to simulate, when continuing from `state`.
#'
#' @return An object of class `seirstrata_run`. Use [get_output()] to extract
#'   results as a table, and `summary()` / `plot()` to look at them.
#' @export
#' @examples
#' m <- build_model(N = 1e5, R0 = 2, I0 = 10, p_hosp = 0.02)
#' res <- run_model(m, t_max = 200)
#' res_s <- run_model(m, "stochastic", t_max = 200, n_sim = 20, seed = 1)
run_model <- function(model, type = c("deterministic", "stochastic"),
                      t_max = 365, n_sim = 100, dt = 0.1, seed = NULL,
                      n_threads = 1, vars = NULL, state = NULL,
                      t_start = 0) {
  if (!inherits(model, "seirstrata_model")) {
    stop("'model' must be created with build_model()")
  }
  type <- match.arg(type)
  check_range(t_max, "t_max", min_open = TRUE)
  check_range(t_start, "t_start")
  if (t_start >= t_max) {
    stop("'t_start' must be before 't_max'")
  }
  times <- seq(ceiling(t_start), ceiling(t_max))
  if (!is.na(model$beta_t_end) && t_max > model$beta_t_end) {
    warning(sprintf(paste(
      "'beta_t' covers days 0 to %g, but the run goes to day %g:",
      "its last value is held after that"), model$beta_t_end, t_max))
  }

  if (type == "deterministic") {
    n_sim <- 1L
    sys <- dust2::dust_system_create(seir_ode(), model$pars)
  } else {
    check_range(n_sim, "n_sim", min = 1)
    check_range(dt, "dt", min_open = TRUE, max = 1)
    if (abs(1 / dt - round(1 / dt)) > 1e-8) {
      stop("'dt' must divide one day evenly (e.g. 1, 0.5, 0.25, 0.1)")
    }
    n_sim <- as.integer(n_sim)
    sys <- dust2::dust_system_create(seir_stoch(), model$pars,
                                     n_particles = n_sim, dt = dt,
                                     seed = seed, n_threads = n_threads)
  }
  inherited <- NULL
  if (is.null(state)) {
    dust2::dust_system_set_state_initial(sys)
  } else {
    # Continue from a state reached earlier (e.g. a particle filter's state at
    # the end of the data). The daily incidence counters hold that day's
    # partial count; left in place they would be added to the first simulated
    # day, so zero them here and put them back afterwards as the first day's
    # output.
    index <- dust2::dust_unpack_state(sys, seq_along(state))
    daily <- grep("^inc_", names(index), value = TRUE)
    inherited <- lapply(index[daily], function(i) state[i])
    state[unlist(index[daily])] <- 0
    dust2::dust_system_set_state(sys, state)
  }
  dust2::dust_system_set_time(sys, times[[1]])
  y <- dust2::dust_system_simulate(sys, times)
  raw <- dust2::dust_unpack_state(sys, y)

  # [group, stratum, (sim), time] -> [age, region, stratum, sim, time].
  # Variables without strata (R and doses) get a single stratum "all".
  dn <- list(age = model$dims$age, region = model$dims$region,
             stratum = model$dims$stratum, sim = seq_len(n_sim), time = times)
  dn1 <- dn
  dn1$stratum <- "all"
  n_stratified <- prod(lengths(dn, use.names = FALSE))
  states <- lapply(raw, function(x) {
    if (length(x) == n_stratified) {
      array(x, lengths(dn, use.names = FALSE), dn)
    } else {
      array(x, lengths(dn1, use.names = FALSE), dn1)
    }
  })
  for (v in names(inherited)) {
    states[[v]][, , , , 1] <- rep(inherited[[v]], times = n_sim)
  }
  states$D <- states$D_hosp + states$D_comm
  states$inc_deaths <- states$inc_deaths_hosp + states$inc_deaths_comm
  if (!is.null(vars)) {
    unknown <- setdiff(vars, names(states))
    if (length(unknown) > 0) {
      stop(sprintf("Unknown variable(s): %s. Available: %s",
                   paste(unknown, collapse = ", "),
                   paste(names(states), collapse = ", ")))
    }
    states <- states[vars]
  }

  structure(list(model = model, type = type, time = times, n_sim = n_sim,
                 dt = if (type == "stochastic") dt else NA_real_,
                 states = states),
            class = "seirstrata_run")
}

#' @export
print.seirstrata_run <- function(x, ...) {
  cat(sprintf("<seirstrata_run> %s, %d day(s)%s\n", x$type, max(x$time),
              if (x$type == "stochastic")
                sprintf(", %d simulation(s), dt = %g", x$n_sim, x$dt) else ""))
  cat(sprintf("  Variables: %s\n", paste(names(x$states), collapse = ", ")))
  cat("  Use get_output(), summary() and plot() to look at the results\n")
  invisible(x)
}
