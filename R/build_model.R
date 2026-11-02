#' Build a model
#'
#' Collects, checks and expands all parameters of the SEIR-HD model with an
#' immunity ladder. Nothing is simulated; pass the result to [run_model()].
#'
#' The population is stratified by age group and region. Group parameters can
#' be given as a single value, a vector by age group, a vector by region
#' (length `n_region`, or named with the region labels) or an
#' `n_age x n_region` matrix.
#'
#' @section Immunity:
#' Susceptible people sit on a ladder of immunity strata, from stratum 1 (no
#' immunity) to stratum `n_imm = length(strata)`. Stratum parameters
#' (`rel_susceptibility`, `rel_transmissibility`, `rel_hosp`, `rel_death`,
#' `waning_rate`, `initial_strata`) can be a single value, a vector by stratum
#' (length `n_imm`), or an `n_age x n_imm` matrix. A vector whose length
#' matches both the number of age groups and of strata is taken to be by
#' stratum.
#'
#' * `rel_susceptibility`: multiplies `susceptibility`.
#' * `rel_transmissibility`: multiplies the infectiousness of people infected
#'   in that stratum.
#' * `rel_hosp`: multiplies `p_hosp`. This is protection *given infection*;
#'   total protection against hospitalisation is
#'   `1 - rel_susceptibility * rel_hosp`.
#' * `rel_death`: multiplies `p_death_hosp` and `p_death_comm`.
#'
#' An infection keeps the stratum of the person infected, so these multipliers
#' apply throughout the infection. Everyone who recovers goes to a single
#' recovered compartment `R`, which is fully protected. People leave `R` at
#' rate `waning_rate_R` and enter the ladder at `R_to_stratum`. On the ladder,
#' people move from stratum k to k - 1 at rate `waning_rate[k]`.
#'
#' Vaccination status is not tracked. Doses from `vacc_schedule` go to people
#' in strata below `vacc_stratum` and to people in `R`, in proportion. Doses to
#' `R` are counted but change nothing.
#'
#' Initially a fraction `initial_R` of each group is in `R`. The rest are spread
#' over the strata by `initial_strata`. `I0` and `E0` are placed in stratum 1
#' and taken from its share.
#'
#' @param N Population: a single number, a (named) vector by age group, or an
#'   `n_age x n_region` matrix whose dimnames give the age group and region
#'   labels. A region-only model is a `1 x n_region` matrix.
#' @param C_age Age contact matrix, `C_age[a, b]` = contacts per day an
#'   age-`a` person has with age-`b` people. Required with more than one age
#'   group. Only its shape matters when `R0` is given, not its scale.
#' @param M_region Region mixing matrix, `M_region[r, s]` = fraction of a
#'   region-`r` resident's contacts made with region-`s` residents (rows sum to
#'   1). Required with more than one region. See [combine_mixing()].
#' @param R0,beta Give exactly one. `R0` calibrates `beta` with the
#'   next-generation matrix for a population entirely in stratum 1 (see
#'   [R0()]).
#' @param latent_period,infectious_period Mean durations (days) of E and I.
#' @param susceptibility Relative susceptibility, by group.
#' @param transmissibility Relative infectiousness when infected, by group.
#' @param p_hosp Probability that an infection in stratum 1 is hospitalised.
#' @param p_death_hosp Probability that a hospitalised case in stratum 1 dies.
#' @param p_death_comm Probability that an infection in stratum 1 that is not
#'   hospitalised dies (in the community).
#' @param dur_hosp Mean length of hospital stay (days), by group.
#' @param strata Labels of the immunity strata, from no immunity upwards. Their
#'   number sets `n_imm`. The default is a single stratum, i.e. an SEIR model
#'   where `R` is the only immunity.
#' @param rel_susceptibility,rel_transmissibility,rel_hosp,rel_death Stratum
#'   multipliers, see the Immunity section.
#' @param waning_rate Rate (per day) of moving from stratum k down to k - 1.
#'   A single value applies to strata 2 and up; as a vector or matrix, the
#'   entry for stratum 1 must be 0.
#' @param waning_rate_R Rate (per day) of leaving `R`, by group. The default 0
#'   means immunity from infection lasts for the whole run.
#' @param R_to_stratum Stratum (label or number) entered on leaving `R`.
#'   Defaults to the top stratum.
#' @param initial_strata Initial fractions over the strata of the population
#'   not in `R`. A vector by stratum, or an `n_age x n_imm` matrix, whose rows
#'   sum to 1. Default: everyone in stratum 1.
#' @param initial_R Fraction of each group initially in `R`.
#' @param vacc_stratum Stratum (label or number) that vaccination moves people
#'   to, from the strata below it. Defaults to the top stratum.
#' @param vacc_schedule Optional vaccination during the epidemic: a data.frame
#'   with columns `time` and `doses` (doses per day), and optionally `age`
#'   and/or `region` to target groups (`NA` = all). Rows work as in
#'   `beta_schedule`: from `time` onwards the matching groups receive `doses`
#'   per day between them, shared in proportion to their population, until a
#'   later row matching them overrides it (`doses = 0` ends a campaign). When a
#'   group's eligible people run out, its vaccination stops.
#' @param I0,E0 Initial number of infectious / exposed people (placed in
#'   stratum 1).
#' @param mixing_schedule Optional changes of the mixing matrices over time: a
#'   list of entries like `list(time = 40, C_age = C_holiday)`, each replacing
#'   `C_age`, `M_region` or both from that time on (the one not given carries
#'   over). `C_age` and `M_region` themselves apply from time 0. `R0` is
#'   calibrated from the mixing at time 0.
#' @param beta_schedule Optional contact multipliers over time, e.g. for
#'   interventions. Either a data.frame with columns `time` and `mult` (applies
#'   to everyone), or with extra columns `age` and/or `region` to target
#'   groups (`NA` = all). Rows are applied in time order: from `time` onwards
#'   the matching groups get multiplier `mult` until a later row matching them
#'   overrides it. `mult` is the fraction of the group's contacts that still
#'   happen: contacts between groups i and j are scaled by
#'   `min(mult_i, mult_j)`.
#' @param beta_t Optional global transmission multiplier over time, for
#'   seasonality, holidays and the like, typically built outside the package.
#'   It multiplies `beta` for everyone, on top of `beta_schedule`. Either a
#'   numeric vector of daily values (`beta_t[1]` at day 0, `beta_t[2]` at day
#'   1, ...), or a data.frame with columns `time` (starting at 0) and `mult`.
#'   The values are knots: the multiplier is linear between them, so a sharp
#'   step is two knots a day apart. (Linear rather than piecewise constant,
#'   because daily steps make the ODE solver several times slower.) After the
#'   last knot it keeps the last value: [run_model()] warns when a run goes
#'   beyond the end of a daily vector. `R0` does not include `beta_t`;
#'   `R_eff(model, time)` does.
#'
#' @return An object of class `seirstrata_model`.
#' @export
#' @examples
#' m <- build_model(N = 1e5, R0 = 2, I0 = 10, p_hosp = 0.02)
#' m
#'
#' # Three strata: none, partial and good protection, waning over ~6 months
#' m3 <- build_model(N = 1e5, R0 = 2, I0 = 10, p_hosp = 0.02,
#'                   strata = c("none", "partial", "good"),
#'                   rel_susceptibility = c(1, 0.7, 0.4),
#'                   rel_hosp = c(1, 0.6, 0.3),
#'                   waning_rate = 1 / 90,
#'                   initial_strata = c(0.5, 0.3, 0.2))
build_model <- function(N,
                        C_age = NULL,
                        M_region = NULL,
                        R0 = NULL,
                        beta = NULL,
                        latent_period = 3,
                        infectious_period = 5,
                        susceptibility = 1,
                        transmissibility = 1,
                        p_hosp = 0,
                        p_death_hosp = 0,
                        p_death_comm = 0,
                        dur_hosp = 7,
                        strata = "1",
                        rel_susceptibility = 1,
                        rel_transmissibility = 1,
                        rel_hosp = 1,
                        rel_death = 1,
                        waning_rate = 0,
                        waning_rate_R = 0,
                        R_to_stratum = NULL,
                        initial_strata = NULL,
                        initial_R = 0,
                        vacc_stratum = NULL,
                        I0 = 0,
                        E0 = 0,
                        beta_schedule = NULL,
                        mixing_schedule = NULL,
                        vacc_schedule = NULL,
                        beta_t = NULL) {
  args <- as.list(environment())

  pop <- parse_population(N)
  N <- pop$N
  dims <- pop$dims
  dims$stratum <- check_strata(strata)
  n_group <- length(N)
  n_imm <- length(dims$stratum)

  C_age <- check_C_age(C_age, length(dims$age), dims$age)
  M_region <- check_M_region(M_region, length(dims$region), dims$region)
  mixing <- build_mixing_schedule(C_age, M_region, mixing_schedule, dims)
  contact <- mixing$value[, , 1]

  if (is.null(R0) == is.null(beta)) {
    stop("Give exactly one of 'R0' and 'beta'")
  }
  check_range(latent_period, "latent_period", min_open = TRUE)
  check_range(infectious_period, "infectious_period", min_open = TRUE)

  # Everything group-specific as an age x region matrix
  g <- function(x, name, ...) {
    ret <- expand_param(x, name, dims)
    check_range(ret, name, ...)
    ret
  }
  susceptibility <- g(susceptibility, "susceptibility")
  transmissibility <- g(transmissibility, "transmissibility")
  p_hosp <- g(p_hosp, "p_hosp", max = 1)
  p_death_hosp <- g(p_death_hosp, "p_death_hosp", max = 1)
  p_death_comm <- g(p_death_comm, "p_death_comm", max = 1)
  dur_hosp <- g(dur_hosp, "dur_hosp", min_open = TRUE)
  waning_rate_R <- g(waning_rate_R, "waning_rate_R")
  initial_R <- g(initial_R, "initial_R", max = 1)
  I0 <- g(I0, "I0")
  E0 <- g(E0, "E0")
  if (any(abs(c(I0, E0) - round(c(I0, E0))) > 1e-8)) {
    stop("'I0' and 'E0' must be whole numbers of people")
  }

  # Everything stratum-specific as an n_age x n_imm matrix
  s <- function(x, name, ...) {
    ret <- expand_stratum_param(x, name, dims)
    check_range(ret, name, ...)
    ret
  }
  rel_susceptibility <- s(rel_susceptibility, "rel_susceptibility")
  rel_transmissibility <- s(rel_transmissibility, "rel_transmissibility")
  rel_hosp <- s(rel_hosp, "rel_hosp")
  rel_death <- s(rel_death, "rel_death")
  if (length(waning_rate) == 1 && !is.matrix(waning_rate)) {
    waning_rate <- c(0, rep(waning_rate, n_imm - 1))
  }
  waning_rate <- s(waning_rate, "waning_rate")
  if (any(waning_rate[, 1] != 0)) {
    stop("'waning_rate' must be 0 for stratum 1 (there is nothing below it)")
  }
  if (is.null(initial_strata)) {
    initial_strata <- c(1, rep(0, n_imm - 1))
  }
  initial_strata <- s(initial_strata, "initial_strata", max = 1)
  if (any(abs(rowSums(initial_strata) - 1) > 1e-8)) {
    stop("'initial_strata' must sum to 1 over the strata")
  }
  R_to <- stratum_index(R_to_stratum %||% n_imm, "R_to_stratum", dims$stratum)
  vacc_to <- stratum_index(vacc_stratum %||% n_imm, "vacc_stratum",
                           dims$stratum)
  if (!is.null(vacc_schedule) && vacc_to == 1) {
    stop("'vacc_schedule' needs a 'vacc_stratum' above stratum 1 ",
         "(give more than one stratum in 'strata')")
  }

  # [group, stratum] matrices for the backend: group parameter times stratum
  # multiplier. Stratum matrices are by age, repeated over regions.
  by_group <- function(x) x[rep(seq_along(dims$age), length(dims$region)), ,
                            drop = FALSE]
  gs <- function(group, stratum) {
    ret <- as.vector(group) * by_group(stratum)
    dimnames(ret) <- NULL
    ret
  }
  sus <- gs(susceptibility, rel_susceptibility)
  trans <- gs(transmissibility, rel_transmissibility)
  p_hosp_s <- gs(p_hosp, rel_hosp)
  p_death_comm_s <- gs(p_death_comm, rel_death)
  p_death_hosp_s <- gs(p_death_hosp, rel_death)
  for (p in list(list(p_hosp_s, "p_hosp * rel_hosp"),
                 list(p_death_comm_s, "p_death_comm * rel_death"),
                 list(p_death_hosp_s, "p_death_hosp * rel_death"))) {
    check_range(p[[1]], p[[2]], max = 1)
  }

  # Initial state
  N_group <- as.vector(N)
  R_ini <- round(N_group * as.vector(initial_R))
  rest <- N_group - R_ini
  S_ini <- round(rest * by_group(initial_strata))
  S_ini[, 1] <- rest - as.vector(I0) - as.vector(E0) -
    rowSums(S_ini[, -1, drop = FALSE])
  dimnames(S_ini) <- NULL
  if (any(S_ini[, 1] < 0)) {
    stop("I0 + E0 exceeds the people in stratum 1 in some groups")
  }
  seed <- function(x) cbind(as.vector(x), matrix(0, n_group, n_imm - 1))

  sched <- build_beta_schedule(beta_schedule, dims)
  bt <- build_beta_t(beta_t)
  vacc <- build_vacc_schedule(vacc_schedule, dims, N)
  one_hot <- function(k) as.numeric(seq_len(n_imm) == k)

  pars <- list(
    n_group = n_group,
    n_imm = n_imm,
    n_beta_time = length(sched$time),
    beta = NA_real_,
    sigma = 1 / latent_period,
    gamma = 1 / infectious_period,
    n_contact_time = length(mixing$time),
    contact_time = mixing$time,
    contact_value = mixing$value,
    sus = sus,
    trans = trans,
    p_hosp = p_hosp_s,
    p_death_comm = p_death_comm_s,
    p_death_hosp = p_death_hosp_s,
    waning_rate = gs(1, waning_rate),
    dur_hosp = as.vector(dur_hosp),
    waning_rate_R = as.vector(waning_rate_R),
    vacc_from = as.numeric(seq_len(n_imm) < vacc_to),
    vacc_to = one_hot(vacc_to),
    R_to = one_hot(R_to),
    S_ini = S_ini,
    E_ini = seed(E0),
    I_ini = seed(I0),
    R_ini = R_ini,
    beta_time = sched$time,
    beta_mult_value = sched$value,
    n_beta_t = length(bt$time),
    beta_t_time = bt$time,
    beta_t_value = bt$value,
    n_vacc_time = length(vacc$time),
    vacc_time = vacc$time,
    vacc_doses_value = vacc$value,
    # Observation model, only used when fitting (the likelihood replaces these)
    n_obs = 1L,
    obs_map = matrix(1, 1, n_group),
    nb_size = 10,
    obs_weekly = 0)

  if (is.null(beta)) {
    beta <- R0 / r0_per_beta(contact, pars$sus, pars$trans, N_group,
                             infectious_period)
  }
  check_range(beta, "beta")
  pars$beta <- beta

  if (sum(I0 + E0) == 0) {
    warning("No initial infections (I0 and E0 are all zero)")
  }

  structure(
    list(pars = pars,
         beta = beta,
         dims = dims,
         N = N,
         N_group = N_group,
         C_age = C_age,
         M_region = M_region,
         mixing_time = mixing$time,
         beta_t_end = bt$end,
         epi = list(susceptibility = susceptibility,
                    transmissibility = transmissibility,
                    p_hosp = p_hosp, p_death_hosp = p_death_hosp,
                    p_death_comm = p_death_comm, initial_R = initial_R),
         strata = list(rel_susceptibility = rel_susceptibility,
                       rel_transmissibility = rel_transmissibility,
                       rel_hosp = rel_hosp, rel_death = rel_death,
                       waning_rate = waning_rate,
                       initial_strata = initial_strata,
                       R_to_stratum = dims$stratum[R_to],
                       vacc_stratum = dims$stratum[vacc_to]),
         args = args),
    class = "seirstrata_model")
}

#' Change parameters of an existing model
#'
#' Rebuilds the model with some arguments of [build_model()] changed. Giving
#' `beta` replaces a previously given `R0` and vice versa.
#'
#' @param model A model from [build_model()].
#' @param ... Arguments of [build_model()] to change.
#' @return A new `seirstrata_model`.
#' @export
#' @examples
#' m <- build_model(N = 1e5, R0 = 2, I0 = 10)
#' m2 <- update_model(m, R0 = 1.5)
update_model <- function(model, ...) {
  changes <- list(...)
  if (length(changes) > 0 && (is.null(names(changes)) || any(names(changes) == ""))) {
    stop("All arguments to update_model() must be named")
  }
  unknown <- setdiff(names(changes), names(formals(build_model)))
  if (length(unknown) > 0) {
    stop(sprintf("Unknown argument(s): %s", paste(unknown, collapse = ", ")))
  }
  args <- model$args
  if ("beta" %in% names(changes)) args$R0 <- NULL
  if ("R0" %in% names(changes)) args$beta <- NULL
  args[names(changes)] <- changes
  do.call("build_model", args)
}

#' @export
print.seirstrata_model <- function(x, ...) {
  cat("<seirstrata_model>\n")
  cat(sprintf("  %d age group(s) x %d region(s) x %d immunity stratum/strata; population %s\n",
              length(x$dims$age), length(x$dims$region),
              length(x$dims$stratum),
              format(sum(x$N), big.mark = ",", scientific = FALSE)))
  cat(sprintf("  Age groups: %s\n", paste(x$dims$age, collapse = ", ")))
  cat(sprintf("  Regions:    %s\n", paste(x$dims$region, collapse = ", ")))
  cat(sprintf("  Strata:     %s\n", paste(x$dims$stratum, collapse = ", ")))
  cat(sprintf("  R0 = %.3g, beta = %.4g, R_eff(0) = %.3g\n",
              R0(x), x$beta, R_eff(x)))
  cat(sprintf("  Latent period %.3g days, infectious period %.3g days\n",
              1 / x$pars$sigma, 1 / x$pars$gamma))
  cat(sprintf("  Initial infected: %s (I0 + E0); in R: %s\n",
              format(sum(x$pars$I_ini + x$pars$E_ini), big.mark = ","),
              format(sum(x$pars$R_ini), big.mark = ",")))
  if (any(x$pars$waning_rate_R > 0) || any(x$pars$waning_rate > 0)) {
    cat(sprintf("  Waning: R -> stratum '%s'; ladder rates by stratum: %s\n",
                x$strata$R_to_stratum,
                paste(signif(colMeans(x$strata$waning_rate), 3),
                      collapse = ", ")))
  }
  if (!is.null(x$args$beta_t)) {
    cat(sprintf("  beta_t: %d value(s) from %.3g to %.3g%s\n", x$pars$n_beta_t - 1L,
                min(x$pars$beta_t_value), max(x$pars$beta_t_value),
                if (is.na(x$beta_t_end)) "" else
                  sprintf(", daily until day %g", x$beta_t_end)))
  }
  if (x$pars$n_beta_time > 1) {
    cat(sprintf("  Contact multipliers change at time(s): %s\n",
                paste(x$pars$beta_time[-1], collapse = ", ")))
  }
  if (any(x$pars$vacc_doses_value > 0)) {
    doses <- colSums(x$pars$vacc_doses_value)
    cat(sprintf("  Vaccination (to stratum '%s'): %s\n", x$strata$vacc_stratum,
                paste(sprintf("%s doses/day from day %g",
                              format(round(doses), big.mark = ",",
                                     scientific = FALSE, trim = TRUE),
                              x$pars$vacc_time), collapse = ", ")))
  }
  cat("  Parameters by age (averaged over regions, population-weighted):\n")
  w <- x$N / rowSums(x$N)
  tab <- vapply(x$epi, function(p) rowSums(p * w), numeric(length(x$dims$age)))
  tab <- as.data.frame(matrix(tab, nrow = length(x$dims$age),
                              dimnames = list(x$dims$age, names(x$epi))))
  print(noquote(apply(t(tab), c(1, 2), function(v) format(signif(v, 3)))))
  if (length(x$dims$stratum) > 1) {
    cat("  Stratum multipliers (mean over age groups):\n")
    st <- x$strata[c("rel_susceptibility", "rel_transmissibility", "rel_hosp",
                     "rel_death", "waning_rate", "initial_strata")]
    st <- t(vapply(st, colMeans, numeric(length(x$dims$stratum))))
    colnames(st) <- x$dims$stratum
    print(noquote(apply(st, c(1, 2), function(v) format(signif(v, 3)))))
  }
  invisible(x)
}
