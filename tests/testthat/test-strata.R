# The immunity ladder, R and waning.

ladder_model <- function(...) {
  args <- list(strata = c("none", "partial", "good"),
               rel_susceptibility = c(1, 0.6, 0.3),
               rel_transmissibility = c(1, 0.8, 0.6),
               rel_hosp = c(1, 0.5, 0.2),
               rel_death = c(1, 0.7, 0.4),
               waning_rate = c(0, 1 / 60, 1 / 90),
               waning_rate_R = 1 / 120,
               initial_strata = c(0.5, 0.3, 0.2),
               initial_R = 0.1)
  args[names(list(...))] <- list(...)
  do.call(small_model, args)
}

# No transmission, so only waning and vaccination move people
no_transmission <- function(...) {
  suppressWarnings(ladder_model(R0 = NULL, beta = 0, I0 = 0, ...))
}

test_that("2 strata without waning reproduce seirvax", {
  ref <- readRDS(test_path("fixtures", "seirvax_reference.rds"))
  # seirvax argument -> seirstrata argument
  translate <- function(a) {
    s <- vax_args(coverage = a$vacc_coverage %||% 0,
                  ve_infection = a$ve_infection %||% 0,
                  ve_transmission = a$ve_transmission %||% 0,
                  ve_hosp = a$ve_hosp %||% 0, ve_death = a$ve_death %||% 0)
    a[c("vacc_coverage", "ve_infection", "ve_transmission", "ve_hosp",
        "ve_death")] <- NULL
    c(a, s)
  }
  for (case in names(ref)) {
    r <- ref[[case]]
    m <- do.call(build_model, translate(r$args))
    expect_equal(m$beta, r$beta, info = case)
    expect_equal(R0(m), r$R0, info = case)
    expect_equal(R_eff(m), r$R_eff, info = case)
    out <- get_output(run_model(m, t_max = 200), unique(r$out$variable),
                      by = c("age", "region"))
    expected <- as.data.table(r$out)
    setorderv(expected, c("variable", "age", "region", "sim", "time"))
    # "coverage" agrees to ODE solver tolerance. In "doses", seirvax moved
    # vaccinated recovered people to a vaccinated R that took no more doses;
    # here R stays eligible (vaccination status is not tracked), so slightly
    # more doses are wasted on R: about 0.1% more infections.
    tol <- if (case == "doses") 1e-2 else 1e-4
    expect_equal(out$value, expected$value, tolerance = tol, info = case)
  }
})

test_that("population is conserved with waning, R and vaccination", {
  m <- ladder_model(vacc_stratum = "partial",
                    vacc_schedule = data.frame(time = 10, doses = 400))
  det <- run_model(m, t_max = 300)
  expect_equal(total_alive_and_dead(det), rep(sum(m$N), 301), tolerance = 1e-6)
  sto <- run_model(m, "stochastic", t_max = 300, n_sim = 5, seed = 1)
  expect_identical(total_alive_and_dead(sto), rep(sum(m$N), 301 * 5))
  # and nothing goes negative
  d <- get_output(sto, c("S", "E", "I", "H", "R"), by = c("age", "region", "stratum"))
  expect_gte(min(d$value), 0)
})

test_that("waning moves people down the ladder at the given rates", {
  w <- 1 / 50
  m <- no_transmission(initial_strata = c(0, 0, 1), initial_R = 0,
                       waning_rate = c(0, w, w))
  d <- get_output(run_model(m, t_max = 100), "S", by = "stratum")
  N <- sum(m$N)
  t <- d[stratum == "good"]$time
  # a chain of two exponential stages
  expect_equal(d[stratum == "good"]$value, N * exp(-w * t), tolerance = 1e-6)
  expect_equal(d[stratum == "partial"]$value, N * w * t * exp(-w * t),
               tolerance = 1e-6)
})

test_that("R wanes into R_to_stratum", {
  r <- 1 / 40
  m <- no_transmission(initial_strata = c(1, 0, 0), initial_R = 0.5,
                       waning_rate = 0, waning_rate_R = r,
                       R_to_stratum = "partial")
  res <- run_model(m, t_max = 100)
  N <- sum(m$N)
  R <- get_output(res, "R")
  expect_equal(R$value, 0.5 * N * exp(-r * R$time), tolerance = 1e-6)
  S <- get_output(res, "S", by = "stratum")
  expect_equal(S[stratum == "partial"]$value, 0.5 * N * (1 - exp(-r * S[stratum == "partial"]$time)),
               tolerance = 1e-6)
  expect_equal(max(S[stratum == "good"]$value), 0)
})

test_that("with waning_rate_R = 0 immunity from infection lasts", {
  m <- ladder_model(waning_rate = 0, waning_rate_R = 0)
  d <- get_output(run_model(m, t_max = 400), "R")
  expect_true(all(diff(d$value) >= -1e-8))
})

test_that("vaccination moves people from the strata below vacc_stratum only", {
  m <- no_transmission(initial_strata = c(0.5, 0, 0.5), initial_R = 0,
                       waning_rate = 0, vacc_stratum = "partial",
                       vacc_schedule = data.frame(time = 0, doses = 100))
  d <- get_output(run_model(m, t_max = 50), c("S", "cum_doses"), by = "stratum")
  S <- d[variable == "S"]
  expect_equal(S[stratum == "good"]$value, rep(sum(m$N) / 2, 51), tolerance = 1e-6)
  expect_equal(tail(S[stratum == "partial"]$value, 1), 100 * 50, tolerance = 1e-6)
  expect_equal(tail(d[variable == "cum_doses"]$value, 1), 100 * 50, tolerance = 1e-6)
})

test_that("doses given to R are counted but change nothing", {
  m <- no_transmission(initial_strata = c(0, 0, 1), initial_R = 0.5,
                       waning_rate = 0, waning_rate_R = 0,
                       vacc_schedule = data.frame(time = 0, doses = 100))
  d <- get_output(run_model(m, t_max = 20), c("R", "cum_doses"))
  expect_equal(unique(round(d[variable == "R"]$value)), sum(m$N) / 2)
  expect_equal(tail(d[variable == "cum_doses"]$value, 1), 100 * 20, tolerance = 1e-6)
})

test_that("infections are recorded in the stratum of the person infected", {
  m <- ladder_model(waning_rate = 0, waning_rate_R = 0, initial_R = 0)
  d <- get_output(run_model(m, t_max = 1000), c("cum_infections", "cum_hosp"),
                  by = c("age", "stratum"))
  d <- dcast(d[time == 1000], age + stratum ~ variable)
  d[, ratio := cum_hosp / cum_infections]
  # protected strata: p_hosp * rel_hosp (stratum 1 also has the seeds)
  expect_equal(d[stratum == "partial"]$ratio, c(0.01, 0.1) * 0.5, tolerance = 1e-6)
  expect_equal(d[stratum == "good"]$ratio, c(0.01, 0.1) * 0.2, tolerance = 1e-6)
})

test_that("R_eff accounts for the initial strata and R", {
  m <- build_model(N = 1e5, R0 = 2, I0 = 1, strata = c("a", "b"),
                   rel_susceptibility = c(1, 0.5),
                   initial_strata = c(0.6, 0.4), initial_R = 0.2)
  expect_equal(R0(m), 2)
  # single group: R_eff = R0 * sum_k S_k * rel_sus_k / N
  expect_equal(R_eff(m), 2 * 0.8 * (0.6 + 0.4 * 0.5), tolerance = 1e-4)
})

test_that("stochastic mean agrees with the deterministic model with strata", {
  m <- build_model(N = 2e5, R0 = 2.5, I0 = 50, p_hosp = 0.05,
                   strata = c("none", "partial", "good"),
                   rel_susceptibility = c(1, 0.6, 0.3), rel_hosp = c(1, 0.5, 0.2),
                   waning_rate = 1 / 60, waning_rate_R = 1 / 90,
                   initial_strata = c(0.5, 0.3, 0.2), initial_R = 0.1)
  det <- summary(run_model(m, t_max = 300))
  sto <- summary(run_model(m, "stochastic", t_max = 300, n_sim = 20,
                           dt = 0.05, seed = 1))
  for (k in c("cum_infections", "cum_hosp")) {
    expect_equal(sto[measure == k]$mean, det[measure == k]$value,
                 tolerance = 0.03, info = k)
  }
})

test_that("strata arguments are validated", {
  expect_error(ladder_model(waning_rate = c(0.1, 0.1, 0.1)), "stratum 1")
  expect_error(ladder_model(initial_strata = c(0.5, 0.5, 0.5)), "sum to 1")
  expect_error(ladder_model(rel_hosp = c(1, 0.5)), "vector by stratum")
  expect_error(ladder_model(R_to_stratum = "best"), "R_to_stratum")
  expect_error(ladder_model(vacc_stratum = 4), "vacc_stratum")
  expect_error(small_model(strata = c("a", "a")), "unique")
  expect_error(small_model(strata = c("all", "b")), "reserved")
  expect_error(small_model(vacc_schedule = data.frame(time = 1, doses = 1)),
               "vacc_stratum")
  expect_error(ladder_model(rel_hosp = c(1, 20, 1)), "p_hosp \\* rel_hosp")
})

test_that("print shows the strata", {
  expect_output(print(ladder_model()), "Strata: +none, partial, good")
  expect_output(print(ladder_model()), "Stratum multipliers")
})
