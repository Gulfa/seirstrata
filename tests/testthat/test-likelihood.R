# Weekly observation model, data preparation, model_pars and the likelihood.

# Noise-free weekly admissions from the deterministic model
weekly_admissions <- function(model, weeks, by = character(0), pars = NULL) {
  d <- get_output(run_model(model, t_max = 7 * weeks, pars = pars),
                  "inc_hosp_week", by = by)
  d <- d[time > 0 & time %% 7 == 0]
  d[, admissions := round(value)]
  as.data.frame(d[, c("time", by, "admissions"), with = FALSE])
}

test_that("weekly admissions are the sum of the daily ones over each week", {
  m <- small_model()
  for (res in list(run_model(m, t_max = 70),
                   run_model(m, "stochastic", t_max = 70, n_sim = 3, seed = 1))) {
    d <- get_output(res, c("inc_hosp", "inc_hosp_week"))
    daily <- d[variable == "inc_hosp" & time > 0]
    daily[, week := ceiling(time / 7)]
    by_week <- daily[, .(value = sum(value)), by = .(sim, week)]
    weekly <- d[variable == "inc_hosp_week" & time > 0 & time %% 7 == 0]
    expect_equal(weekly$value, by_week$value, tolerance = 1e-6)
  }
})

test_that("model_data checks times, labels and builds the series", {
  m <- small_model()
  expect_error(model_data(data.frame(time = c(7, 10), admissions = 1), m),
               "7, 14")
  expect_silent(model_data(data.frame(time = c(1, 2), admissions = 1), m,
                           period = "day"))
  expect_error(model_data(data.frame(time = 7, admissions = 1.5), m),
               "whole numbers")
  expect_error(model_data(data.frame(time = 7, admissions = 1, age = "mid"), m),
               "Unknown age")
  expect_error(model_data(data.frame(time = c(7, 7), admissions = 1), m),
               "more than one row")

  d <- model_data(data.frame(time = c(7, 7, 14, 14), age = c("young", "old"),
                             admissions = c(1, 2, NA, 4)), m)
  expect_equal(d$series$label, c("young", "old"))
  # groups are age-fastest: young A, old A, young B, old B
  expect_equal(d$obs_map, rbind(c(1, 0, 1, 0), c(0, 1, 0, 1)))
  expect_equal(d$data_dust$admissions[[2]], c(NA, 4))
})

test_that("model age groups can map onto coarser data age groups", {
  N <- matrix(c(1e4, 2e4, 3e4), 3, dimnames = list(c("0-4", "5-14", "15+"), "all"))
  m <- build_model(N = N, C_age = diag(3) + 1, R0 = 2, I0 = c(5, 0, 0),
                   p_hosp = 0.05)
  map <- c("0-4" = "0-14", "5-14" = "0-14", "15+" = "15+")
  d <- model_data(data.frame(time = 7, age = c("0-14", "15+"),
                             admissions = c(3, 4)), m, age_groups = map)
  expect_equal(d$obs_map, rbind(c(1, 1, 0), c(0, 0, 1)))
  expect_error(model_data(data.frame(time = 7, admissions = 1), m,
                          age_groups = c("0-4" = "x")), "every model age group")
})

test_that("model_pars replaces parameters and keeps their sizes", {
  m <- small_model()
  p <- model_pars(m, list(beta = 0.1, p_hosp = m$pars$p_hosp * 2))
  expect_equal(p$beta, 0.1)
  expect_equal(p$p_hosp, m$pars$p_hosp * 2)
  expect_error(model_pars(m, list(nope = 1)), "Unknown model parameter")
  expect_error(model_pars(m, list(p_hosp = 1)), "keep its size")
  expect_error(model_pars(m, list(beta_t_value = c(1, 1, 1))), "keep its size")
  # run_model(pars = ) is the same as rebuilding with that beta
  a <- get_output(run_model(m, t_max = 60, pars = list(beta = m$beta * 0.8)), "H")
  b <- get_output(run_model(update_model(m, beta = m$beta * 0.8), t_max = 60), "H")
  expect_equal(a$value, b$value)
})

test_that("the likelihood is highest at the parameters that made the data", {
  m <- small_model(p_hosp = c(0.05, 0.2))
  data <- weekly_admissions(m, 15)
  lik <- model_likelihood(m, data, pars = "beta", nb_size = 50)
  f <- function(b) monty::monty_model_density(lik$likelihood, b)
  opt <- stats::optimize(f, c(0.5, 1.5) * m$beta, maximum = TRUE)
  expect_equal(opt$maximum, m$beta, tolerance = 0.01)
})

test_that("process maps fitted values onto model parameters", {
  m <- small_model(p_hosp = c(0.05, 0.2))
  truth <- list(beta = m$beta * 1.1, p_hosp = m$pars$p_hosp * 0.7)
  data <- weekly_admissions(m, 15, by = "age", pars = truth)
  process <- function(x) {
    list(beta = m$beta * x$beta_scale, p_hosp = m$pars$p_hosp * x$severity)
  }
  lik <- model_likelihood(m, data, pars = c("beta_scale", "severity", "nb_size"),
                          process = process,
                          initial = c(beta_scale = 1, severity = 1, nb_size = 20))
  expect_equal(nrow(lik$data$series), 2)
  f <- function(x) monty::monty_model_density(lik$likelihood, x)
  opt <- stats::optim(c(1, 1, 20), function(x) -f(x), method = "L-BFGS-B",
                      lower = c(0.5, 0.2, 1), upper = c(2, 3, 1e4))
  expect_equal(opt$par[1:2], c(1.1, 0.7), tolerance = 0.01)

  # the model parameters at the estimate reproduce the data
  p <- likelihood_pars(lik, c(beta_scale = 1.1, severity = 0.7, nb_size = 20))
  expect_setequal(names(p), c("beta", "p_hosp"))
  again <- weekly_admissions(m, 15, by = "age", pars = p)
  expect_equal(again$admissions, data$admissions)
})

test_that("seasonality can be fitted through beta_t", {
  days <- 0:140
  shape <- function(amp) 1 + amp * cos(2 * pi * (days - 120) / 365)
  m <- small_model(beta_t = shape(0))
  data <- weekly_admissions(m, 20, pars = list(beta_t_value = c(shape(0.3), shape(0.3)[141])))
  process <- function(x) {
    b <- shape(x$amp)
    list(beta_t_value = c(b, b[length(b)]))
  }
  lik <- model_likelihood(m, data, pars = "amp", process = process,
                          initial = c(amp = 0), nb_size = 50)
  f <- function(a) monty::monty_model_density(lik$likelihood, a)
  expect_equal(stats::optimize(f, c(0, 0.6), maximum = TRUE)$maximum, 0.3,
               tolerance = 0.02)
})

test_that("model_likelihood validates its inputs", {
  m <- small_model()
  data <- data.frame(time = c(7, 14), admissions = c(1, 2))
  expect_error(model_likelihood(m, data, pars = "scale"), "must be model parameters")
  expect_error(model_likelihood(m, data, pars = "p_hosp"), "scalars")
  expect_error(model_likelihood(m, data, pars = "s", process = function(x) list()),
               "'initial'")
  expect_error(model_likelihood(m, data, pars = "s", initial = c(s = 1),
                                process = function(x) list(nope = 1)),
               "Unknown model parameter")
  expect_error(model_likelihood(m, data, pars = "beta", initial = c(beta = 1),
                                process = function(x) list(beta = 1)),
               "must not return the fitted")
})

test_that("missing data and the stochastic filter work", {
  m <- small_model(p_hosp = c(0.05, 0.2))
  data <- weekly_admissions(m, 10)
  data$admissions[3] <- NA
  det <- model_likelihood(m, data, pars = "beta")
  expect_true(is.finite(monty::monty_model_density(det$likelihood, m$beta)))
  sto <- model_likelihood(m, data, pars = "beta", type = "stochastic",
                          n_particles = 50)
  ll <- monty::monty_model_density(sto$likelihood, m$beta)
  expect_true(is.finite(ll))
  expect_output(print(sto), "stochastic, fitting: beta")
})
