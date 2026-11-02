# The global transmission multiplier beta_t.

cum_inf <- function(m, t_max = 150, by = character(0)) {
  get_output(run_model(m, t_max = t_max), "cum_infections", by = by)$value
}

test_that("a constant beta_t scales beta", {
  m <- small_model()
  scaled <- update_model(m, beta = m$beta * 0.8)
  expect_equal(cum_inf(update_model(m, beta_t = rep(0.8, 200))),
               cum_inf(scaled), tolerance = 1e-6)
  expect_equal(cum_inf(update_model(m, beta_t = data.frame(time = 0, mult = 0.8))),
               cum_inf(scaled), tolerance = 1e-6)
})

test_that("beta_t multiplies with beta_schedule", {
  m <- small_model()
  both <- update_model(m, beta_t = data.frame(time = 0, mult = 0.8),
                       beta_schedule = data.frame(time = 0, mult = 0.9))
  expect_equal(cum_inf(both), cum_inf(update_model(m, beta = m$beta * 0.72)),
               tolerance = 1e-6)
})

test_that("daily values are knots, linear in between", {
  m <- small_model()
  # transmission goes from full at day 29 to none at day 30: no new
  # infections after day 30
  stop30 <- update_model(m, beta_t = c(rep(1, 30), rep(0, 171)))
  inc <- get_output(run_model(stop30, t_max = 200), "inc_infections")
  expect_equal(max(inc[time > 30]$value), 0)
  expect_gt(max(inc[time <= 30]$value), 0)
  # identical up to day 29
  base <- get_output(run_model(m, t_max = 200), "inc_infections")
  expect_equal(inc[time <= 29]$value, base[time <= 29]$value, tolerance = 1e-6)
  # the same as a data.frame with the same knots
  df <- update_model(m, beta_t = data.frame(time = c(0, 29, 30),
                                            mult = c(1, 1, 0)))
  expect_equal(cum_inf(df, 200), cum_inf(stop30, 200), tolerance = 1e-6)
})

test_that("R0 excludes beta_t and R_eff includes it", {
  m <- small_model(R0 = 2, beta_t = c(rep(1.2, 10), rep(0.5, 10)))
  expect_equal(R0(m), 2)
  expect_equal(R_eff(m, 0), 2 * 1.2, tolerance = 1e-3)
  expect_equal(R_eff(m, 15), 2 * 0.5, tolerance = 1e-3)
  expect_equal(R_eff(m, 9.5), 2 * 0.85, tolerance = 1e-3)   # between knots
  expect_equal(R_eff(m, 100), 2 * 0.5, tolerance = 1e-3)    # last value held
})

test_that("the stochastic model uses beta_t", {
  m <- build_model(N = 2e5, R0 = 2, I0 = 50,
                   beta_t = 1 + 0.3 * cos(2 * pi * (0:299) / 365))
  det <- summary(run_model(m, t_max = 299))
  sto <- summary(run_model(m, "stochastic", t_max = 299, n_sim = 20, dt = 0.05,
                           seed = 1))
  k <- "cum_infections"
  expect_equal(sto[measure == k]$mean, det[measure == k]$value, tolerance = 0.02)
})

test_that("running beyond a daily beta_t warns", {
  m <- small_model(beta_t = rep(1, 50))
  expect_warning(run_model(m, t_max = 60), "covers days 0 to 49")
  expect_silent(run_model(m, t_max = 49))
  # a data.frame holds its last value on purpose
  expect_silent(run_model(small_model(beta_t = data.frame(time = c(0, 10),
                                                          mult = c(1, 0.5))),
                          t_max = 60))
})

test_that("beta_t is validated", {
  expect_error(small_model(beta_t = c(1, -0.1)), "non-negative")
  expect_error(small_model(beta_t = c(1, NA)), "non-negative")
  expect_error(small_model(beta_t = "a"), "numeric vector")
  expect_error(small_model(beta_t = data.frame(time = 1, value = 1)),
               "'time' and 'mult'")
  expect_error(small_model(beta_t = data.frame(time = c(5, 2), mult = 1)),
               "strictly increasing")
  expect_error(small_model(beta_t = data.frame(time = 5, mult = 1)),
               "start at 0")
})

test_that("print shows beta_t", {
  expect_output(print(small_model(beta_t = c(1, 0.8, 1.1))),
                "beta_t: 3 value\\(s\\) from 0.8 to 1.1, daily until day 2")
})
