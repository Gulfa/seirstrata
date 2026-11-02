# R0 and R_eff against the growth rate of simulated epidemics.
#
# For SEIR with exponentially distributed stages, the early exponential growth
# rate r and the reproduction number satisfy R = (1 + r / sigma)(1 + r / gamma).
# In a large population, susceptible depletion is negligible early on, so the
# growth rate measured from simulated incidence gives R independently of the
# next-generation matrix.

growth_R <- function(model, from = 60, to = 90) {
  inc <- get_output(run_model(model, t_max = to), "inc_infections")
  inc <- inc[time >= from & time <= to]
  r <- unname(stats::coef(stats::lm(log(value) ~ time, data = inc))[2])
  p <- model$pars
  (1 + r / p$sigma) * (1 + r / p$gamma)
}

test_that("R0 matches the growth rate of a single-group epidemic", {
  m <- build_model(N = 1e9, R0 = 1.5, I0 = 10)
  expect_equal(growth_R(m), R0(m), tolerance = 0.005)
})

test_that("R_eff matches the growth rate with age, regions and strata", {
  N <- matrix(c(2e8, 3e8, 1e8, 1.5e8), 2,
              dimnames = list(c("young", "old"), c("A", "B")))
  m <- build_model(N = N, C_age = matrix(c(3, 1, 1, 2), 2),
                   M_region = matrix(c(0.9, 0.1, 0.1, 0.9), 2), R0 = 3,
                   I0 = matrix(c(10, 0, 0, 0), 2),
                   susceptibility = c(0.8, 1),
                   strata = c("none", "partial", "good"),
                   rel_susceptibility = matrix(c(1, 1, 0.6, 0.4, 0.3, 0.1), 2),
                   rel_transmissibility = c(1, 0.7, 0.5),
                   initial_strata = matrix(c(0.5, 0.3, 0.3, 0.3, 0.2, 0.4), 2),
                   initial_R = c(0.2, 0.1),
                   beta_t = data.frame(time = 0, mult = 0.9))
  R <- R_eff(m)
  expect_gt(R, 1.05)    # still growing, so the test means something
  expect_lt(R, R0(m))
  expect_equal(growth_R(m), R, tolerance = 0.01)
  # and a declining epidemic
  low <- update_model(m, R0 = 2)
  expect_lt(R_eff(low), 1)
  expect_equal(growth_R(low), R_eff(low), tolerance = 0.01)
})

test_that("R_eff follows beta_t over time", {
  m <- build_model(N = 1e9, R0 = 1.5, I0 = 10,
                   beta_t = data.frame(time = c(0, 59, 60), mult = c(1, 1, 0.8)))
  expect_equal(growth_R(m, 80, 110), R_eff(m, 80), tolerance = 0.005)
  expect_equal(R_eff(m, 80), 1.5 * 0.8)
})
