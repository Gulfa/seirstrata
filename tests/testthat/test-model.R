test_that("population is conserved", {
  m <- vax_model(coverage = 0.4, ve_infection = 0.5)
  det <- run_model(m, t_max = 200)
  expect_equal(total_alive_and_dead(det), rep(sum(m$N), 201), tolerance = 1e-6)
  sto <- run_model(m, "stochastic", t_max = 200, n_sim = 5, seed = 1)
  expect_identical(total_alive_and_dead(sto), rep(sum(m$N), 201 * 5))
})

test_that("incidence sums to the cumulative counts", {
  m <- small_model()
  for (res in list(run_model(m, t_max = 150),
                   run_model(m, "stochastic", t_max = 150, n_sim = 3, seed = 1))) {
    d <- get_output(res, c("inc_infections", "cum_infections", "inc_hosp", "cum_hosp"))
    tot <- dcast(d, sim + time ~ variable)
    tot[, `:=`(ci = cumsum(inc_infections), ch = cumsum(inc_hosp)), by = sim]
    expect_equal(tot$ci, tot$cum_infections, tolerance = 1e-6)
    expect_equal(tot$ch, tot$cum_hosp, tolerance = 1e-6)
  }
})

test_that("with all VE = 0 a second stratum changes nothing", {
  base <- run_model(small_model(), t_max = 200)
  vacc <- run_model(vax_model(coverage = 0.6), t_max = 200)
  a <- get_output(base, c("cum_infections", "cum_hosp", "D", "H"), by = c("age", "region"))
  b <- get_output(vacc, c("cum_infections", "cum_hosp", "D", "H"), by = c("age", "region"))
  expect_equal(a$value, b$value, tolerance = 1e-6)
})

test_that("a vaccine with ve_infection = 1 prevents all vaccinated infections", {
  m <- vax_model(coverage = 0.5, ve_infection = 1)
  for (res in list(run_model(m, t_max = 200),
                   run_model(m, "stochastic", t_max = 200, n_sim = 5, seed = 2))) {
    d <- get_output(res, "cum_infections", by = "stratum")
    expect_equal(max(d[stratum == "vaccinated"]$value), 0)
    expect_gt(max(d[stratum == "unvaccinated"]$value), 0)
  }
})

test_that("ve_hosp reduces hospitalisations per vaccinated infection", {
  m <- vax_model(coverage = 0.5, ve_hosp = 0.8)
  res <- run_model(m, t_max = 300)
  d <- get_output(res, c("cum_infections", "cum_hosp"), by = c("age", "stratum"))
  d <- dcast(d[time == 300], age + stratum ~ variable)
  d[, ratio := cum_hosp / cum_infections]
  # (not compared with the unvaccinated ratio: the seeded I0 are hospitalised
  # without being counted in cum_infections)
  expect_equal(d[stratum == "vaccinated"]$ratio, c(0.01, 0.1) * 0.2,
               tolerance = 1e-6)
})

test_that("isolated regions do not infect each other", {
  m <- small_model(M_region = diag(2))
  for (res in list(run_model(m, t_max = 150),
                   run_model(m, "stochastic", t_max = 150, n_sim = 5, seed = 3))) {
    d <- get_output(res, "cum_infections", by = "region")
    expect_equal(max(d[region == "B"]$value), 0)
  }
})

test_that("identical regions reproduce the single-region model", {
  N_age <- c(young = 20000, old = 10000)
  C_age <- matrix(c(3, 1, 1, 2), 2)
  one <- build_model(N = N_age, C_age = C_age, R0 = 2, I0 = c(5, 0), p_hosp = 0.1)
  two <- build_model(N = cbind(A = N_age, B = N_age), C_age = C_age,
                     M_region = matrix(0.5, 2, 2), R0 = 2,
                     I0 = matrix(c(5, 0), 2, 2), p_hosp = 0.1)
  expect_equal(one$beta, two$beta)
  a <- get_output(run_model(one, t_max = 100), "H", by = "age")
  b <- get_output(run_model(two, t_max = 100), "H", by = c("age", "region"))
  expect_equal(b[region == "A"]$value, a$value, tolerance = 1e-6)
})

test_that("single group final size matches the final size equation", {
  R0 <- 1.8
  m <- build_model(N = 1e6, R0 = R0, I0 = 10)
  res <- run_model(m, t_max = 600)
  z <- tail(get_output(res, "cum_infections")$value, 1) / 1e6
  z_exact <- uniroot(function(z) z - (1 - exp(-R0 * z)), c(0.1, 1))$root
  expect_equal(z, z_exact, tolerance = 1e-3)
})

test_that("R0 calibration and R_eff", {
  m <- small_model(R0 = 2.5)
  expect_equal(R0(m), 2.5)
  expect_equal(R_eff(m), 2.5, tolerance = 1e-3) # only 10 infected
  v <- do.call(update_model, c(list(m), vax_args(coverage = 0.5, ve_infection = 1)))
  expect_equal(R_eff(v), 1.25, tolerance = 1e-3)
  # single group: R0 = beta * sus * C * D
  s <- build_model(N = 1000, R0 = 2, I0 = 1, susceptibility = 0.5,
                   infectious_period = 4)
  expect_equal(s$beta, 2 / (0.5 * 4))
})

test_that("stochastic mean is close to the deterministic model", {
  m <- build_model(N = 2e5, R0 = 2, I0 = 50, p_hosp = 0.05, p_death_hosp = 0.1)
  det <- summary(run_model(m, t_max = 300))
  sto <- summary(run_model(m, "stochastic", t_max = 300, n_sim = 20,
                           dt = 0.05, seed = 1))
  for (k in c("cum_infections", "cum_hosp")) {
    expect_equal(sto[measure == k]$mean, det[measure == k]$value,
                 tolerance = 0.02)
  }
})

test_that("stochastic runs are reproducible with a seed", {
  m <- small_model()
  a <- get_output(run_model(m, "stochastic", t_max = 50, n_sim = 3, seed = 42), "I")
  b <- get_output(run_model(m, "stochastic", t_max = 50, n_sim = 3, seed = 42), "I")
  expect_identical(a, b)
})

test_that("continuing from a state reproduces an uninterrupted run", {
  m <- small_model()
  full <- run_model(m, t_max = 120)
  # the state halfway through, as a filter would hand it over
  sys <- dust2::dust_system_create(seirstrata:::seir_ode(), m$pars)
  dust2::dust_system_set_state_initial(sys)
  y <- dust2::dust_system_simulate(sys, 0:60)
  state <- y[, ncol(y)]

  cont <- run_model(m, t_max = 120, state = state, t_start = 60)
  for (v in c("cum_hosp", "H", "inc_hosp", "inc_infections")) {
    expect_equal(get_output(cont, v)$value,
                 get_output(full, v)[time >= 60]$value,
                 info = v)
  }
})
