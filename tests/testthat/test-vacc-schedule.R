campaign <- function(time = 20, doses = 500, ...) {
  data.frame(time = c(time, time + 30), doses = c(doses, 0), ...)
}

test_that("population is conserved with vaccination over time", {
  m <- vax_model(vacc_schedule = campaign(), ve_infection = 0.5)
  det <- run_model(m, t_max = 200)
  expect_equal(total_alive_and_dead(det), rep(sum(m$N), 201), tolerance = 1e-6)
  sto <- run_model(m, "stochastic", t_max = 200, n_sim = 5, seed = 1)
  expect_identical(total_alive_and_dead(sto), rep(sum(m$N), 201 * 5))
})

test_that("the scheduled doses are given, and only while the schedule runs", {
  m <- vax_model(vacc_schedule = campaign(time = 20, doses = 500))
  d <- get_output(run_model(m, t_max = 100), c("cum_doses", "inc_doses"))
  doses <- d[variable == "cum_doses"]$value
  # (the ODE solver can step a hair across the start: allow a rounding error)
  expect_equal(doses[d[variable == "cum_doses"]$time <= 20], rep(0, 21),
               tolerance = 1e-3)
  expect_equal(max(doses), 500 * 30, tolerance = 1e-6)
  inc <- d[variable == "inc_doses"]
  # (the solver can step a hair across the end too)
  expect_equal(inc[time %in% 21:50]$value, rep(500, 30), tolerance = 1e-4)
  expect_equal(sum(inc[time > 50]$value), 0)

  # doses have no stratum
  by_stratum <- get_output(run_model(m, t_max = 100), "cum_doses", by = "stratum")
  expect_equal(unique(as.character(by_stratum$stratum)), "all")
})

test_that("the epidemic before the campaign is unchanged", {
  base <- run_model(vax_model(ve_infection = 0.7), t_max = 100)
  vacc <- run_model(vax_model(ve_infection = 0.7,
                                vacc_schedule = campaign(time = 40)), t_max = 100)
  a <- get_output(base, "cum_infections")
  b <- get_output(vacc, "cum_infections")
  expect_equal(a[time <= 40]$value, b[time <= 40]$value, tolerance = 1e-6)
  expect_lt(max(b$value), max(a$value))
})

test_that("with all VE = 0 vaccination over time changes nothing", {
  base <- run_model(vax_model(), t_max = 200)
  vacc <- run_model(vax_model(vacc_schedule = campaign()), t_max = 200)
  a <- get_output(base, c("cum_infections", "cum_hosp", "D"), by = c("age", "region"))
  b <- get_output(vacc, c("cum_infections", "cum_hosp", "D"), by = c("age", "region"))
  expect_equal(a$value, b$value, tolerance = 1e-6)
})

test_that("targeted doses go to the targeted groups, shared by population", {
  m <- vax_model(vacc_schedule = campaign(doses = 200, age = "old"))
  d <- get_output(run_model(m, t_max = 60), "cum_doses",
                  by = c("age", "region"))[time == 60]
  expect_equal(sum(d[age == "young"]$value), 0)
  old <- d[age == "old"]
  # old: 30000 in A, 15000 in B
  expect_equal(old$value, 200 * 30 * c(2, 1) / 3, tolerance = 1e-6)

  m <- vax_model(vacc_schedule = campaign(doses = 200, region = "B"))
  d <- get_output(run_model(m, t_max = 60), "cum_doses", by = "region")[time == 60]
  expect_equal(d[region == "A"]$value, 0)
  expect_equal(d[region == "B"]$value, 200 * 30, tolerance = 1e-6)
})

test_that("vaccination of S stops when the strata below run out", {
  # far more doses than there are young people in region A. Vaccination
  # status is not tracked, so people in R stay eligible and take the doses
  # left over: cum_doses is not capped by the population.
  m <- vax_model(vacc_schedule = data.frame(time = 0, doses = 1e5,
                                            age = "young", region = "A"))
  for (res in list(run_model(m, t_max = 30),
                   run_model(m, "stochastic", t_max = 30, n_sim = 3, seed = 1))) {
    d <- get_output(res, c("S", "R"), by = c("age", "region", "stratum"))
    expect_gte(min(d$value), -1e-6)
    S <- d[variable == "S" & age == "young" & region == "A"]
    expect_lt(min(S[stratum == "unvaccinated"]$value), 1)
    expect_lte(max(S[stratum == "vaccinated"]$value), 20000 + 1e-6)
    expect_gt(max(S[stratum == "vaccinated"]$value), 0.99 * (20000 - 10))
  }
})

test_that("stochastic and deterministic vaccination agree", {
  m <- vax_model(vacc_schedule = campaign(doses = 300), ve_infection = 0.5)
  det <- get_output(run_model(m, t_max = 100), "cum_doses")[time == 100]$value
  sto <- get_output(run_model(m, "stochastic", t_max = 100, n_sim = 50,
                              seed = 1), "cum_doses")[time == 100]$value
  expect_equal(mean(sto), det, tolerance = 0.01)
})

test_that("a later vaccination schedule acts on a continued state", {
  # as when continuing from a filtered state: continue a run from day 30 with
  # a model that starts vaccinating at day 40
  m <- vax_model(ve_infection = 0.8)
  sys <- dust2::dust_system_create(seir_ode(), m$pars)
  dust2::dust_system_set_state_initial(sys)
  dust2::dust_system_run_to_time(sys, 30)
  state <- dust2::dust_system_state(sys)

  vacc <- update_model(m, vacc_schedule = campaign(time = 40, doses = 1000))
  a <- run_model(m, t_max = 150, state = state, t_start = 30)
  b <- run_model(vacc, t_max = 150, state = state, t_start = 30)
  ia <- get_output(a, "cum_infections")
  ib <- get_output(b, "cum_infections")
  expect_equal(ia[time <= 40]$value, ib[time <= 40]$value, tolerance = 1e-6)
  expect_lt(max(ib$value), max(ia$value))
  expect_equal(max(get_output(b, "cum_doses")$value), 30000, tolerance = 1e-6)
})

test_that("vacc_schedule is validated", {
  expect_error(vax_model(vacc_schedule = data.frame(time = 1, mult = 1)),
               "'time' and 'doses'")
  expect_error(vax_model(vacc_schedule = data.frame(time = 1, doses = -1)),
               "doses")
  expect_error(vax_model(vacc_schedule = data.frame(time = 1, doses = 1,
                                                      age = "middle")),
               "Unknown age")
  expect_error(vax_model(vacc_schedule = data.frame(time = 1, doses = 1,
                                                      when = 2)),
               "Unknown column")
})

test_that("print shows the vaccination schedule", {
  m <- vax_model(vacc_schedule = campaign(time = 20, doses = 1500))
  expect_output(print(m), "1,500 doses/day from day 20, 0 doses/day from day 50")
})
