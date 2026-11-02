test_that("compare_scenarios accepts models and runs", {
  m <- small_model()
  m2 <- update_model(m, R0 = 1.2)

  # models are run for you
  cmp <- compare_scenarios(list(base = m, low = m2), t_max = 100)
  expect_s3_class(cmp, "seirstrata_comparison")
  expect_equal(cmp$scenarios, c("base", "low"))
  expect_equal(max(cmp$runs$base$time), 100)

  # already-run scenarios are taken as they are
  runs <- compare_scenarios(list(base = run_model(m, t_max = 50),
                                 low = run_model(m2, t_max = 50)))
  expect_equal(max(runs$runs$base$time), 50)

})

test_that("compare_scenarios validates its input", {
  m <- small_model()
  expect_error(compare_scenarios(list(a = m)), "at least two")
  expect_error(compare_scenarios(list(m, m)), "unique names")
  expect_error(compare_scenarios(list(a = m, a = m)), "unique names")
  expect_error(compare_scenarios(list(a = m, b = "x")), "must be a model or a run")
  expect_error(compare_scenarios(list(a = m, b = build_model(N = 1e4, R0 = 2,
                                                             I0 = 5))),
               "same age groups and regions")
})

test_that("get_output, summary and plot work on a comparison", {
  m <- small_model()
  cmp <- compare_scenarios(list(base = m, low = update_model(m, R0 = 1.2)),
                           type = "stochastic", t_max = 80, n_sim = 5, seed = 1)

  d <- get_output(cmp, "H", by = "age")
  expect_setequal(names(d), c("scenario", "time", "sim", "age", "variable",
                              "value"))
  expect_equal(levels(d$scenario), c("base", "low"))
  # the same numbers as asking each run on its own
  expect_equal(d[scenario == "base"]$value,
               get_output(cmp$runs$base, "H", by = "age")$value)

  s <- summary(cmp)
  expect_equal(levels(s$scenario), c("base", "low"))
  expect_setequal(names(s), c("scenario", "measure", "mean", "median",
                              "lower", "upper"))
  peak <- s[measure == "peak_hosp_occupancy"]
  expect_lt(peak[scenario == "low"]$median, peak[scenario == "base"]$median)
  expect_true("age" %in% names(summary(cmp, by = "age")))

  expect_s3_class(plot(cmp, var = "H"), "ggplot")
  expect_s3_class(plot(cmp, var = "H", style = "spaghetti"), "ggplot")
  expect_s3_class(plot(cmp, var = "inc_hosp", by = "region"), "ggplot")
  expect_output(print(cmp), "2 scenarios")
  expect_error(plot(cmp, var = c("H", "I")), "single variable")
})

test_that("deterministic and stochastic scenarios can be mixed", {
  m <- small_model()
  cmp <- compare_scenarios(list(det = run_model(m, t_max = 60),
                                sto = run_model(m, "stochastic", t_max = 60,
                                                n_sim = 4, seed = 1)))
  s <- summary(cmp)
  # the deterministic scenario has no interval, the stochastic one does
  expect_true(all(is.na(s[scenario == "det"]$lower)))
  expect_false(any(is.na(s[scenario == "sto"]$lower)))
  expect_s3_class(plot(cmp, var = "H"), "ggplot")
})
