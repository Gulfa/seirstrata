test_that("get_output aggregates consistently", {
  m <- vax_model(coverage = 0.3)
  res <- run_model(m, t_max = 50)
  full <- get_output(res, "I", by = c("age", "region", "stratum"))
  tot <- get_output(res, "I")
  expect_equal(full[, .(v = sum(value)), by = time]$v, tot$value)
  expect_equal(levels(get_output(res, "I", by = "region")$region), c("A", "B"))
  expect_error(get_output(res, "nope"), "Unknown variable")
  expect_error(get_output(res, "I", by = "county"), "subset")
})

test_that("summary and plot work for both model types", {
  m <- small_model()
  det <- run_model(m, t_max = 100)
  sto <- run_model(m, "stochastic", t_max = 100, n_sim = 4, seed = 1)
  s1 <- summary(det, by = "age")
  expect_setequal(names(s1), c("age", "measure", "value"))
  s2 <- summary(sto)
  expect_setequal(names(s2), c("measure", "mean", "median", "lower", "upper"))
  expect_s3_class(plot(det, "H", by = c("age", "region")), "ggplot")
  expect_s3_class(plot(sto, "inc_hosp", by = "age", n_traces = 2), "ggplot")
})

test_that("stochastic runs are spaghetti by default, ribbon on request", {
  m <- small_model()
  sto <- run_model(m, "stochastic", t_max = 60, n_sim = 5, seed = 1)
  layers <- function(p) vapply(p$layers, function(l) class(l$geom)[1], "")

  spaghetti <- plot(sto, "I")
  expect_true("GeomLine" %in% layers(spaghetti))
  expect_false("GeomRibbon" %in% layers(spaghetti))
  # one line per simulation, plus the median
  expect_equal(length(unique(spaghetti$layers[[1]]$data$.trace)), 5)
  expect_equal(sum(layers(spaghetti) == "GeomLine"), 2)
  expect_equal(sum(layers(plot(sto, "I", show_median = FALSE)) == "GeomLine"), 1)
  expect_equal(length(unique(plot(sto, "I", n_traces = 2)$layers[[1]]$data$.trace)), 2)

  ribbon <- plot(sto, "I", style = "ribbon")
  expect_true("GeomRibbon" %in% layers(ribbon))
  expect_equal(sum(layers(ribbon) == "GeomLine"), 1)      # no traces by default
  expect_equal(sum(layers(plot(sto, "I", style = "ribbon",
                               n_traces = 3)) == "GeomLine"), 2)

  # deterministic runs ignore the style
  det <- run_model(m, t_max = 60)
  expect_equal(unname(layers(plot(det, "I"))), "GeomLine")
  expect_error(plot(sto, "I", style = "wiggly"), "should be one of")
})

test_that("run_model can keep a subset of variables", {
  res <- run_model(small_model(), t_max = 20, vars = c("H", "inc_deaths"))
  expect_equal(names(res$states), c("H", "inc_deaths"))
  expect_error(run_model(small_model(), t_max = 20, vars = "nope"), "Unknown")
})

test_that("summary reports what the kept variables support", {
  m <- small_model()
  full <- summary(run_model(m, t_max = 80))
  expect_setequal(full$measure, c("cum_infections", "cum_hosp", "D",
                                  "peak_hosp_occupancy", "peak_hosp_day"))

  # a run that kept only some variables gives the measures it can
  peaks <- summary(run_model(m, t_max = 80, vars = c("inc_hosp", "H")))
  expect_setequal(peaks$measure, c("peak_hosp_occupancy", "peak_hosp_day"))
  expect_equal(peaks$value, full[measure %in% peaks$measure]$value)

  totals <- summary(run_model(m, t_max = 80, vars = c("cum_hosp", "D")))
  expect_setequal(totals$measure, c("cum_hosp", "D"))

  # and the same through a comparison
  cmp <- compare_scenarios(list(a = m, b = update_model(m, R0 = 1.2)),
                           t_max = 80, vars = c("inc_hosp", "H"))
  expect_setequal(summary(cmp)$measure,
                  c("peak_hosp_occupancy", "peak_hosp_day"))

  # nothing usable is a clear error, not a cryptic one
  expect_error(summary(run_model(m, t_max = 50, vars = "inc_hosp")),
               "needs at least one of")
})

test_that("plots use the Dark2 palette for groups", {
  m <- small_model()
  sto <- run_model(m, "stochastic", t_max = 60, n_sim = 5, seed = 1)
  dark2 <- scales::brewer_pal(palette = "Dark2")(8)

  scale_values <- function(p) {
    unlist(lapply(p$scales$scales, function(x) x$palette(2)))
  }
  cols <- scale_values(plot(sto, "H", by = "region"))
  expect_true(all(dark2[1:2] %in% cols))

  # more groups than the palette has: interpolated, not recycled or NA
  many <- dark2_colours(20)
  expect_length(many, 20)
  expect_equal(anyDuplicated(many), 0L)
  expect_equal(dark2_colours(3), dark2[1:3])

  cmp <- compare_scenarios(list(a = m, b = update_model(m, R0 = 1.2)),
                           t_max = 50)
  expect_true(all(dark2[1:2] %in% scale_values(plot(cmp, var = "H"))))
})

test_that("spaghetti lines are opaque enough to follow", {
  m <- small_model()
  sto <- run_model(m, "stochastic", t_max = 60, n_sim = 100, seed = 1)
  line_alpha <- function(p) p$layers[[1]]$aes_params$alpha

  expect_gte(line_alpha(plot(sto, "H", n_traces = 100)), 0.25)
  # fewer lines, more opaque
  expect_gt(line_alpha(plot(sto, "H", n_traces = 10)),
            line_alpha(plot(sto, "H", n_traces = 100)))
  # and never invisible, however many
  expect_gte(trace_alpha(1000), 0.3)
  expect_lte(trace_alpha(1), 0.85)
  # explicit alpha wins
  expect_equal(line_alpha(plot(sto, "H", alpha = 0.5)), 0.5)
})
