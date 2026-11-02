dims <- list(age = c("a1", "a2", "a3"), region = c("r1", "r2"))

test_that("expand_param recycles scalars, age vectors, region vectors and matrices", {
  expect_equal(unname(expand_param(2, "x", dims)), matrix(2, 3, 2))
  expect_equal(unname(expand_param(1:3, "x", dims)), matrix(1:3, 3, 2))
  expect_equal(unname(expand_param(c(5, 6), "x", dims)),
               matrix(c(5, 6), 3, 2, byrow = TRUE))
  expect_equal(unname(expand_param(c(r2 = 6, r1 = 5), "x", dims)),
               matrix(c(5, 6), 3, 2, byrow = TRUE))
  m <- matrix(1:6, 3, 2)
  expect_equal(unname(expand_param(m, "x", dims)), m)
  expect_error(expand_param(1:4, "x", dims), "vector by age")
  expect_error(expand_param(matrix(1, 2, 3), "x", dims), "3 x 2")
})

test_that("as.vector of an age x region matrix is in backend group order", {
  m <- expand_param(matrix(1:6, 3, 2), "x", dims)
  # group (r - 1) * n_age + a
  expect_equal(as.vector(m)[(2 - 1) * 3 + 3], m["a3", "r2"])
})

test_that("combine_mixing is the Kronecker product with age fastest", {
  C_age <- matrix(c(1, 2, 3, 4), 2)
  M <- matrix(c(0.7, 0.2, 0.3, 0.8), 2)
  C <- combine_mixing(C_age, M)
  # age b = 2 in region s = 1 seen from age a = 1 in region r = 2
  expect_equal(C[(2 - 1) * 2 + 1, (1 - 1) * 2 + 2], C_age[1, 2] * M[2, 1])
})

test_that("build_model validates its inputs", {
  expect_error(build_model(N = 100, I0 = 1), "exactly one")
  expect_error(build_model(N = 100, R0 = 2, beta = 1, I0 = 1), "exactly one")
  expect_error(build_model(N = c(a = 10, b = 10), R0 = 2, I0 = 1), "C_age")
  expect_error(build_model(N = 100, R0 = 2, I0 = 1, p_hosp = 1.5), "p_hosp")
  expect_error(build_model(N = 100, R0 = 2, I0 = 200), "exceeds")
  expect_error(small_model(M_region = matrix(0.4, 2, 2)), "sum to 1")
  expect_warning(build_model(N = 100, R0 = 2), "No initial infections")
})

test_that("stratum parameters are group parameters times stratum multipliers", {
  m <- vax_model(coverage = 0.5, ve_infection = 0.6, ve_hosp = 0.5,
                 ve_transmission = 0.2, ve_death = 0.25)
  p <- m$pars
  expect_equal(p$sus[, 2], p$sus[, 1] * 0.4)
  expect_equal(p$trans[, 2], p$trans[, 1] * 0.8)
  expect_equal(p$p_hosp[, 2], p$p_hosp[, 1] * 0.5)
  expect_equal(p$p_death_hosp[, 2], p$p_death_hosp[, 1] * 0.75)
  expect_equal(p$p_death_comm[, 2], p$p_death_comm[, 1] * 0.75)
  expect_equal(rowSums(p$S_ini + p$I_ini) + p$R_ini, as.vector(m$N))
})

test_that("stratum parameters can vary by age", {
  rel <- matrix(c(1, 1, 0.5, 0.2), 2, dimnames = list(NULL, c("u", "v")))
  m <- small_model(strata = c("u", "v"), rel_hosp = rel)
  # groups are age-fastest: young A, old A, young B, old B
  expect_equal(m$pars$p_hosp[, 2], c(0.01 * 0.5, 0.1 * 0.2, 0.01 * 0.5, 0.1 * 0.2))
})

test_that("update_model changes arguments and swaps R0/beta", {
  m <- small_model()
  m2 <- do.call(update_model, c(list(m), vax_args(coverage = 0.5,
                                                   ve_infection = 0.3)))
  expect_equal(m2$pars$sus[, 2], m2$pars$sus[, 1] * 0.7)
  m3 <- update_model(m, beta = 0.1)
  expect_equal(m3$beta, 0.1)
  m4 <- update_model(m3, R0 = 3)
  expect_equal(R0(m4), 3)
  expect_error(update_model(m, not_an_arg = 1), "Unknown")
})

test_that("mixing_schedule switches the contact matrix over time", {
  C2 <- matrix(c(6, 2, 2, 4), 2)        # a different age matrix
  iso <- diag(2)
  m <- small_model(mixing_schedule = list(list(time = 30, C_age = C2),
                                          list(time = 60, M_region = iso)))
  p <- m$pars
  expect_equal(p$contact_time, c(0, 30, 60))
  expect_equal(dim(p$contact_value), c(4, 4, 3))
  # each slice is the combination in force then; the matrix not given carries over
  M0 <- matrix(c(0.9, 0.1, 0.1, 0.9), 2)
  C0 <- matrix(c(3, 1, 1, 2), 2)
  expect_equal(p$contact_value[, , 1], combine_mixing(C0, M0), ignore_attr = TRUE)
  expect_equal(p$contact_value[, , 2], combine_mixing(C2, M0), ignore_attr = TRUE)
  expect_equal(p$contact_value[, , 3], combine_mixing(C2, iso), ignore_attr = TRUE)

  # without a schedule there is a single slice
  expect_equal(small_model()$pars$n_contact_time, 1)
})

test_that("a schedule entry at time 0 replaces the starting mixing", {
  C2 <- matrix(c(6, 2, 2, 4), 2)
  a <- small_model(mixing_schedule = list(list(time = 0, C_age = C2)))
  b <- small_model(C_age = C2)
  expect_equal(a$pars$contact_value[, , 1], b$pars$contact_value[, , 1])
  expect_equal(R0(a), R0(b))
  expect_equal(a$beta, b$beta)
})

test_that("R0 uses the mixing at time 0 and R_eff follows the schedule", {
  half <- matrix(c(3, 1, 1, 2), 2) / 2
  m <- small_model(R0 = 2, mixing_schedule = list(list(time = 50, C_age = half)))
  expect_equal(R0(m), 2)
  expect_equal(R_eff(m, 0), 2, tolerance = 1e-3)
  expect_equal(R_eff(m, 60), 1, tolerance = 1e-3)   # contacts halved
})

test_that("changed mixing changes the epidemic", {
  iso <- diag(2)
  I0 <- matrix(c(10, 0, 0, 0), 2)        # seeded in region A only
  mixed <- small_model(I0 = I0)
  early <- small_model(I0 = I0, mixing_schedule = list(list(time = 0, M_region = iso)))
  late <- small_model(I0 = I0, mixing_schedule = list(list(time = 15, M_region = iso)))
  B <- function(m) {
    d <- get_output(run_model(m, t_max = 120), "cum_infections", by = "region")
    d[region == "B" & time == 120]$value
  }
  expect_equal(B(early), 0)              # never mixed, so B stays clean
  expect_gt(B(mixed), 0)
  expect_lt(B(late), B(mixed))           # isolating early helps B
})

test_that("mixing_schedule is validated", {
  C2 <- matrix(c(6, 2, 2, 4), 2)
  expect_error(small_model(mixing_schedule = list(list(C_age = C2))), "single 'time'")
  expect_error(small_model(mixing_schedule = list(list(time = 10))), "must change")
  expect_error(small_model(mixing_schedule = list(list(time = -1, C_age = C2))),
               "mixing_schedule time")
  expect_error(small_model(mixing_schedule = list(list(time = 10, Cage = C2))),
               "Unknown element")
  expect_error(small_model(mixing_schedule = list(time = 10, C_age = C2)),
               "list of lists")
  expect_error(small_model(mixing_schedule = list(list(time = 10, M_region = matrix(0.4, 2, 2)))),
               "sum to 1")
  expect_error(small_model(mixing_schedule = list(list(time = 10, C_age = matrix(1, 3, 3)))),
               "2 x 2 matrix")
})
