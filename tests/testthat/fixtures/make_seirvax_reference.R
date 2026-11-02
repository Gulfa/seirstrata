# Reference output from seirvax (modelling_prep_ex @64face5) for the
# equivalence tests in seirstrata. Run once; not part of the package.
suppressMessages(devtools::load_all("/fhi/modelling_prep_ex", quiet = TRUE))
N <- matrix(c(20000, 30000, 10000, 15000), 2,
            dimnames = list(c("young", "old"), c("A", "B")))
base <- list(N = N, C_age = matrix(c(3, 1, 1, 2), 2),
             M_region = matrix(c(0.9, 0.1, 0.1, 0.9), 2), R0 = 2,
             I0 = matrix(c(10, 0, 0, 0), 2), p_hosp = c(0.01, 0.1),
             p_death_hosp = c(0.02, 0.2), p_death_comm = c(0, 0.01))
vars <- c("S", "E", "I", "H", "D_hosp", "D_comm", "cum_infections", "cum_hosp",
          "inc_hosp", "cum_doses")
run_ref <- function(args) {
  m <- do.call(build_model, args)
  out <- get_output(run_model(m, t_max = 200), vars, by = c("age", "region"))
  list(R0 = R0(m), R_eff = R_eff(m), beta = m$beta, out = as.data.frame(out))
}
# 1: initial coverage with all VE types, no doses over time
a1 <- modifyList(base, list(vacc_coverage = c(0.2, 0.6), ve_infection = 0.5,
                            ve_transmission = 0.2, ve_hosp = 0.6,
                            ve_death = 0.3,
                            beta_schedule = data.frame(time = 40, mult = 0.7)))
# 2: doses over time with ve_infection = 1, so every R is unvaccinated
a2 <- modifyList(base, list(ve_infection = 1, ve_hosp = 0.5,
                            vacc_schedule = data.frame(time = c(10, 60),
                                                       doses = c(300, 0))))
ref <- list(coverage = c(run_ref(a1), list(args = a1)),
            doses = c(run_ref(a2), list(args = a2)))
saveRDS(ref, "/fhi/seirstrata/tests/testthat/fixtures/seirvax_reference.rds")
cat("saved\n")
