# A small two-age, two-region model used across tests.
small_model <- function(...) {
  N <- matrix(c(20000, 30000, 10000, 15000), 2,
              dimnames = list(c("young", "old"), c("A", "B")))
  C_age <- matrix(c(3, 1, 1, 2), 2)
  M_region <- matrix(c(0.9, 0.1, 0.1, 0.9), 2)
  args <- list(N = N, C_age = C_age, M_region = M_region, R0 = 2,
               I0 = matrix(c(10, 0, 0, 0), 2), p_hosp = c(0.01, 0.1),
               p_death_hosp = c(0.02, 0.2), p_death_comm = c(0, 0.01))
  args[names(list(...))] <- list(...)
  do.call("build_model", args)
}

total_alive_and_dead <- function(res) {
  d <- get_output(res, c("S", "E", "I", "H", "R", "D"))
  d[, .(total = sum(value)), by = .(sim, time)]$total
}

# Noise-free daily admissions from the deterministic model
det_admissions <- function(model, t_max, by = character(0)) {
  d <- get_output(run_model(model, t_max = t_max), "inc_hosp", by = by)[time > 0]
  d[, admissions := round(value)]
  d[, c("time", by, "admissions"), with = FALSE]
}

# Two strata playing the role of seirvax's unvaccinated / vaccinated: VE
# becomes 1 - the stratum multiplier, coverage the initial share in stratum 2.
vax_args <- function(coverage = 0, ve_infection = 0, ve_transmission = 0,
                     ve_hosp = 0, ve_death = 0) {
  list(strata = c("unvaccinated", "vaccinated"),
       rel_susceptibility = c(1, 1 - ve_infection),
       rel_transmissibility = c(1, 1 - ve_transmission),
       rel_hosp = c(1, 1 - ve_hosp),
       rel_death = c(1, 1 - ve_death),
       initial_strata = if (length(coverage) > 1) cbind(1 - coverage, coverage)
                        else c(1 - coverage, coverage))
}

vax_model <- function(coverage = 0, ve_infection = 0, ve_transmission = 0,
                      ve_hosp = 0, ve_death = 0, ...) {
  do.call(small_model, c(vax_args(coverage, ve_infection, ve_transmission,
                                  ve_hosp, ve_death), list(...)))
}
