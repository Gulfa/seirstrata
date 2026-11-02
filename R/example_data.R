#' Made-up example data
#'
#' Illustrative (not real) inputs for four age groups and three regions,
#' used in the examples and tests. The numbers are loosely plausible but
#' invented; replace them with real data for any actual analysis.
#'
#' @return A list with
#' * `N`: population, 4 x 3 (age x region) matrix.
#' * `C_age`: contacts per day by age, 4 x 4. Reciprocal with respect to the
#'   total population by age (`N_a C_age[a, b] = N_b C_age[b, a]`).
#' * `M_region`: fraction of contacts made in each region, 3 x 3, rows sum
#'   to 1.
#' * `epi`: data.frame of age-specific parameters (`susceptibility`,
#'   `p_hosp`, `p_death_hosp`, `p_death_comm`, `vacc_coverage`).
#' @export
#' @examples
#' ex <- example_data()
#' m <- build_model(N = ex$N, C_age = ex$C_age, M_region = ex$M_region,
#'                  R0 = 2, I0 = 5,
#'                  susceptibility = ex$epi$susceptibility,
#'                  p_hosp = ex$epi$p_hosp)
example_data <- function() {
  age <- c("0-19", "20-44", "45-64", "65+")
  region <- c("North", "Central", "South")

  N <- matrix(c(60000, 90000, 70000, 50000,
                150000, 250000, 160000, 90000,
                80000, 110000, 90000, 70000),
              nrow = 4, dimnames = list(age, region))

  contacts <- matrix(c(8.0, 2.5, 1.5, 0.7,
                       3.0, 6.0, 3.0, 1.2,
                       1.5, 2.5, 4.0, 1.5,
                       0.5, 0.7, 1.0, 2.5),
                     nrow = 4, dimnames = list(age, age))
  # Make reciprocal: total contacts between a and b are the same both ways
  N_age <- rowSums(N)
  total <- contacts * N_age
  C_age <- (total + t(total)) / 2 / N_age

  M_region <- matrix(c(0.85, 0.04, 0.03,
                       0.12, 0.90, 0.12,
                       0.03, 0.06, 0.85),
                     nrow = 3, dimnames = list(region, region))

  epi <- data.frame(
    age = age,
    susceptibility = c(0.7, 1, 1, 1),
    p_hosp = c(0.002, 0.01, 0.03, 0.12),
    p_death_hosp = c(0.01, 0.03, 0.08, 0.25),
    p_death_comm = c(0, 0.0001, 0.0005, 0.01),
    vacc_coverage = c(0, 0.3, 0.6, 0.9))

  list(N = N, C_age = C_age, M_region = M_region, epi = epi)
}
