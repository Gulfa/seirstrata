# Next-generation matrix (NGM) calculations.
#
# With the force of infection used in the models, one infected person in
# group j, stratum w, infects on average
#
#   K[(i,v), (j,w)] = beta * min(m_i, m_j) * sus[i,v] * C[i,j] * trans[j,w]
#                     * S[i,v] / N_j * infectious_period
#
# people in group i, stratum v. Types are ordered group-fastest, i.e. as
# as.vector() of a [group, stratum] matrix.
ngm <- function(beta, contact, sus, trans, S, N_group, mult, infectious_period) {
  n_imm <- ncol(sus)
  a <- as.vector(sus * S)
  b <- as.vector(trans / N_group)
  C2 <- kronecker(matrix(1, n_imm, n_imm), contact * outer(mult, mult, pmin))
  beta * infectious_period * (a * C2) * rep(b, each = length(a))
}

spectral_radius <- function(K) {
  max(Mod(eigen(K, only.values = TRUE)$values))
}

# R0 per unit beta: everyone susceptible in stratum 1 (no immunity), no
# contact multipliers.
r0_per_beta <- function(contact, sus, trans, N_group, infectious_period) {
  S <- matrix(0, length(N_group), ncol(sus))
  S[, 1] <- N_group
  spectral_radius(ngm(1, contact, sus, trans, S, N_group,
                      rep(1, length(N_group)), infectious_period))
}

#' Reproduction numbers of a model
#'
#' `R0` is the basic reproduction number: the spectral radius of the
#' next-generation matrix for a population entirely in stratum 1 (no
#' immunity) without contact multipliers. `R_eff` uses the initial
#' distribution over the immunity strata (people in `R` cannot be infected)
#' and the contact multipliers at time `time`. It does not include later
#' changes: infections, waning, or doses given by `vacc_schedule`.
#'
#' Note that the NGM R0 gives the right *generation-to-generation* growth once
#' the epidemic has settled into the dominant eigenvector's age/region
#' distribution. If the initial infections are distributed differently (for
#' example seeded in one region), early growth can be faster or slower than
#' R0 suggests.
#'
#' @param model A model from [build_model()].
#' @param time Time at which to evaluate the contact multipliers (R_eff only).
#' @return A single number.
#' @export
R0 <- function(model) {
  p <- model$pars
  model$beta * r0_per_beta(contact_at(model, 0), p$sus, p$trans,
                           model$N_group, 1 / p$gamma)
}

# The contact matrix in force at `time`
contact_at <- function(model, time) {
  p <- model$pars
  p$contact_value[, , max(1, findInterval(time, p$contact_time))]
}

#' @rdname R0
#' @export
R_eff <- function(model, time = 0) {
  p <- model$pars
  k <- findInterval(time, p$beta_time)
  S <- p$S_ini + p$E_ini + p$I_ini
  spectral_radius(ngm(p$beta, contact_at(model, time), p$sus, p$trans, S,
                      model$N_group, p$beta_mult_value[, k], 1 / p$gamma))
}
