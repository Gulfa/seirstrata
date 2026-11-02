## Deterministic SEIR-HD model with an immunity ladder (ODE).
##
## S, E, I and H are [n_group, n_imm] arrays: column k is immunity stratum k,
## from 1 (no immunity) up to n_imm. Each stratum has its own relative
## susceptibility, infectiousness and severity (built in R). An infection
## keeps the stratum of the person infected, so breakthrough infections use
## that stratum's infectiousness and severity.
##
## R is a single compartment per group: recovered people are fully protected
## until their immunity wanes, when they enter the ladder at stratum
## R_to_stratum. The model knows nothing about age or region: the R layer
## (build_model) flattens age x region into n_group groups and builds the
## n_group x n_group contact matrix.
##
##   S[k] -> E[k] -> I[k] -> R
##                   I[k] -> H[k] -> R
##                   I[k] -> D_comm   (at the end of the infectious period)
##                           H[k] -> D_hosp
##   S[k] -> S[k - 1]               (waning down the ladder)
##   R    -> S[R_to_stratum]        (waning of infection-acquired immunity)
##   S[k] -> S[vacc_stratum]        (vaccination, for k < vacc_stratum)
##
## Vaccination status is not tracked: doses go to people in S below the
## vaccine stratum and to people in R (where they are wasted), in proportion.
##
## The parameter interface and the state declarations must stay identical to
## seir_stoch.R (see docs/modifying_the_model.md).

## ---- Parameters -----------------------------------------------------------
n_group     <- parameter(type = "integer", constant = TRUE)
n_imm       <- parameter(type = "integer", constant = TRUE)
n_beta_time <- parameter(type = "integer", constant = TRUE)
n_beta_t    <- parameter(type = "integer", constant = TRUE)

beta  <- parameter()   # transmission rate (calibrated from R0 in R)
sigma <- parameter()   # 1 / latent period
gamma <- parameter()   # 1 / infectious period

# contact[i, j]: contacts per day a person in group i has with group j.
# Piecewise constant in time, so the mixing pattern can change (school
# holidays, travel restrictions); the last dimension of contact_value is time.
n_contact_time <- parameter(type = "integer", constant = TRUE)
contact_time   <- parameter()
contact_value  <- parameter()
dim(contact_time)  <- n_contact_time
dim(contact_value) <- c(n_group, n_group, n_contact_time)
contact <- interpolate(contact_time, contact_value, "constant")
dim(contact) <- c(n_group, n_group)

# [group, stratum] parameters
sus          <- parameter()   # relative susceptibility
trans        <- parameter()   # relative infectiousness when infected
p_hosp       <- parameter()   # P(hospitalised | infected)
p_death_comm <- parameter()   # P(community death | infected, not hospitalised)
p_death_hosp <- parameter()   # P(death | hospitalised)
waning_rate  <- parameter()   # rate of moving from stratum k to k - 1
dim(sus, trans, p_hosp, p_death_comm, p_death_hosp, waning_rate) <-
  c(n_group, n_imm)

dur_hosp      <- parameter()  # mean length of hospital stay (days)
waning_rate_R <- parameter()  # rate of leaving R
dim(dur_hosp, waning_rate_R) <- n_group
# 0/1 masks over strata (built in R): which strata vaccination takes people
# from, which stratum it moves them to, and which stratum R wanes into.
vacc_from <- parameter()
vacc_to   <- parameter()
R_to      <- parameter()
dim(vacc_from, vacc_to, R_to) <- n_imm

S_ini <- parameter()
E_ini <- parameter()
I_ini <- parameter()
dim(S_ini, E_ini, I_ini) <- c(n_group, n_imm)
R_ini <- parameter()
dim(R_ini) <- n_group

# Piecewise-constant contact multiplier per group: the fraction of the group's
# contacts that still happen. Contacts between groups i and j are scaled by
# min(beta_mult[i], beta_mult[j]) -- the more restricted group decides -- so a
# multiplier m applied to everyone scales all transmission by m.
beta_time       <- parameter()
beta_mult_value <- parameter()
dim(beta_time)       <- n_beta_time
dim(beta_mult_value) <- c(n_group, n_beta_time)
beta_mult <- interpolate(beta_time, beta_mult_value, "constant")
dim(beta_mult) <- n_group

# Global transmission multiplier over time, e.g. seasonality and holidays
# built outside the package. Linear between knots (daily values change
# every day, and a piecewise-constant input would make the solver stop at
# every day boundary); multiplies beta for everyone, on top of the targeted
# contact multipliers above. The R side adds a final knot so the last value
# holds.
beta_t_time  <- parameter()
beta_t_value <- parameter()
dim(beta_t_time, beta_t_value) <- n_beta_t
beta_t <- interpolate(beta_t_time, beta_t_value, "linear")

# Vaccination over time: doses per day given to each group, linear between
# knots (built in R from vacc_schedule: each change of the schedule is a
# one-day ramp, so totals are exact and the solver is not slowed by steps).
n_vacc_time      <- parameter(type = "integer", constant = TRUE)
vacc_time        <- parameter()
vacc_doses_value <- parameter()
dim(vacc_time)        <- n_vacc_time
dim(vacc_doses_value) <- c(n_group, n_vacc_time)
vacc_doses <- interpolate(vacc_time, vacc_doses_value, "linear")
dim(vacc_doses) <- n_group

## ---- Force of infection ---------------------------------------------------
N[] <- sum(S[i, ]) + sum(E[i, ]) + sum(I[i, ]) + sum(H[i, ]) + R[i]
dim(N) <- n_group

# Infectious pressure exerted by group i (per member of group i)
trans_I[, ] <- trans[i, j] * I[i, j]
dim(trans_I) <- c(n_group, n_imm)
prev[] <- sum(trans_I[i, ]) / N[i]
dim(prev) <- n_group

contact_prev[, ] <- contact[i, j] * min(beta_mult[i], beta_mult[j]) * prev[j]
dim(contact_prev) <- c(n_group, n_group)

lambda[] <- beta * beta_t * sum(contact_prev[i, ])
dim(lambda) <- n_group

## ---- Vaccination ----------------------------------------------------------
# Per-capita rate among the eligible that delivers vacc_doses per day. Eligible
# are S below the vaccine stratum and R (whose doses are wasted). Capped at
# 10 per day, so a group whose eligible people run out is vaccinated to
# near-exhaustion within hours rather than going negative.
S_elig[, ] <- vacc_from[j] * S[i, j]
dim(S_elig) <- c(n_group, n_imm)
vacc_eligible[] <- sum(S_elig[i, ]) + R[i]
vacc_rate[] <- min(vacc_doses[i] / max(vacc_eligible[i], 1e-10), 10)
dim(vacc_eligible, vacc_rate) <- n_group

## ---- Flows ----------------------------------------------------------------
n_SE[, ] <- lambda[i] * sus[i, j] * S[i, j]
n_S_wane[, ] <- waning_rate[i, j] * S[i, j]   # waning_rate[, 1] = 0
n_SV[, ] <- vacc_from[j] * vacc_rate[i] * S[i, j]
n_EI[, ] <- sigma * E[i, j]
n_I_out[, ] <- gamma * I[i, j]
n_IH[, ] <- n_I_out[i, j] * p_hosp[i, j]
n_ID[, ] <- (n_I_out[i, j] - n_IH[i, j]) * p_death_comm[i, j]
n_H_out[, ] <- H[i, j] / dur_hosp[i]
n_HD[, ] <- n_H_out[i, j] * p_death_hosp[i, j]
# recoveries, into R
n_IR[, ] <- n_I_out[i, j] - n_IH[i, j] - n_ID[i, j]
n_HR[, ] <- n_H_out[i, j] - n_HD[i, j]
dim(n_SE, n_S_wane, n_SV, n_EI, n_I_out, n_IH, n_ID, n_H_out, n_HD, n_IR,
    n_HR) <- c(n_group, n_imm)

n_R_wane[] <- waning_rate_R[i] * R[i]
n_R_vacc[] <- vacc_rate[i] * R[i]   # doses given to R: counted, no effect
vacc_total[] <- sum(n_SV[i, ])
dim(n_R_wane, n_R_vacc, vacc_total) <- n_group

# arriving in S[, k]: from the stratum above, from vaccination, from R
n_S_down[, ] <- if (j < n_imm) n_S_wane[i, j + 1] else 0
n_S_vacc_in[, ] <- vacc_to[j] * vacc_total[i]
n_S_R_in[, ] <- R_to[j] * n_R_wane[i]
dim(n_S_down, n_S_vacc_in, n_S_R_in) <- c(n_group, n_imm)

## ---- States ---------------------------------------------------------------
initial(S[, ]) <- S_ini[i, j]
initial(E[, ]) <- E_ini[i, j]
initial(I[, ]) <- I_ini[i, j]
initial(H[, ]) <- 0
initial(D_hosp[, ]) <- 0
initial(D_comm[, ]) <- 0
initial(cum_infections[, ]) <- 0
initial(cum_hosp[, ]) <- 0
# Daily incidence: flow during (t - 1, t]
initial(inc_infections[, ], zero_every = 1) <- 0
initial(inc_hosp[, ], zero_every = 1) <- 0
initial(inc_deaths_hosp[, ], zero_every = 1) <- 0
initial(inc_deaths_comm[, ], zero_every = 1) <- 0
dim(S, E, I, H, D_hosp, D_comm, cum_infections, cum_hosp,
    inc_infections, inc_hosp, inc_deaths_hosp, inc_deaths_comm) <-
  c(n_group, n_imm)

initial(R[]) <- R_ini[i]
# Doses given by vacc_schedule (to S and R alike)
initial(cum_doses[]) <- 0
initial(inc_doses[], zero_every = 1) <- 0
# Weekly admissions: flow during (t - 7, t] at t = 7, 14, ...; model day 0 is
# the first day of a reporting week
initial(inc_hosp_week[], zero_every = 7) <- 0
dim(R, cum_doses, inc_doses, inc_hosp_week) <- n_group

deriv(S[, ]) <- -n_SE[i, j] - n_S_wane[i, j] - n_SV[i, j] + n_S_down[i, j] +
  n_S_vacc_in[i, j] + n_S_R_in[i, j]
deriv(E[, ]) <- n_SE[i, j] - n_EI[i, j]
deriv(I[, ]) <- n_EI[i, j] - n_I_out[i, j]
deriv(H[, ]) <- n_IH[i, j] - n_H_out[i, j]
deriv(R[]) <- sum(n_IR[i, ]) + sum(n_HR[i, ]) - n_R_wane[i]
deriv(D_hosp[, ]) <- n_HD[i, j]
deriv(D_comm[, ]) <- n_ID[i, j]
deriv(cum_infections[, ]) <- n_SE[i, j]
deriv(cum_hosp[, ]) <- n_IH[i, j]
deriv(inc_infections[, ]) <- n_SE[i, j]
deriv(inc_hosp[, ]) <- n_IH[i, j]
deriv(inc_deaths_hosp[, ]) <- n_HD[i, j]
deriv(inc_deaths_comm[, ]) <- n_ID[i, j]
deriv(cum_doses[]) <- sum(n_SV[i, ]) + n_R_vacc[i]
deriv(inc_doses[]) <- sum(n_SV[i, ]) + n_R_vacc[i]
deriv(inc_hosp_week[]) <- sum(n_IH[i, ])

## ---- Comparison to data (used only when fitting) --------------------------
# Observed hospital admissions (daily or weekly), as n_obs series.
# obs_map[k, g] is the weight of group g in series k (built in R from the
# data's region / age columns; by default one national series).
n_obs  <- parameter(type = "integer", constant = TRUE)
obs_map <- parameter()
dim(obs_map) <- c(n_obs, n_group)
nb_size <- parameter()   # negative binomial dispersion (larger = less noise)
# 1: compare weekly admissions (data at days 7, 14, ...); 0: daily
obs_weekly <- parameter()

adm_group[] <- obs_weekly * inc_hosp_week[i] +
  (1 - obs_weekly) * sum(inc_hosp[i, ])
dim(adm_group) <- n_group
adm_contrib[, ] <- obs_map[i, j] * adm_group[j]
dim(adm_contrib) <- c(n_obs, n_group)
model_admissions[] <- sum(adm_contrib[i, ])
dim(model_admissions) <- n_obs

admissions <- data()
dim(admissions) <- n_obs
admissions[] ~ NegativeBinomial(size = nb_size, mu = model_admissions[i] + 1e-6)
