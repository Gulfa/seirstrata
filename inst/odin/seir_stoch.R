## Stochastic SEIR-HD model with an immunity ladder (tau-leaping).
##
## Discrete-time version of seir_ode.R: the same parameters and states, with
## binomial draws over a step of `dt` days (set when the system is created).
## Outflows from S, I, H and R are split with chained binomials, so the
## competing outcomes of a compartment never exceed the number of people
## leaving it. See seir_ode.R for the model structure.
##
## The parameter interface and the state declarations must stay identical to
## seir_ode.R (see docs/modifying_the_model.md).

## ---- Parameters -----------------------------------------------------------
n_group     <- parameter(type = "integer", constant = TRUE)
n_imm       <- parameter(type = "integer", constant = TRUE)
n_beta_time <- parameter(type = "integer", constant = TRUE)

beta  <- parameter()   # transmission rate (calibrated from R0 in R)
sigma <- parameter()   # 1 / latent period
gamma <- parameter()   # 1 / infectious period

# contact[i, j]: contacts per day a person in group i has with group j.
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

# Piecewise-constant contact multiplier per group (see seir_ode.R)
beta_time       <- parameter()
beta_mult_value <- parameter()
dim(beta_time)       <- n_beta_time
dim(beta_mult_value) <- c(n_group, n_beta_time)
beta_mult <- interpolate(beta_time, beta_mult_value, "constant")
dim(beta_mult) <- n_group

# Vaccination over time: doses per day given to each group
n_vacc_time      <- parameter(type = "integer", constant = TRUE)
vacc_time        <- parameter()
vacc_doses_value <- parameter()
dim(vacc_time)        <- n_vacc_time
dim(vacc_doses_value) <- c(n_group, n_vacc_time)
vacc_doses <- interpolate(vacc_time, vacc_doses_value, "constant")
dim(vacc_doses) <- n_group

## ---- Force of infection ---------------------------------------------------
N[] <- sum(S[i, ]) + sum(E[i, ]) + sum(I[i, ]) + sum(H[i, ]) + R[i]
dim(N) <- n_group

trans_I[, ] <- trans[i, j] * I[i, j]
dim(trans_I) <- c(n_group, n_imm)
prev[] <- sum(trans_I[i, ]) / N[i]
dim(prev) <- n_group

contact_prev[, ] <- contact[i, j] * min(beta_mult[i], beta_mult[j]) * prev[j]
dim(contact_prev) <- c(n_group, n_group)

lambda[] <- beta * sum(contact_prev[i, ])
dim(lambda) <- n_group

## ---- Vaccination ----------------------------------------------------------
S_elig[, ] <- vacc_from[j] * S[i, j]
dim(S_elig) <- c(n_group, n_imm)
vacc_eligible[] <- sum(S_elig[i, ]) + R[i]
vacc_rate[] <- min(vacc_doses[i] / max(vacc_eligible[i], 1e-10), 10)
dim(vacc_eligible, vacc_rate) <- n_group

## ---- Flows (binomial draws over one step of length dt) ---------------------
p_EI <- 1 - exp(-sigma * dt)
p_I_out <- 1 - exp(-gamma * dt)

# S is left by infection, waning or vaccination, split with chained binomials
rate_SE[, ] <- lambda[i] * sus[i, j]
rate_SW[, ] <- waning_rate[i, j]                 # waning_rate[, 1] = 0
rate_SV[, ] <- vacc_from[j] * vacc_rate[i]
p_S_out[, ] <- 1 - exp(-(rate_SE[i, j] + rate_SW[i, j] + rate_SV[i, j]) * dt)
p_SE_given_out[, ] <- rate_SE[i, j] /
  max(rate_SE[i, j] + rate_SW[i, j] + rate_SV[i, j], 1e-30)
p_SW_given_not_SE[, ] <- rate_SW[i, j] / max(rate_SW[i, j] + rate_SV[i, j], 1e-30)
p_H_out[] <- 1 - exp(-dt / dur_hosp[i])
dim(rate_SE, rate_SW, rate_SV, p_S_out, p_SE_given_out, p_SW_given_not_SE) <-
  c(n_group, n_imm)
dim(p_H_out) <- n_group

n_S_out[, ] <- Binomial(S[i, j], p_S_out[i, j])
n_SE[, ] <- Binomial(n_S_out[i, j], p_SE_given_out[i, j])
n_S_wane[, ] <- Binomial(n_S_out[i, j] - n_SE[i, j], p_SW_given_not_SE[i, j])
n_SV[, ] <- n_S_out[i, j] - n_SE[i, j] - n_S_wane[i, j]
n_EI[, ] <- Binomial(E[i, j], p_EI)
n_I_out[, ] <- Binomial(I[i, j], p_I_out)
n_IH[, ] <- Binomial(n_I_out[i, j], p_hosp[i, j])
n_ID[, ] <- Binomial(n_I_out[i, j] - n_IH[i, j], p_death_comm[i, j])
n_H_out[, ] <- Binomial(H[i, j], p_H_out[i])
n_HD[, ] <- Binomial(n_H_out[i, j], p_death_hosp[i, j])
n_IR[, ] <- n_I_out[i, j] - n_IH[i, j] - n_ID[i, j]
n_HR[, ] <- n_H_out[i, j] - n_HD[i, j]
dim(n_S_out, n_SE, n_S_wane, n_SV, n_EI, n_I_out, n_IH, n_ID, n_H_out, n_HD,
    n_IR, n_HR) <- c(n_group, n_imm)

# R is left by waning; doses given to R are counted but change nothing
p_R_out[] <- 1 - exp(-(waning_rate_R[i] + vacc_rate[i]) * dt)
p_RW_given_out[] <- waning_rate_R[i] / max(waning_rate_R[i] + vacc_rate[i], 1e-30)
n_R_out[] <- Binomial(R[i], p_R_out[i])
n_R_wane[] <- Binomial(n_R_out[i], p_RW_given_out[i])
n_R_vacc[] <- n_R_out[i] - n_R_wane[i]
vacc_total[] <- sum(n_SV[i, ])
dim(p_R_out, p_RW_given_out, n_R_out, n_R_wane, n_R_vacc, vacc_total) <- n_group

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
initial(cum_doses[]) <- 0
initial(inc_doses[], zero_every = 1) <- 0
dim(R, cum_doses, inc_doses) <- n_group

update(S[, ]) <- S[i, j] - n_S_out[i, j] + n_S_down[i, j] + n_S_vacc_in[i, j] +
  n_S_R_in[i, j]
update(E[, ]) <- E[i, j] + n_SE[i, j] - n_EI[i, j]
update(I[, ]) <- I[i, j] + n_EI[i, j] - n_I_out[i, j]
update(H[, ]) <- H[i, j] + n_IH[i, j] - n_H_out[i, j]
update(R[]) <- R[i] + sum(n_IR[i, ]) + sum(n_HR[i, ]) - n_R_wane[i]
update(D_hosp[, ]) <- D_hosp[i, j] + n_HD[i, j]
update(D_comm[, ]) <- D_comm[i, j] + n_ID[i, j]
update(cum_infections[, ]) <- cum_infections[i, j] + n_SE[i, j]
update(cum_hosp[, ]) <- cum_hosp[i, j] + n_IH[i, j]
update(inc_infections[, ]) <- inc_infections[i, j] + n_SE[i, j]
update(inc_hosp[, ]) <- inc_hosp[i, j] + n_IH[i, j]
update(inc_deaths_hosp[, ]) <- inc_deaths_hosp[i, j] + n_HD[i, j]
update(inc_deaths_comm[, ]) <- inc_deaths_comm[i, j] + n_ID[i, j]
update(cum_doses[]) <- cum_doses[i] + vacc_total[i] + n_R_vacc[i]
update(inc_doses[]) <- inc_doses[i] + vacc_total[i] + n_R_vacc[i]

## ---- Comparison to data (used only when fitting) --------------------------
n_obs  <- parameter(type = "integer", constant = TRUE)
obs_map <- parameter()
dim(obs_map) <- c(n_obs, n_group)
nb_size <- parameter()   # negative binomial dispersion (larger = less noise)

adm_contrib[, ] <- obs_map[i, j] * sum(inc_hosp[j, ])
dim(adm_contrib) <- c(n_obs, n_group)
model_admissions[] <- sum(adm_contrib[i, ])
dim(model_admissions) <- n_obs

admissions <- data()
dim(admissions) <- n_obs
admissions[] ~ NegativeBinomial(size = nb_size, mu = model_admissions[i] + 1e-6)
