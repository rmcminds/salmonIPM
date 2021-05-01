# Test spawner-smolt-spawner IPM with simulated data
library(salmonIPM)

# Simulate data 
set.seed(123)
N_pop <- 20
N_year <- 30
N <- N_pop * N_year

test_data <- data.frame(pop = rep(1:N_pop, each = N_year), year = rep(1:N_year, N_pop), 
                        A = 1, p_HOS = 0, F_rate = runif(N, 0.3, 0.6), B_rate = 0,
                        n_age_obs = 50, n_HW_obs = 0, x = runif(N))

sim_out <- salmonIPM_sim(pars = list(mu_alpha = 5, beta_alpha = , sigma_alpha = 1, 
                                     mu_Rmax = 9, beta_Rmax = -1, sigma_Rmax = 1, rho_alphaRmax = 0,
                                     beta_M = 0, rho_M = 0.4, sigma_year_M = 0.3, sigma_M = 0.2,
                                     mu_MS = 0.03, beta_MS = 0, rho_MS = 0.6, sigma_year_MS = 0.5, sigma_MS = 0.2,
                                     mu_p = c(0.1, 0.5, 0.4), sigma_p = c(0.1, 0.2), R_p = diag(2),
                                     sigma_pop_p = c(0.2, 0.4), R_pop_p = diag(2), 
                                     tau_M = 0.3, tau_S = 0.5, S_init_K = 0.2), 
                         fish_data = test_data, life_cycle = "SMS", 
                         N_age = 3, max_age = 5, ages = list(M = 2), SR_fun = "BH")

# Fit model with no pooling across populations
fit_SMS_np <- salmonIPM(stan_model = "IPM_SMS_np", par_models = list(Mmax ~ x),
                        ages = list(M = 2), fish_data = sim_out$sim_dat, 
                        chains = 3, cores = 3, iter = 1500, warmup = 500, 
                        control = list(adapt_delta = 0.95, max_treedepth = 12))

print(fit_SMS_np, probs = c(0.05,0.5,0.95),
      pars = c("p","p_HOS","B_rate","S","M","s_MS","q","LL"), include = FALSE)

# Fit model with hierarchical partial pooling across populations
fit_SMS_pp <- salmonIPM(stan_model = "IPM_SMS_pp", par_models = list(Mmax ~ x),
                        ages = list(M = 2), fish_data = sim_out$sim_dat, 
                        chains = 3, cores = 3, iter = 1500, warmup = 500, 
                        control = list(adapt_delta = 0.95, max_treedepth = 12))

print(fit_SMS_pp, probs = c(0.05,0.5,0.95),
      pars = c("alpha","Mmax","eta_year_M","eta_year_MS","mu_pop_alr_p","p",
               "p_HOS","B_rate","S","M","s_MS","q","LL"), include = FALSE)
