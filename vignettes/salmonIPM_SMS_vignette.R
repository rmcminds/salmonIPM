# Test spawner-smolt-spawner IPM with simulated data
library(salmonIPM)

# Simulate data 
set.seed(123)
N_pop <- 20
N_year <- 30
N <- N_pop * N_year

test_data <- data.frame(pop = rep(1:N_pop, each = N_year), year = rep(1:N_year, N_pop), 
                        A = 1, p_HOS = 0, F_rate = runif(N, 0.3, 0.6), B_rate = 0,
                        n_age_obs = 50, n_HW_obs = 0)

sim_out <- IPM_sim(pars = list(mu_alpha = 5, sigma_alpha = 1, 
                               mu_Rmax = 9, sigma_Rmax = 1, rho_alphaRmax = 0,
                               beta_phi_M = 0, rho_phi_M = 0.4, sigma_phi_M = 0.3, sigma_M = 0.2,
                               mu_MS = 0.03, beta_phi_MS = 0, rho_phi_MS = 0.6, sigma_phi_MS = 0.5, sigma_MS = 0.2,
                               mu_p = c(0.1, 0.5, 0.4), sigma_p = c(0.1, 0.2), R_p = diag(2),
                               sigma_gamma = c(0.2, 0.4), R_gamma = diag(2), 
                               tau_M = 0.3, tau_S = 0.5, S_init_K = 0.2), 
                   fish_data = test_data, life_cycle = "SMS", 
                   N_age = 3, max_age = 5, ages = list(M = 2), SR_fun = "BH")

# Fit model with no pooling across populations
fit_SMS_np <- salmonIPM(fish_data = sim_out$sim_dat, 
                     ages = list(M = 2), stan_model = "IPM_SMS_np", 
                     pars = c("alpha","Rmax","beta_M","rho_M","sigma_M",
                              "mu_MS","beta_MS","rho_MS","sigma_MS","s_MS",
                              "mu_p","sigma_p","R_p","p",
                              "tau_M","tau_S","p_HOS","M","S","q","LL"),
                     chains = 3, cores = 3, iter = 1500, warmup = 500, 
                     control = list(adapt_delta = 0.95, max_treedepth = 12))

print(fit_SMS_np, probs = c(0.05,0.5,0.95),
      pars = c("p","p_HOS","S","M","s_MS","q","LL"), include = FALSE)

# Fit model with hierarchical partial pooling across populations
fit_SMS_pp <- salmonIPM(fish_data = sim_out$sim_dat, 
                        ages = list(M = 2), stan_model = "IPM_SMS_pp", 
                        pars = c("mu_alpha","sigma_alpha","alpha",
                                 "mu_Rmax","sigma_Rmax","Rmax","rho_alphaRmax",
                                 "beta_phi_M","rho_phi_M","sigma_phi_M","phi_M","sigma_M",
                                 "mu_MS","beta_phi_MS","rho_phi_MS","sigma_phi_MS","phi_MS","sigma_MS","s_MS",
                                 "mu_p","sigma_gamma","R_gamma","sigma_p","R_p","p",
                                 "tau_M","tau_S","p_HOS","M","S","q","LL"),
                        chains = 3, cores = 3, iter = 2000, warmup = 1000, 
                        control = list(adapt_delta = 0.95, max_treedepth = 12))

print(fit_SMS_pp, probs = c(0.05,0.5,0.95),
      pars = c("alpha","Rmax","phi_M","phi_MS","p","p_HOS","S","M","s_MS","q","LL"), 
      include = FALSE)
