# Test spawner-smolt-spawner IPM with simulated data
library(salmonIPM)

# Simulate data 
set.seed(123)
N <- 50
set.seed(123)
test_data <- data.frame(pop = rep(1,N), year =  1:N, A = 1, p_HOS = 0, 
                        F_rate = runif(N, 0.3, 0.6), B_rate = 0,
                        n_age_obs = 50, n_HW_obs = 0)
sim_out <- IPM_sim(pars = list(mu_alpha = 5, sigma_alpha = 1, mu_Rmax = 9, sigma_Rmax = 1,
                               rho_alphaRmax = 0,
                               beta_M = 0, rho_M = 0.4, sigma_M = 0.3,
                               mu_MS = 0.03, beta_MS = 0, rho_MS = 0.6, sigma_MS = 0.5,
                               mu_p = c(0.1, 0.5, 0.4), sigma_p = c(0.1, 0.2), R_p = diag(2),
                               sigma_gamma = c(0.2, 0.4), R_gamma = diag(2), 
                               tau_M = 0.3, tau_S = 0.5, 
                               S_init_K = 0.2), 
                   fish_data = test_data, life_cycle = "SMS", 
                   N_age = 3, max_age = 5, ages = list(M = 2), SR_fun = "BH")

# Fit model
ipm_fit <- salmonIPM(fish_data = sim_out$sim_dat, 
                     ages = list(M = 2), stan_model = "IPM_SMS_np", 
                     chains = 3, cores = 3, iter = 1500, warmup = 500, 
                     control = list(adapt_delta = 0.99))

print(ipm_fit, pars = c("p","B_rate_all","S","M","s_MS","q"), include = FALSE)
launch_shinystan(ipm_fit)
