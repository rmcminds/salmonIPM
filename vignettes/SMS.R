# graphics device for this script
if(.Platform$OS.type == "windows") options(device = "windows")

#------------------------------
# Load packages
#------------------------------

## @knitr packages 
library(salmonIPM)
library(dplyr)           # data wrangling
library(tidyr)
library(distributional)  # plotting priors
library(posterior)       # working with posterior samples
library(ggplot2)         # alpha function
library(vioplot)         # posterior violin plots
library(here)            # file system paths

## @knitr unused
theme_set(theme_bw(base_size = 16))  # customize ggplot theme
theme_update(panel.grid = element_blank(),
             strip.background = element_rect(fill = NA),
             strip.text = element_text(margin = margin(b = 3, t = 3)),
             legend.background = element_blank())
library(bayesplot)      # Bayesian graphics

## @knitr multicore 
options(mc.cores = parallel::detectCores(logical = FALSE))
## @knitr

#------------------------------
# Data dimensions
#------------------------------

## @knitr data_setup
set.seed(123)
N <- 40
N_age <- 3
max_age <- 5
## @knitr

#------------------------------
# True parameter values 
#------------------------------

## @knitr pars
pars <- list(mu_alpha = 6, sigma_alpha = 0, mu_Mmax = 10, sigma_Mmax = 0, rho_alphaMmax = 0, 
             rho_M = 0.2, sigma_year_M = 0.1, sigma_M = 0,
             mu_MS = 0.02, rho_MS = 0.6, sigma_year_MS = 0.3, sigma_MS = 0,
             mu_p = c(0.05, 0.55, 0.4), sigma_pop_p = rep(0,2), R_pop_p = diag(2), 
             sigma_p = c(0.2, 0.2), R_p = matrix(c(1, 0.5, 0.5, 1), 2, 2), 
             tau_M = 0.3, tau_S = 0.2, S_init_K = 0.3)
## @knitr

#------------------------------
# Data structure
# - single population
# - habitat area
# - p_HOS
# - broodstock removal rate
# - fishing mortality
# - sample sizes
#------------------------------

## @knitr data_struct
df <- data.frame(pop = 1, year = 1:N + 2020 - N, A = 1, 
                 p_HOS = rbeta(N, 2, 2), 
                 B_rate = rbeta(N, 1, 9), 
                 F_rate = rbeta(N, 2, 3),
                 n_age_obs = runif(N, 10, 100),
                 n_HW_obs = runif(N, 10, 100))

#------------------------------
# Simulate data 
#------------------------------

## @knitr data
sim <- simIPM(life_cycle = "SMS", SR_fun = "BH", 
              N_age = N_age, max_age = max_age, ages = list(M = 2),
              pars = pars, fish_data = df)
names(sim$pars_out)
sim$pars_out[c("alpha","Mmax")]
format(head(sim$sim_dat, 10), digits = 2)
## @knitr

#-----------------------------------------------------
# Fit IPM
#-----------------------------------------------------

## @knitr fit
fit <- salmonIPM(life_cycle = "SMS", pool_pops = FALSE, 
                 SR_fun = "BH", ages = list(M = 2), fish_data = sim$sim_dat, 
                 chains = 4, iter = 2000, warmup = 1000, 
                 seed = 123)

print(fit, pars = stan_pars("IPM_SMS_np", "hyper"), prob = c(0.025, 0.5, 0.975))
## @knitr
