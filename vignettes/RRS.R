#===========================================================================
# SETUP
#===========================================================================

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
library(shinystan)       # interactive exploration of posterior
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

#===========================================================================
# SINGLE POPULATION
#===========================================================================

#------------------------------
# Data dimensions
#------------------------------

## @knitr singlepop_data_setup
set.seed(12321)
N <- 30
N_age <- 3
max_age <- 5
## @knitr

#------------------------------
# True parameter values 
#------------------------------

## @knitr singlepop_pars
pars1pop <- list(mu_alpha = 2, sigma_alpha = 0, mu_Rmax = 7, sigma_Rmax = 0, 
                 rho_alphaRmax = 0, rho_R = 0.6, sigma_year_R = 0.2, sigma_R = 0,
                 mu_p = c(0.05, 0.55, 0.4), sigma_pop_p = rep(0,2), R_pop_p = diag(2), 
                 sigma_p = c(0.2, 0.2), R_p = matrix(c(1, 0.5, 0.5, 1), 2, 2), 
                 tau = 0.3, S_init_K = 0.3)
## @knitr

#------------------------------
# Data structure
# - habitat area
# - p_HOS
# - broodstock removal rate
# - fishing mortality
# - sample sizes
#------------------------------

## @knitr singlepop_data_struct
df1pop <- data.frame(pop = 1, year = 1:N + 2020 - N, A = 1, 
                     p_HOS = runif(N/2, 0, 0.5), 
                     B_rate = 0, F_rate = rbeta(N, 3, 2),
                     n_age_obs = runif(N, 10, 100),
                     n_HW_obs = runif(N, 10, 100))

#------------------------------
# Simulate data 
#------------------------------

## @knitr singlepop_data
sim1pop <- simIPM(life_cycle = "SS", SR_fun = "BH", pars = pars1pop, 
                  fish_data = df1pop, N_age = N_age, max_age = max_age)
names(sim1pop$pars_out)
sim1pop$pars_out[c("alpha","Rmax")]
format(head(sim1pop$sim_dat, 10), digits = 2)
## @knitr

#-----------------------------------------------------
# Fit IPM
#-----------------------------------------------------

## @knitr singlepop_fit
fit1pop <- salmonIPM(life_cycle = "SS", pool_pops = FALSE, 
                     SR_fun = "BH", RRS = c("alpha","Rmax"),
                     fish_data = sim1pop$sim_dat, 
                     chains = 4, iter = 2000, warmup = 1000, 
                     seed = 123)

print(fit1pop, pars = stan_pars("IPM_SS_np", "hyper"), prob = c(0.025, 0.5, 0.975))
## @knitr

