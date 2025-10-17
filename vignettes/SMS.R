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

#===========================================================================
# SINGLE POPULATION
#===========================================================================

#------------------------------
# Data dimensions
#------------------------------

## @knitr singlepop_data_setup
set.seed(123)
N <- 40
N_age <- 3
max_age <- 5
## @knitr

#------------------------------
# True parameter values
#------------------------------

## @knitr singlepop_pars
pars1pop <- list(mu_alpha = 6, sigma_alpha = 0, mu_Mmax = 10, sigma_Mmax = 0,
                 rho_alphaMmax = 0, rho_M = 0.2, sigma_year_M = 0.1, sigma_M = 0,
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

## @knitr singlepop_data_struct
df1pop <- data.frame(pop = 1, year = 1:N + 2020 - N, A = 1,
                     p_HOS = rbeta(N, 2, 2),
                     B_rate = rbeta(N, 1, 9),
                     F_rate = rbeta(N, 2, 3),
                     n_age_obs = runif(N, 10, 100),
                     n_HW_obs = runif(N, 10, 100))

#------------------------------
# Simulate data
#------------------------------

## @knitr singlepop_data
sim1pop <- simIPM(life_cycle = "SMS", SR_fun = "BH",
                  N_age = N_age, max_age = max_age, ages = list(M = 2),
                  pars = pars1pop, fish_data = df1pop)
names(sim1pop$pars_out)
sim1pop$pars_out[c("alpha","Mmax")]
format(head(sim1pop$sim_dat, 10), digits = 2)
## @knitr

#-----------------------------------------------------
# Fit IPM
#-----------------------------------------------------

## @knitr singlepop_fit
fit1pop <- salmonIPM(life_cycle = "SMS", SR_fun = "BH", ages = list(M = 2),
                     fish_data = sim1pop$sim_dat, seed = 123)

print(fit1pop)
## @knitr


#===========================================================================
# MULTIPLE POPULATIONS
#===========================================================================

#------------------------------
# Data dimensions
#------------------------------

## @knitr multipop_data_setup
N_pop <- 8
N_year <- 20
N <- N_pop*N_year
## @knitr

#------------------------------
# True hyperparameter values
#------------------------------

## @knitr multipop_pars
parsNpop <- list(mu_alpha = 6, sigma_alpha = 0.3, mu_Mmax = 10, sigma_Mmax = 0.3,
                 rho_alphaMmax = 0.5, rho_M = 0.2, sigma_year_M = 0.2, sigma_M = 0.1,
                 mu_MS = 0.02, rho_MS = 0.6, sigma_year_MS = 0.3, sigma_MS = 0.1,
                 mu_p = c(0.05, 0.55, 0.4), sigma_pop_p = c(0.2, 0.3),
                 R_pop_p = matrix(c(1, 0.3, 0.3, 1), 2, 2),
                 sigma_p = c(0.3, 0.5), R_p = matrix(c(1, 0.5, 0.5, 1), 2, 2),
                 tau_M = 0.3, tau_S = 0.2, S_init_K = 0.3)
## @knitr

#---------------------------------------------
# Data structure
# - habitat area
# - p_HOS
# - broodstock removal rate
# - fishing mortality
# - sample sizes (some age samples missing)
#---------------------------------------------

## @knitr multipop_data_struct
dfNpop <- data.frame(pop = rep(LETTERS[1:N_pop], each = N_year),
                     year = rep(1:N_year + 2020 - N_year, N_pop),
                     A = rep(runif(N_pop, 10, 100), each = N_year),
                     p_HOS = rbeta(N, 2, 2),
                     B_rate = rbeta(N, 1, 9),
                     F_rate = rbeta(N, 2, 3),
                     n_age_obs = runif(N, 10, 100),
                     n_HW_obs = runif(N, 10, 100))
## @knitr

#------------------------------
# Simulate data
# - some spawner obs missing
#------------------------------

## @knitr multipop_data
simNpop <- simIPM(life_cycle = "SMS", SR_fun = "BH",
                  N_age = N_age, max_age = max_age, ages = list(M = 2),
                  pars = parsNpop, fish_data = dfNpop)
names(simNpop$pars_out)
simNpop$pars_out[c("alpha","Mmax")]
format(head(simNpop$sim_dat, 10), digits = 2)
## @knitr

#-----------------------------------------------------
# IPM with no pooling across populations
#-----------------------------------------------------

## @knitr fit_np
fitNnp <- salmonIPM(life_cycle = "SMS", pool_pops = FALSE, SR_fun = "BH",
                    ages = list(M = 2), fish_data = simNpop$sim_dat,
                    seed = 321)

print(fitNnp)
## @knitr

#-----------------------------------------------------
# IPM with partial population pooling
#-----------------------------------------------------

## @knitr fit_pp
fitNpp <- salmonIPM(life_cycle = "SMS", SR_fun = "BH", ages = list(M = 2),
                    fish_data = simNpop$sim_dat, seed = 321)

print(fitNpp)
## @knitr

#----------------------------------------------------------
# Plot pop-level S-R parameter posteriors and true values
# - No population pooling vs. partial population pooling
#----------------------------------------------------------

## @knitr multipop_np_vs_pp
par(mfcol = c(2,1), mar = c(1,4,0,1), oma = c(3,1,0,0))

# intrinsic productivity
vioplot(log(as.matrix(fitNnp, "alpha")),
        col = "slategray4", border = NA, drawRect = FALSE, side = "left",
        las = 1, xlab = "", names = NA, ylab = "", cex.axis = 1.2, cex.lab = 1.5)
vioplot(log(as.matrix(fitNpp, "alpha")),
        col = "salmon", border = NA, drawRect = FALSE, side = "right", add = TRUE)
points(1:N_pop, log(simNpop$pars_out$alpha), pch = 16, cex = 1.5)
mtext("log(alpha)", side = 2, line = 2.5, cex = 1.5)
legend("top", c("true","np","pp"), cex = 1.2, bty = "n", horiz = TRUE,
       pch = c(16,NA,NA), pt.cex = 1.5, fill = c(NA,"slategray4","salmon"), border = NA)

# maximum recruitment
vioplot(log(as.matrix(fitNnp, "Mmax")),
        col = "slategray4", border = NA, drawRect = FALSE, side = "left",
        las = 1, xlab = "", names = LETTERS[1:N_pop], ylab = "",
        cex.axis = 1.2, cex.lab = 1.5)
vioplot(log(as.matrix(fitNpp, "Mmax")),
        col = "salmon", border = NA, drawRect = FALSE, side = "right", add = TRUE)
points(1:N_pop, log(simNpop$pars_out$Mmax), pch = 16, cex = 1.5)
mtext(c("Population","log(Mmax)"), side = 1:2, line = 2.5, cex = 1.5)
## @knitr

