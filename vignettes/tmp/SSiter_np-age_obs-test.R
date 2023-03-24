#===========================================================================
# SETUP
#===========================================================================

# graphics device for this script
options(device = ifelse(.Platform$OS.type == "windows", "windows", "quartz"))

#------------------------------
# Load packages
#------------------------------

library(salmonIPM)
library(dplyr)           # data wrangling
library(tidyr)
library(matrixStats)
library(Hmisc)           # binomial CI function
library(shinystan)       # interactive exploration of posterior
library(distributional)  # plotting priors
library(posterior)       # working with posterior samples
library(ggplot2)         # alpha function
library(viridis)         # plot colors
library(vioplot)         # posterior violin plots
library(here)            # file system paths

theme_set(theme_bw(base_size = 16))  # customize ggplot theme
theme_update(panel.grid = element_blank(),
             strip.background = element_rect(fill = NA),
             strip.text = element_text(margin = margin(b = 3, t = 3)),
             legend.background = element_blank())
library(bayesplot)      # Bayesian graphics

options(mc.cores = parallel::detectCores(logical = FALSE))

#===========================================================================
# MULTIPLE POPULATIONS
#===========================================================================

#------------------------------
# Data dimensions
#------------------------------

set.seed(54321)
N_pop <- 8
N_year <- 20
N <- N_pop*N_year
N_age <- 5   # number of maiden ages
max_age <- 7 # oldest maiden spawners

#------------------------------
# True hyperparameter values 
#------------------------------

parsNpop <- list(mu_alpha = 2, sigma_alpha = 0.3, mu_Rmax = 3, sigma_Rmax = 0.3, 
                 rho_alphaRmax = 0.5, rho_R = 0.6, sigma_year_R = 0.3, sigma_R = 0.1,
                 mu_p = c(0.05, 0.4, 0.4, 0.1, 0.05), 
                 sigma_pop_p = c(0.3, 0.5, 0.5, 0.3), R_pop_p = 1 - 0.7*(1 - diag(4)), 
                 sigma_p = c(0.1, 0.2, 0.2, 0.1), R_p = 1 - 0.5*(1 - diag(4)), 
                 mu_SS = 0.1, rho_SS = 0.6, sigma_year_SS = 0.2, sigma_SS = 0.1,
                 tau = 0.3, S_init_K = 0.3)

#---------------------------------------------
# Data structure
# - habitat area
# - p_HOS
# - broodstock removal rate
# - fishing mortality
# - sample sizes (some age samples missing)
#---------------------------------------------

dfNpop <- data.frame(pop = rep(LETTERS[1:N_pop], each = N_year), 
                     year = rep(1:N_year + 2020 - N_year, N_pop),
                     A = rep(runif(N_pop, 10, 100), each = N_year), 
                     p_HOS = c(rep(0, N/2), runif(N/2, 0, 0.5)), 
                     B_rate = c(rep(0, N/2), rbeta(N/2, 1, 19)),
                     F_rate = rbeta(N, 3, 2),
                     n_age_obs = replace(runif(N, 10, 100), sample(N, N/10), 0),
                     n_HW_obs = c(rep(0, N/2), runif(N/2, 10, 100)))

#------------------------------
# Simulate data 
# - some spawner obs missing
#------------------------------

simNpop <- simIPM(life_cycle = "SSiter", SR_fun = "BH", pars = parsNpop, 
                  fish_data = dfNpop, N_age = N_age, max_age = max_age)
simNpop$sim_dat$S_obs[sample(N, N/10)] <- NA
simNpop$sim_dat <- simNpop$sim_dat %>% 
  mutate(across(starts_with("n_age"), ~ replace(.x, pop %in% LETTERS[5:8], 0)))

#-----------------------------------------------------
# IPM with no pooling across populations
#-----------------------------------------------------

fitNnp <- salmonIPM(life_cycle = "SSiter", pool_pops = FALSE, SR_fun = "BH", 
                    fish_data = simNpop$sim_dat,
                    chains = 4, iter = 2000, warmup = 1000, 
                    seed = 321)

print(fitNnp, pars = c("mu_p","sigma_p","R_p","p","s_SS","p_HOS","S","R","q"), 
      include = FALSE, prob = c(c(0.025, 0.5, 0.975)))

#-----------------------------------------------------
# IPM with partial population pooling
#-----------------------------------------------------

fitNpp <- salmonIPM(life_cycle = "SSiter", pool_pops = TRUE, SR_fun = "BH",
                    fish_data = simNpop$sim_dat,
                    chains = 4, iter = 2000, warmup = 1000, 
                    control = list(adapt_delta = 0.95, max_treedepth = 12),
                    seed = 321)

print(fitNpp, pars = c("alpha","Rmax","eta_year_R","mu_pop_alr_p","R_pop_p",
                       "R_p","p","eta_year_SS","s_SS","p_HOS","S","R","q"), 
      include = FALSE, prob = c(c(0.025, 0.5, 0.975)))

#----------------------------------------------------------
# Plot pop-level S-R parameter posteriors and true values
# - No population pooling vs. partial population pooling
#----------------------------------------------------------

windows()
par(mfcol = c(2,1), mar = c(1,4,0,1), oma = c(3,0,0,0))

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
vioplot(log(as.matrix(fitNnp, "Rmax")), 
        col = "slategray4", border = NA, drawRect = FALSE, side = "left", 
        las = 1, xlab = "", names = LETTERS[1:N_pop], ylab = "", 
        cex.axis = 1.2, cex.lab = 1.5)
vioplot(log(as.matrix(fitNpp, "Rmax")), 
        col = "salmon", border = NA, drawRect = FALSE, side = "right", add = TRUE)
points(1:N_pop, log(simNpop$pars_out$Rmax), pch = 16, cex = 1.5)
mtext(c("Population","log(Rmax)"), side = 1:2, line = 2.5, cex = 1.5)

#----------------------------------------------------------
# Plot hyperparameter posteriors, priors, and true values
# for partial population pooling model
#----------------------------------------------------------

par_names <- c("mu_alpha","mu_Rmax","sigma_year_R","rho_R","sigma_R",
               "mu_p","sigma_pop_p","rho_pop_p","sigma_p","rho_p",
               "mu_SS","rho_SS","sigma_year_SS","sigma_SS","tau")

# specify priors using distributional package
log_R_obs <- log(simNpop$sim_dat$S_obs / (1 - simNpop$sim_dat$F_rate))
priorNpop <- c(mu_alpha = dist_normal(2,5),
               mu_Rmax = dist_normal(quantile(log_R_obs, 0.8, na.rm = TRUE, names = FALSE), 
                                     sd(log_R_obs, na.rm = TRUE)),
               sigma_year_R = 2*dist_normal(0,3),
               rho_R = dist_wrap("gnorm", 0, 0.85, 50),
               sigma_R = 2*dist_normal(0,3),
               mu_p = dist_beta(1, N_age - 1),
               sigma_pop_p = 2*dist_normal(0,2),
               rho_pop_p = 2*dist_beta((N_age - 1)/2, (N_age - 1)/2) - 1, # LKJ
               sigma_p = 2*dist_normal(0,2),
               rho_p = 2*dist_beta((N_age - 1)/2, (N_age - 1)/2) - 1, # LKJ
               mu_SS = dist_uniform(0,1),
               rho_SS = dist_wrap("gnorm", 0, 0.85, 20),
               sigma_year_SS = 2*dist_normal(0,3),
               sigma_SS = 2*dist_normal(0,3),
               tau = 2*dist_normal(0,1))

# true parameter values
trueNpop <- simNpop$pars_out %>%
  c(rho_pop_p = .$R_pop_p[2,1], rho_p = .$R_p[2,1]) %>%
  .[par_names] %>% unlist() %>% setNames(gsub("(\\d)", "[\\1]", names(.)))

# extract and transform draws using posterior package
postNpop <- as_draws_rvars(fitNpp) %>%
  mutate_variables(mu_p = as.vector(mu_p),
                   sigma_pop_p = as.vector(sigma_pop_p), rho_pop_p = as.vector(R_pop_p[2,1]),
                   sigma_p = as.vector(sigma_p), rho_p = as.vector(R_p[2,1])) %>%
  as_draws_matrix(.[par_names])

# plot
windows()
par(mfrow = c(5,5), mar = c(5,1,0,1))
for(j in names(trueNpop)) {
  hist(postNpop[,j], 20, prob = TRUE, col = alpha("slategray4", 0.5), border = "white",
       xlab = j, ylab = "", yaxt = "n", main = "", cex.axis = 1.2, cex.lab = 1.5)
  curve(density(priorNpop[gsub("\\[.\\]", "", j)], at = x)[[1]], lwd = 0.5, add = TRUE)
  abline(v = trueNpop[[j]], lwd = 2, lty = 3)
}

