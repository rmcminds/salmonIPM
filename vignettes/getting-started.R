#===========================================================================
# SETUP
#===========================================================================

# graphics device for this script
if(.Platform$OS.type == "windows") options(device = "windows")

#------------------------------
# Installing salmonIPM
#------------------------------

## @knitr install
if(!require("devtools")) install.packages("devtools")
try(devtools::install_github("ebuhle/salmonIPM", auth_token = "my_PAT"))

#------------------------------
# Load packages
#------------------------------

## @knitr packages 
library(salmonIPM)
library(dplyr)           # data wrangling
library(tidyr)
library(posterior)       # working with posterior samples
library(ggplot2)         # plotting
library(distributional)  # plotting priors
library(ggdist)
library(ggh4x)
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
df1pop <- data.frame(pop = 1, year = 1:N + 2020 - N,
                     A = 1, p_HOS = 0, B_rate = 0,
                     F_rate = rbeta(N, 3, 2),
                     n_age_obs = runif(N, 10, 100),
                     n_HW_obs = 0)

#------------------------------
# Simulate data 
#------------------------------

## @knitr singlepop_data
sim1pop <- simIPM(life_cycle = "SS", SR_fun = "BH", 
                  N_age = N_age, max_age = max_age, 
                  pars = pars1pop, fish_data = df1pop)
names(sim1pop$pars_out)
sim1pop$pars_out[c("alpha","Rmax")]
format(head(sim1pop$sim_dat, 10), digits = 2)
## @knitr

#-----------------------------------------------------
# Fit IPM
#-----------------------------------------------------

## @knitr singlepop_fit
fit1pop <- salmonIPM(life_cycle = "SS", pool_pops = FALSE, SR_fun = "BH", 
                     fish_data = sim1pop$sim_dat, 
                     chains = 4, iter = 2000, warmup = 1000, 
                     seed = 123)

print(fit1pop, pars = stan_pars("IPM_SS_np", "hyper"), prob = c(0.025, 0.5, 0.975))
## @knitr

#-----------------------------------------------------
# Plot posteriors, priors, and true values
#-----------------------------------------------------

## @knitr singlepop_posteriors
par_names <- c("log(alpha)","log(Rmax)","rho_R","sigma_R","mu_p","sigma_p","rho_p","tau")

# extract and transform draws using posterior package
post1pop <- as_draws_rvars(fit1pop) %>% 
  mutate_variables(`log(alpha)` = log(alpha), `log(Rmax)` = log(Rmax),
                   mu_p = as.vector(mu_p), sigma_p = as.vector(sigma_p), 
                   rho_p = as.vector(R_p[1,2,1])) %>% 
  .[par_names]

postdf1pop <- post1pop %>% as_draws_df() %>% as.data.frame() %>% 
  pivot_longer(cols = !starts_with("."), names_to = "param_indx", values_to = "value") %>% 
  mutate(param_indx = factor(param_indx, levels = unique(param_indx)))

# specify priors using distributional package
log_S_obs <- log(sim1pop$sim_dat$S_obs)
prior1pop <- c(`log(alpha)` = dist_normal(2,2),
               `log(Rmax)` = dist_normal(max(log_S_obs), sd(log_S_obs)),
               rho_R = dist_wrap("gnorm", 0, 0.85, 20),
               sigma_R = dist_normal(0,5),
               mu_p = dist_beta(1, N_age - 1),
               sigma_p = dist_normal(0,5),
               rho_p = 2*dist_beta((N_age - 1)/2, (N_age - 1)/2) - 1, # LKJ
               tau = dist_wrap("gnorm", 1, 0.85, 30))

priordf1pop <- data.frame(param_indx = unique(postdf1pop$param_indx)) %>% 
  mutate(param = gsub("\\[.\\]", "", param_indx), prior = prior1pop[param])

# true parameter values
true1pop <- sim1pop$pars_out %>% replace("sigma_R", .$sigma_year_R) %>% 
  c(`log(alpha)` = log(.$alpha), `log(Rmax)` = log(.$Rmax), rho_p = .$R_p[2,1]) %>% 
  .[par_names] %>% unlist() %>% setNames(gsub("(\\d)", "[\\1]", names(.))) %>% 
  data.frame(param_indx = factor(names(.), levels = names(.)), value = .)

# plot
scales <- post1pop %>% as_draws_df() %>% as.data.frame() %>% select(!starts_with(".")) %>% 
  lapply(FUN = function(x) scale_x_continuous(limits = range(x)))

postdf1pop %>%   
  ggplot(aes(x = value)) + 
  geom_histogram(aes(y = after_stat(density)), color = "white", fill = "slategray4", alpha = 0.5) +
  stat_slab(data = priordf1pop, aes(xdist = prior), inherit.aes = FALSE,
            normalize = "none", col = "black", lwd = 0.8, fill = NA) +
  geom_vline(data = true1pop, aes(xintercept = value), lwd = 0.8, lty = 2) +
  facet_wrap(~ param_indx, scales = "free", strip.position = "bottom") + 
  facetted_pos_scales(x = scales) + theme_classic() + 
  theme(strip.placement = "outside", strip.background = element_blank(), 
        strip.text = element_text(size = 11, margin = margin(b = 3, t = 0)), 
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text = element_text(size = 10), axis.text.y = element_blank()) +
  labs(x = "", y = "")

#-----------------------------------------------------
# Plot true S-R curve, obs, states and fitted draws
#-----------------------------------------------------

## @knitr singlepop_SR_base
SR <- as_draws_rvars(as.matrix(fit1pop, c("S","R")))
RR <- run_recon(sim1pop$sim_dat)
SRdat <- cbind(RR, S_true = sim1pop$pars_out$S, R_true = sim1pop$pars_out$R,
               S = SR$S, R = SR$R)
alphaRmax <- as.data.frame(fit1pop, c("alpha", "Rmax")) %>% 
  rename(alpha = `alpha[1]`, Rmax = `Rmax[1]`)

curve(SR(SR_fun = "BH", alpha = sim1pop$pars_out$alpha,
         Rmax = sim1pop$pars_out$Rmax, S = x),
      from = 0, to = max(SRdat$S_true, SRdat$S_obs, quantile(SR$S, 0.975)), 
      ylim = range(0, SRdat$R_true, SRdat$R_obs, quantile(SR$R, 0.975), na.rm=TRUE)*1.02,
      xaxs = "i", yaxs = "i", lty = 3, lwd = 3, xlab = "Spawners", ylab = "Recruits", 
      las = 1, cex.axis = 1.2, cex.lab = 1.5)
for(i in sample(4000, 200))
  curve(SR(SR_fun = "BH", alpha = alphaRmax$alpha[i], Rmax = alphaRmax$Rmax[i], S = x),
        col = alpha("slategray4", 0.2), from = par("usr")[1], to = par("usr")[2], 
        add = TRUE)
segments(x0 = SRdat$S_true, x1 = SRdat$S_obs, y0 = SRdat$R_true, y1 = SRdat$R_obs,
         col = alpha("black", 0.3))
segments(x0 = SRdat$S_true, x1 = median(SRdat$S), y0 = SRdat$R_true, y1 = median(SRdat$R),
         col = alpha("black", 0.3))
points(R_true ~ S_true, data = SRdat, pch = 21, bg = "white", cex = 1.2)
points(R_obs ~ S_obs, data = SRdat, pch = 16, cex = 1.2)
points(median(R) ~ median(S), data = SRdat, pch = 16, cex = 1.2, col = "slategray4")
segments(x0 = quantile(SRdat$S, 0.025), x1 = quantile(SRdat$S, 0.975),
         y0 = median(SRdat$R), col = "slategray4")
segments(x0 = median(SRdat$S), y0 = quantile(SRdat$R, 0.025), 
         y1 = quantile(SRdat$R, 0.975), col = "slategray4")
legend("topleft", c("true","obs","states","fit"), cex = 1.2, bty = "n",
       pch = c(21,16,16,NA), pt.cex = 1.2, pt.bg = c("white",NA,NA,NA), 
       pt.lwd = 1, lty = c(3,NA,1,1), lwd = c(3,NA,1,1),
       col = c("black", "black", "slategray4", alpha("slategray4", 0.5)))


## @knitr singlepop_SR_ggplot
SR <- as_draws_rvars(as.matrix(fit1pop, c("S","R")))
alphaRmax <- as.data.frame(fit1pop, c("alpha", "Rmax")) %>% 
  rename(alpha = `alpha[1]`, Rmax = `Rmax[1]`)
RR <- run_recon(sim1pop$sim_dat)

cbind(RR, S_true = sim1pop$pars_out$S, R_true = sim1pop$pars_out$R,
      S = draws1pop$S, R = draws1pop$R) %>% 
  ggplot(aes(x = median(S), y = median(R))) +
  geom_function(fun = ~ SR(SR_fun = "BH", alpha = sim1pop$pars_out$alpha,
                           Rmax = sim1pop$pars_out$Rmax, S = .x),
                aes(lty = "true", col = "true"), lwd = 1) +
  lapply(sample(4000, 100), function(i) {
    geom_function(fun = ~ SR(SR_fun = "BH", alpha = alphaRmax$alpha[i],
                             Rmax = alphaRmax$Rmax[i], S = .x),
                  aes(lty = "fit", col = "fit"))
  }) +
  geom_segment(aes(x = S_true, xend = S_obs, y = R_true, yend = R_obs),
               col = "slategray4", alpha = 0.3) +
  geom_segment(aes(x = S_true, xend = median(S), y = R_true, yend = median(R)),
               col = "slategray4", alpha = 0.3) +
  geom_point(aes(x = S_true, y = R_true, pch = "true", col = "true"), 
             fill = "white", size = 2.5) +
  geom_point(aes(x = S_obs, y = R_obs, pch = "obs", col = "obs"), size = 2.5) +
  geom_point(aes(pch = "states", col = "states"), size = 2.5) +
  geom_segment(aes(x = quantile(S, 0.025), xend = quantile(S, 0.975),
                   y = median(R), yend = median(R), 
                   lty = "states", col = "states")) +
  geom_segment(aes(x = median(S), xend = median(S),
                   y = quantile(R, 0.025), yend = quantile(R, 0.975),
                   lty = "states", col = "states")) +
  scale_x_continuous(limits = c(0,NA), expand = c(0,1.05)) +
  scale_y_continuous(limits = c(0,NA), expand = c(0,1.05)) +
  scale_shape_manual(values = c(true = 21, obs = 16, states = 16, fit = NA)) +
  scale_linetype_manual(values = c(true = "dotted", obs = NA, 
                                   states = "solid", fit = "solid")) +
  scale_color_manual(values = c(true = "black", obs = "black", 
                                states = "slategray4",
                                fit = alpha("slategray4", 0.3))) +
  labs(x = "Spawners", y = "Recruits", shape = "", linetype = "", color = "") +
  theme(legend.position = c(0.1,0.93))
## @knitr


#-----------------------------------------------------
# Spawner time series plots
#-----------------------------------------------------

## @knitr singlepop_spawners_ggplot
draws1pop <- as_draws_rvars(fit1pop) %>% 
  mutate_variables(S_ppd = rvar_rng(rlnorm, N, log(S), tau))

sim1pop$sim_dat %>% cbind(S = draws1pop$S, S_ppd = draws1pop$S_ppd) %>% 
  ggplot(aes(x = year, y = S_obs)) + 
  geom_line(aes(y = sim1pop$pars_out$S, color = "true", lty = "true"), lwd = 1) +
  geom_ribbon(aes(ymin = quantile(S, 0.025), ymax = quantile(S, 0.975), fill = "states")) +
  geom_ribbon(aes(ymin = quantile(S_ppd, 0.025), ymax = quantile(S_ppd, 0.975), fill = "PPD")) +
  geom_line(aes(y = median(draws1pop$S), lty = "states", col = "states"), lwd = 1) +
  geom_point(aes(col = "obs", pch = "obs"), size = 3) + 
  scale_y_continuous(limits = c(0,NA), expand = c(0,0)) +
  scale_shape_manual(values = c(true = NA, obs = 16, states = NA, PPD = NA)) +
  scale_color_manual(values = c(true = "black", obs = "black",
                                states = "slategray4", PPD = "white")) +
  scale_fill_manual(values = c(true = "white", obs = "white", 
                               states = alpha("slategray4", 0.3),
                               PPD = alpha("slategray4", 0.2))) +
  scale_linetype_manual(values = c(true = "dotted", obs = NA,  
                                   states = "solid", PPD = NA)) +
  labs(x = "Year", y = "Spawners", shape = "", color = "", fill = "", linetype = "") +
  theme(legend.position = c(0.9,0.93))


## @knitr singlepop_spawners_base
draws1pop <- as_draws_rvars(fit1pop) %>% 
  mutate_variables(S_ppd = rvar_rng(rlnorm, N, log(S), tau))
ppd1pop <- sim1pop$sim_dat %>% cbind(S = draws1pop$S, S_ppd = draws1pop$S_ppd) 

plot(ppd1pop$year, sim1pop$pars_out$S, type = "l", lty = 3, lwd = 3,
     ylim = c(0, max(quantile(ppd1pop$S_ppd, 0.975))), yaxs = "i",
     xlab = "Year", ylab = "Spawners", cex.axis = 1.2, cex.lab = 1.5, las = 1)
polygon(c(ppd1pop$year, rev(ppd1pop$year)),
        c(quantile(ppd1pop$S, 0.025), rev(quantile(ppd1pop$S, 0.975))),
        col = alpha("slategray4", 0.3), border = NA)
polygon(c(ppd1pop$year, rev(ppd1pop$year)),
        c(quantile(ppd1pop$S_ppd, 0.025), rev(quantile(ppd1pop$S_ppd, 0.975))),
        col = alpha("slategray4", 0.2), border = NA)
lines(median(S) ~ year, data = ppd1pop, col = "slategray4", lwd = 3)
points(S_obs ~ year, data = ppd1pop, pch = 16, cex = 1.5)
legend("topright", c("true","obs","states","PPD"), cex = 1.2, y.intersp = 1.2,
       pch = c(NA,16,NA,NA), pt.cex = 1.5, lty = c(3,NA,1,1), lwd = c(3,NA,15,15), 
       col = c("black", "black", alpha("slategray4", c(0.5,0.2))), bty = "n")
## @knitr


#===========================================================================
# ENVIRONMENTAL COVARIATES
#===========================================================================

#------------------------------
# Data structure
# - add covariates
#------------------------------

## @knitr singlepop_covariate_data_struct
df1pop$X1 <- rnorm(N,0,1)
df1pop$X2 <- rnorm(N,0,1)
## @knitr

#------------------------------
# True parameter values 
# - add regression coefs
#------------------------------

## @knitr singlepop_covariate_pars
pars1pop$beta_Rmax <- 0.5
pars1pop$beta_R <- -0.5
## @knitr

#------------------------------
# Simulate data 
#------------------------------

## @knitr singlepop_covariate_data
simX1pop <- simIPM(life_cycle = "SS", SR_fun = "BH", 
                   N_age = N_age, max_age = max_age,
                   pars = pars1pop, par_models = list(Rmax ~ X1, R ~ X2), 
                   fish_data = df1pop)
                   
format(head(simX1pop$sim_dat, 10), digits = 2)
## @knitr

#-----------------------------------------------------
# Fit IPM
#-----------------------------------------------------

## @knitr singlepop_covariate_fit
fitX1pop <- salmonIPM(life_cycle = "SS", pool_pops = FALSE, SR_fun = "BH", 
                      par_models = list(Rmax ~ X1, R ~ X2),
                      fish_data = simX1pop$sim_dat, 
                      chains = 4, iter = 2000, warmup = 1000, 
                      seed = 123)

print(fitX1pop, pars = c("beta_Rmax","beta_R"), prob = c(0.025, 0.5, 0.975))
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
parsNpop <- list(mu_alpha = 2, sigma_alpha = 0.3, mu_Rmax = 3, sigma_Rmax = 0.3, 
                 rho_alphaRmax = 0.5, rho_R = 0.6, sigma_year_R = 0.3, sigma_R = 0.1,
                 mu_p = c(0.05, 0.55, 0.4), sigma_pop_p = c(0.2, 0.3), 
                 R_pop_p = matrix(c(1, 0.3, 0.3, 1), 2, 2), 
                 sigma_p = c(0.3, 0.5), R_p = matrix(c(1, 0.5, 0.5, 1), 2, 2), 
                 tau = 0.3, S_init_K = 0.3)
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
                     p_HOS = c(rep(0, N/2), runif(N/2, 0, 0.5)), 
                     B_rate = c(rep(0, N/2), rbeta(N/2, 1, 19)),
                     F_rate = rbeta(N, 3, 2),
                     n_age_obs = replace(runif(N, 10, 100), sample(N, N/10), 0),
                     n_HW_obs = c(rep(0, N/2), runif(N/2, 10, 100)))
## @knitr

#------------------------------
# Simulate data 
# - some spawner obs missing
#------------------------------

## @knitr multipop_data
simNpop <- simIPM(life_cycle = "SS", SR_fun = "BH", 
                  N_age = N_age, max_age = max_age,
                  pars = parsNpop, fish_data = dfNpop)
simNpop$sim_dat$S_obs[sample(N, N/10)] <- NA
names(simNpop$pars_out)
simNpop$pars_out[c("alpha","Rmax")]
format(head(simNpop$sim_dat, 10), digits = 2)
## @knitr

#-----------------------------------------------------
# IPM with no pooling across populations
#-----------------------------------------------------

## @knitr fit_np
fitNnp <- salmonIPM(life_cycle = "SS", pool_pops = FALSE, SR_fun = "BH", 
                    fish_data = simNpop$sim_dat,
                    chains = 4, iter = 2000, warmup = 1000, 
                    seed = 321)

print(fitNnp, pars = stan_pars("IPM_SS_np", "hyper"), prob = c(0.025, 0.5, 0.975))
## @knitr

#-----------------------------------------------------
# IPM with partial population pooling
#-----------------------------------------------------

## @knitr fit_pp
fitNpp <- salmonIPM(life_cycle = "SS", pool_pops = TRUE, SR_fun = "BH",
                    fish_data = simNpop$sim_dat,
                    chains = 4, iter = 2000, warmup = 1000, 
                    seed = 321)

print(fitNpp, pars = stan_pars("IPM_SS_pp", "hyper"), prob = c(0.025, 0.5, 0.975))
## @knitr

#----------------------------------------------------------
# Plot pop-level S-R parameter posteriors and true values
# - No population pooling vs. partial population pooling
#----------------------------------------------------------

## @knitr multipop_np_vs_pp
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
## @knitr

#----------------------------------------------------------
# Plot hyperparameter posteriors, priors, and true values
# for partial population pooling model
#----------------------------------------------------------

## @knitr multipop_posteriors
par_names <- c("mu_alpha","mu_Rmax","sigma_year_R","rho_R","sigma_R",
               "mu_p","sigma_pop_p","rho_pop_p","sigma_p","rho_p","tau")

# extract and transform draws using posterior package
postNpop <- as_draws_rvars(fitNpp) %>%
  mutate_variables(mu_p = as.vector(mu_p), sigma_pop_p = as.vector(sigma_pop_p),
                   rho_pop_p = as.vector(R_pop_p[2,1]), sigma_p = as.vector(sigma_p),
                   rho_p = as.vector(R_p[2,1])) %>%
  .[par_names]

postdfNpop <- postNpop %>% as_draws_df() %>% as.data.frame() %>% 
  pivot_longer(cols = !starts_with("."), names_to = "param_indx", values_to = "value") %>% 
  mutate(param_indx = factor(param_indx, levels = unique(param_indx)))

# specify priors using distributional package
log_S_obs <- na.omit(log(simNpop$sim_dat$S_obs))
priorNpop <- c(mu_alpha = dist_normal(2,2),
               mu_Rmax = dist_normal(max(log_S_obs), sd(log_S_obs)),
               sigma_year_R = dist_normal(0,3),
               rho_R = dist_wrap("gnorm", 0, 0.85, 20),
               sigma_R = dist_normal(0,3),
               mu_p = dist_beta(1, N_age - 1),
               sigma_pop_p = dist_normal(0,2),
               rho_pop_p = 2*dist_beta((N_age - 1)/2, (N_age - 1)/2) - 1, # LKJ
               sigma_p = dist_normal(0,2),
               rho_p = 2*dist_beta((N_age - 1)/2, (N_age - 1)/2) - 1, # LKJ
               tau = dist_normal(0,1))

priordfNpop <- data.frame(param_indx = unique(postdfNpop$param_indx)) %>% 
  mutate(param = gsub("\\[.\\]", "", param_indx), prior = priorNpop[param])

# true parameter values
trueNpop <- simNpop$pars_out %>%
  c(rho_pop_p = .$R_pop_p[2,1], rho_p = .$R_p[2,1]) %>%
  .[par_names] %>% unlist() %>% setNames(gsub("(\\d)", "[\\1]", names(.))) %>% 
  data.frame(param_indx = factor(names(.), levels = names(.)), value = .)

# plot
scales <- postNpop %>% as_draws_df() %>% as.data.frame() %>% select(!starts_with(".")) %>% 
  lapply(FUN = function(x) scale_x_continuous(limits = range(x)))

postdfNpop %>%   
  ggplot(aes(x = value)) + 
  geom_histogram(aes(y = after_stat(density)), color = "white", fill = "slategray4", alpha = 0.5) +
  stat_slab(data = priordfNpop, aes(xdist = prior), inherit.aes = FALSE,
            normalize = "none", col = "black", lwd = 0.8, fill = NA) +
  geom_vline(data = trueNpop, aes(xintercept = value), lwd = 0.8, lty = 2) +
  facet_wrap(~ param_indx, scales = "free", strip.position = "bottom") + 
  facetted_pos_scales(x = scales) + theme_classic() + 
  theme(strip.placement = "outside", strip.background = element_blank(), 
        strip.text = element_text(size = 11, margin = margin(b = 3, t = 0)), 
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text = element_text(size = 10), axis.text.y = element_blank()) +
  labs(x = "", y = "")
## @knitr
