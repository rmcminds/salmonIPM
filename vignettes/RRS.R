# graphics device for this script
if(.Platform$OS.type == "windows") options(device = "windows")

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
library(rgl)             # 3-D graphics
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
set.seed(123)
N <- 50
N_age <- 3
max_age <- 5
## @knitr

#------------------------------
# True parameter values 
#------------------------------

## @knitr singlepop_pars
pars1pop <- list(mu_alpha_W = 1.5, mu_alpha_H = 0.5, sigma_alpha = 0, 
                 mu_Rmax_W = 7, mu_Rmax_H = 6, sigma_Rmax = 0, 
                 rho_alphaRmax = 0, rho_WH = 0, 
                 rho_R = 0.5, sigma_year_R = 0.3, sigma_R = 0,
                 mu_p = c(0.05, 0.55, 0.4), sigma_pop_p = rep(0,2), R_pop_p = diag(2), 
                 sigma_p = c(0.2, 0.2), R_p = matrix(c(1, 0.5, 0.5, 1), 2, 2), 
                 tau = 0.2, S_init_K = 0.2)
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
                     F_rate = rbeta(N, 3, 2),
                     n_age_obs = runif(N, 10, 100),
                     n_HW_obs = runif(N, 10, 100))

#------------------------------
# Simulate data 
#------------------------------

## @knitr singlepop_data
sim1pop <- simIPM(life_cycle = "SS", SR_fun = "Ricker", pars = pars1pop, 
                  fish_data = df1pop, N_age = N_age, max_age = max_age)
names(sim1pop$pars_out)
sim1pop$pars_out[c("alpha_W","alpha_H","Rmax_W","Rmax_H")]
format(head(sim1pop$sim_dat, 10), digits = 2)
## @knitr

#-----------------------------------------------------
# Fit IPM
#-----------------------------------------------------

## @knitr singlepop_fit
fit1pop <- salmonIPM(life_cycle = "SS", pool_pops = FALSE, 
                     SR_fun = "Ricker", RRS = c("alpha","Rmax"),
                     fish_data = sim1pop$sim_dat, 
                     seed = 123)

print(fit1pop)
## @knitr

#-----------------------------------------------------
# Plot posteriors, priors, and true values
#-----------------------------------------------------

## @knitr singlepop_posteriors
par_names <- c("log(alpha_W)","log(alpha_H)","log(Rmax_W)","log(Rmax_H)")

# extract and transform draws using posterior package
post1pop <- as_draws_rvars(fit1pop) %>% 
  mutate_variables(`log(alpha_W)` = log(alpha_W), `log(alpha_H)` = log(alpha_H), 
                   `log(Rmax_W)` = log(Rmax_W), `log(Rmax_H)` = log(Rmax_H)) %>% 
  .[par_names]

postdf1pop <- post1pop %>% as_draws_df() %>% as.data.frame() %>% 
  pivot_longer(cols = !starts_with("."), names_to = "param_indx", values_to = "value") %>% 
  mutate(param_indx = factor(param_indx, levels = unique(param_indx)))

# specify priors using distributional package
prior_Rmax <- stan_data("IPM_SS_np", fish_data = sim1pop$sim_dat)$prior_Rmax
prior1pop <- list(`log(alpha_W)` = dist_normal(2,2), `log(alpha_H)` = dist_normal(2,2),
                  `log(Rmax_W)` = dist_normal(prior_Rmax[1], prior_Rmax[2]),
                  `log(Rmax_H)` = dist_normal(prior_Rmax[1], prior_Rmax[2]))

priordf1pop <- data.frame(param_indx = unique(postdf1pop$param_indx)) %>% 
  mutate(param = gsub("\\[.\\]", "", param_indx), prior = prior1pop[param])

# true parameter values
true1pop <- sim1pop$pars_out %>% 
  c(`log(alpha_W)` = log(.$alpha_W), `log(alpha_H)` = log(.$alpha_H), 
    `log(Rmax_W)` = log(.$Rmax_W), `log(Rmax_H)` = log(.$Rmax_H)) %>% 
  .[par_names] %>% unlist() %>% 
  data.frame(param_indx = gsub("(\\d)", "[\\1]", names(.)), value = ., row.names = NULL) %>% 
  mutate(param_indx = factor(param_indx, levels = levels(postdf1pop$param_indx)))

# plot
scales <- post1pop %>% as_draws_df() %>% as.data.frame() %>% select(!starts_with(".")) %>% 
  lapply(function(x) scale_x_continuous(limits = range(x)))

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
## @knitr

#-----------------------------------------------------
# Plot pure W or H S-R curves and obs, states
#-----------------------------------------------------

## @knitr singlepop_SR_ggplot
# states and observations including reconstructed recruits
RR <- run_recon(sim1pop$sim_dat)
states <- as.matrix(fit1pop, c("S","p_HOS","R")) %>% as_draws_rvars()
states_obs <- RR %>% select(pop, year, A, S_obs, R_obs) %>% 
  data.frame(states) %>% mutate(p_HOS = median(p_HOS))

# W and H spawner abundances at which to evaluate S-R function
S_grid <- states_obs %>% 
  reframe(S = rep(seq(0, max(quantile(S, 0.9), S_obs), length = 100), 2)) %>% 
  mutate(p_HOS = rep(c(0,1), each = 100), S_W = S*(1 - p_HOS), S_H = S*p_HOS)

# posteriors of S-R fit
SR_fit <- as.matrix(fit1pop, c("alpha_W","alpha_H","Rmax_W","Rmax_H")) %>% 
  as_draws_rvars() %>% 
  mutate_variables(S_W = as_rvar(S_grid$S_W), S_H = as_rvar(S_grid$S_H),
                   alpha_W = rep(alpha_W, each = 200), 
                   alpha_H = rep(alpha_H, each = 200),
                   Rmax_W = rep(Rmax_W, each = 200),
                   Rmax_H = rep(Rmax_H, each = 200),
                   R = SR(SR_fun = "Ricker",
                          alpha_W = alpha_W, alpha_H = alpha_H,
                          Rmax_W = Rmax_W, Rmax_H = Rmax_H, 
                          S_W = S_W, S_H = S_H))

SR_fit <- data.frame(S = S_grid$S, p_HOS = S_grid$p_HOS, R = SR_fit$R)

# plot
SR_fit %>% ggplot(aes(x = S, y = median(R), color = p_HOS, fill = p_HOS, 
                      group = factor(p_HOS))) +
  geom_ribbon(aes(ymin = t(quantile(R, 0.25)), ymax = t(quantile(R, 0.75))), 
              color = NA, alpha = 0.3) +
  geom_line(lwd = 1) +
  geom_segment(aes(x = t(quantile(S, 0.25)), xend = t(quantile(S, 0.75)),
                   y = median(R), yend = median(R)),
               data = states_obs, alpha = 0.8) +
  geom_segment(aes(x = median(S), xend = median(S),
                   y = t(quantile(R, 0.25)), yend = t(quantile(R, 0.75))),
               data = states_obs, alpha = 0.8) +
  geom_segment(aes(x = S_obs, xend = median(S), y = R_obs, yend = median(R)),
               data = states_obs, alpha = 0.5) +
  geom_point(aes(x = S_obs, y = R_obs), data = states_obs,
             pch = 16, size = 2.5, alpha = 0.8) +
  geom_point(aes(x = median(S), y = median(R)), data = states_obs,
             pch = 21, size = 2.5, fill = "white") +
  labs(x = "Spawners", y = "Recruits", 
       color = bquote(italic(p)[HOS]), fill = bquote(italic(p)[HOS])) + 
  scale_x_continuous(expand = c(0.01,0)) + scale_y_continuous(expand = c(0.01,0)) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_fill_gradient(low = "blue", high = "red") +
  theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11),
        panel.grid = element_blank(), strip.background = element_rect(fill = NA),
        strip.text = element_text(margin = margin(b = 2, t = 1)))
## @knitr

#---------------------------------------------------------------------
# 3D surface plot of R vs (S_W, S_H) at median S-R parameter values
#---------------------------------------------------------------------

## @knitr singlepop_SR3D_rgl
# states and observations
states <- as.matrix(fit1pop, c("S","p_HOS","R")) %>% 
  as_draws_rvars() %>% mutate_variables(S_W = S*(1 - p_HOS), S_H = S*p_HOS)

# S-R parameters
SR_pars <- as.matrix(fit1pop, c("alpha_W","alpha_H","Rmax_W","Rmax_H")) %>% 
  as_draws_rvars()

# W and H spawner abundances at which to evaluate S-R function
S_grid <- states %>% data.frame() %>% 
  reframe(S_W = seq(0, max(median(S_W), median(S_H)), length = 100),
          S_H = seq(0, max(median(S_W), median(S_H)), length = 100))

# recruits
R_fit <- outer(S_grid$S_W, S_grid$S_H, function(x, y) 
  SR(SR_fun = "Ricker",
     alpha_W = median(SR_pars$alpha_W), alpha_H = median(SR_pars$alpha_H), 
     Rmax_W = median(SR_pars$Rmax_W), Rmax_H = median(SR_pars$Rmax_H),
     S_W = x, S_H = y))

# map p_HOS onto color gradient
cfun <- function(S_W, S_H)
  rgb(colorRamp(c("blue","red"))(S_H/(S_W + S_H)), maxColorValue = 255)
csurf <- outer(c(0.01, S_grid$S_W[-1]), c(0.01, S_grid$S_H[-1]), cfun)
cpts <- cfun(median(states$S_W), median(states$S_H))

# plot
par3d(windowRect = c(200,100,1000,900))
persp3d(S_grid$S_W, S_grid$S_H, R_fit, col = csurf, alpha = 0.7,
        xlab = bquote(S[W]), ylab = bquote(S[H]), zlab = "R")
spheres3d(median(states$S_W), median(states$S_H), median(states$R), 
          radius = 40, col = cpts)
view3d(userMatrix = matrix(c(-0.82,0.37,-0.42,0,-0.56,-0.42,0.72,0,
                             0.09,0.83,0.55,0,0,0,0,1), 4, 4))
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
parsNpop <- list(mu_alpha_W = 1.5, mu_alpha_H = 0.5, sigma_alpha = 0.3,
                 mu_Rmax_W = 7, mu_Rmax_H = 5, sigma_Rmax = 0.3,
                 rho_alphaRmax = 0.5, rho_WH = 0.3,
                 rho_R = 0.3, sigma_year_R = 0.1, sigma_R = 0.1,
                 mu_p = c(0.05, 0.55, 0.4), sigma_pop_p = c(0.2, 0.3),
                 R_pop_p = matrix(c(1, 0.3, 0.3, 1), 2, 2),
                 sigma_p = c(0.1, 0.2), R_p = matrix(c(1, 0.5, 0.5, 1), 2, 2),
                 tau = 0.2, S_init_K = 0.3)
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
                     p_HOS = runif(N, 0.2, 0.8), 
                     B_rate = rbeta(N, 1, 9), 
                     F_rate = rbeta(N, 6, 4),
                     n_age_obs = runif(N, 10, 100),
                     n_HW_obs = runif(N, 10, 100))
## @knitr

#------------------------------
# Simulate data 
#------------------------------

## @knitr multipop_data
simNpop <- simIPM(life_cycle = "SS", SR_fun = "BH", pars = parsNpop, 
                  fish_data = dfNpop, N_age = N_age, max_age = max_age)
## @knitr

#-----------------------------------------------------
# IPM with no pooling across populations
#-----------------------------------------------------

## @knitr fit_np
fitNnp <- salmonIPM(life_cycle = "SS", pool_pops = FALSE, 
                    SR_fun = "BH", RRS = c("alpha", "Rmax"),
                    fish_data = simNpop$sim_dat,
                    control = list(max_treedepth = 13),
                    seed = 321)

print(fitNnp, pars = c("alpha_W","alpha_H","delta_alpha","Rmax_W","Rmax_H","delta_Rmax"))
## @knitr

#-----------------------------------------------------
# IPM with partial population pooling
#-----------------------------------------------------

## @knitr fit_pp
fitNpp <- salmonIPM(life_cycle = "SS", pool_pops = TRUE, 
                    SR_fun = "BH", RRS = c("mu_alpha", "mu_Rmax"),
                    fish_data = simNpop$sim_dat,
                    seed = 321)

print(fitNpp, 
      pars = c("mu_alpha_W","mu_alpha_H","delta_mu_alpha",
               "mu_Rmax_W","mu_Rmax_H","delta_mu_Rmax"))
## @knitr

#----------------------------------------------------------
# Plot pop-level H/W discounts posteriors and true values
# - No population pooling vs. partial population pooling
#----------------------------------------------------------

## @knitr multipop_np_vs_pp
par(mfcol = c(2,1), mar = c(1,4,0,1), oma = c(3,0,0,0))

# intrinsic productivity
vioplot(as.matrix(fitNnp, "delta_alpha"), col = "slategray4", 
        border = NA, drawRect = FALSE, side = "left", las = 1, 
        xlab = "", names = NA, ylab = "", cex.axis = 1.2, cex.lab = 1.5)
vioplot(as.matrix(fitNpp, "delta_alpha"), col = "salmon", 
        border = NA, drawRect = FALSE, side = "right", add = TRUE)
points(1:N_pop, log(simNpop$pars_out$alpha_H) - log(simNpop$pars_out$alpha_W), 
       pch = 16, cex = 1.5)
mtext("delta_alpha", side = 2, line = 2.5, cex = 1.5)
legend("top", c("true","np","pp"), cex = 1.2, bty = "n", horiz = TRUE, 
       pch = c(16,NA,NA), pt.cex = 1.5, fill = c(NA,"slategray4","salmon"), border = NA)

# maximum recruitment
vioplot(as.matrix(fitNnp, "delta_Rmax"), col = "slategray4", 
        border = NA, drawRect = FALSE, side = "left", las = 1, 
        xlab = "", names = LETTERS[1:N_pop], ylab = "", cex.axis = 1.2, cex.lab = 1.5)
vioplot(as.matrix(fitNpp, "delta_Rmax"), col = "salmon", 
        border = NA, drawRect = FALSE, side = "right", add = TRUE)
points(1:N_pop, log(simNpop$pars_out$Rmax_W) - log(simNpop$pars_out$Rmax_W), 
       pch = 16, cex = 1.5)
mtext(c("Population","delta_Rmax"), side = 1:2, line = 2.5, cex = 1.5)
## @knitr
