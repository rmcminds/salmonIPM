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
library(rgl)             # 3-D graphics
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
N <- 50
N_age <- 3
max_age <- 5
## @knitr

#------------------------------
# True parameter values 
#------------------------------

## @knitr pars
pars <- list(mu_alpha_W = 1.5, mu_alpha_H = 0.5, sigma_alpha = 0, 
             mu_Rmax_W = 7, mu_Rmax_H = 6, sigma_Rmax = 0, 
             rho_alphaRmax = 0, rho_WH = 0, rho_R = 0.5, sigma_year_R = 0.3, sigma_R = 0,
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

## @knitr data_struct
df <- data.frame(pop = 1, year = 1:N + 2020 - N, A = 1, 
                 p_HOS = rbeta(N, 2, 2), 
                 B_rate = rbeta(N, 1, 9), 
                 F_rate = rbeta(N, 3, 2),
                 n_age_obs = runif(N, 10, 100),
                 n_HW_obs = runif(N, 10, 100))

#------------------------------
# Simulate data 
#------------------------------

## @knitr data
sim <- simIPM(life_cycle = "SS", SR_fun = "Ricker", pars = pars, 
              fish_data = df, N_age = N_age, max_age = max_age)
names(sim$pars_out)
sim$pars_out[c("alpha_W","alpha_H","Rmax_W","Rmax_H")]
format(head(sim$sim_dat, 10), digits = 2)
## @knitr

#-----------------------------------------------------
# Fit IPM
#-----------------------------------------------------

## @knitr fit
fit <- salmonIPM(life_cycle = "SS", pool_pops = FALSE, 
                 SR_fun = "Ricker", RRS = c("alpha","Rmax"),
                 fish_data = sim$sim_dat, 
                 chains = 4, iter = 2000, warmup = 1000, 
                 seed = 123)

print(fit, pars = stan_pars("IPM_SS_np", "hyper"), prob = c(0.025, 0.5, 0.975))
## @knitr

#-----------------------------------------------------
# Plot posteriors, priors, and true values
#-----------------------------------------------------

## @knitr posteriors
par_names <- c("log(alpha_W)","log(alpha_H)","log(Rmax_W)","log(Rmax_H)")

# specify priors using distributional package
prior_Rmax <- stan_data("IPM_SS_np", fish_data = sim$sim_dat)$prior_Rmax
prior <- list(`log(alpha_W)` = dist_normal(2,2), `log(alpha_H)` = dist_normal(2,2),
              `log(Rmax_W)` = dist_normal(prior_Rmax[1], prior_Rmax[2]),
              `log(Rmax_H)` = dist_normal(prior_Rmax[1], prior_Rmax[2]))

# true parameter values
true <- sim$pars_out %>% 
  c(`log(alpha_W)` = log(.$alpha_W), `log(alpha_H)` = log(.$alpha_H), 
    `log(Rmax_W)` = log(.$Rmax_W), `log(Rmax_H)` = log(.$Rmax_H)) %>% 
  .[par_names] %>% unlist()

# extract and transform draws using posterior package
post <- as_draws_rvars(fit) %>% 
  mutate_variables(`log(alpha_W)` = log(alpha_W), `log(alpha_H)` = log(alpha_H), 
                   `log(Rmax_W)` = log(Rmax_W), `log(Rmax_H)` = log(Rmax_H)) %>% 
  .[par_names] %>% as_draws_matrix()

# plot
par(mfcol = c(2,2), mar = c(5,1,1,1))
for(j in names(true)) {
  hist(post[,j], 20, prob = TRUE, col = alpha("slategray4", 0.5), border = "white",
       xlab = j, ylab = "", yaxt = "n", main = "", cex.axis = 1.2, cex.lab = 1.4)
  curve(density(prior[[j]], at = x)[[1]], lwd = 0.5, add = TRUE)
  abline(v = true[[j]], lwd = 2, lty = 3)
}
legend("right", c("true","prior","posterior"), cex = 1.4, text.col = "white", 
       fill = c(NA, NA, alpha("slategray4", 0.5)), border = NA,
       inset = c(-1,0), xpd = NA, bty = "n")
legend("right", c("true","prior","posterior"), cex = 1.4, 
       lty = c(3,1,NA), lwd = c(2,1,NA), col = c(rep("black",2), "slategray4"),
       inset = c(-1,0), xpd = NA, bty = "n")
## @knitr

#-----------------------------------------------------
# Plot pure W or H S-R curves and obs, states
#-----------------------------------------------------

## @knitr SR_ggplot
# states and observations including reconstructed recruits
RR <- run_recon(sim$sim_dat)
states <- as.matrix(fit, c("S","p_HOS","R")) %>% as_draws_rvars()
states_obs <- RR %>% select(pop, year, A, S_obs, R_obs) %>% 
  data.frame(states) %>% mutate(p_HOS = median(p_HOS))

# W and H spawner abundances at which to evaluate S-R function
S_grid <- states_obs %>% 
  reframe(S = rep(seq(0, max(quantile(S, 0.9), S_obs), length = 100), 2)) %>% 
  mutate(p_HOS = rep(c(0,1), each = 100), S_W = S*(1 - p_HOS), S_H = S*p_HOS)

# posteriors of S-R fit
SR_fit <- as.matrix(fit, c("alpha_W","alpha_H","Rmax_W","Rmax_H")) %>% 
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

## @knitr SR3D_rgl
# states and observations
states <- as.matrix(fit, c("S","p_HOS","R")) %>% 
  as_draws_rvars() %>% mutate_variables(S_W = S*(1 - p_HOS), S_H = S*p_HOS)

# S-R parameters
SR_pars <- as.matrix(fit, c("alpha_W","alpha_H","Rmax_W","Rmax_H")) %>% 
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




