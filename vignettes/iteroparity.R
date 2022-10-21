#===========================================================================
# SETUP
#===========================================================================

# graphics device for this script
options(device = ifelse(.Platform$OS.type == "windows", "windows", "quartz"))

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
N_age <- 5
max_age <- 7
## @knitr

#------------------------------
# True parameter values 
#------------------------------

## @knitr singlepop_pars
pars1pop <- list(mu_alpha = 2, sigma_alpha = 0, mu_Rmax = 7, sigma_Rmax = 0, 
                 rho_alphaRmax = 0, rho_R = 0.6, sigma_year_R = 0.2, sigma_R = 0,
                 mu_p = c(0.1, 0.3, 0.4, 0.1, 0.1), sigma_pop_p = rep(0,4), R_pop_p = diag(4), 
                 sigma_p = c(0.1, 0.2, 0.2, 0.1), R_p = diag(4), 
                 mu_SS = 0.2, rho_SS = 0.6, sigma_year_SS = 0.2, sigma_SS = 0,
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
sim1pop <- simIPM(life_cycle = "SSsthd", SR_fun = "BH", pars = pars1pop, 
                  fish_data = df1pop, N_age = N_age, max_age = max_age)
names(sim1pop$pars_out)
sim1pop$pars_out[c("alpha","Rmax")]
format(head(sim1pop$sim_dat, 10), digits = 2)
## @knitr

#-----------------------------------------------------
# Fit IPM
#-----------------------------------------------------

## @knitr singlepop_fit
fit1pop <- salmonIPM(life_cycle = "SSsthd", pool_pops = FALSE, SR_fun = "BH", 
                     fish_data = sim1pop$sim_dat, 
                     chains = 4, iter = 2000, warmup = 1000, 
                     control = list(adapt_delta = 0.95),
                     seed = 123)

print(fit1pop, pars = c("p","R_p","s_SS","p_HOS","S","R","q_MR","LL"), 
      include = FALSE, prob = c(c(0.025, 0.5, 0.975)))
## @knitr

# #-----------------------------------------------------
# # Plot posteriors, priors, and true values
# #-----------------------------------------------------
# 
# ## @knitr singlepop_posteriors
# par_names <- c("log(alpha)","log(Rmax)","rho_R","sigma_R","mu_p","sigma_p","rho_p","tau")
# 
# # specify priors using distributional package
# log_S_obs <- log(sim1pop$sim_dat$S_obs)
# prior1pop <- c(`log(alpha)` = dist_normal(2,2),
#                `log(Rmax)` = dist_normal(max(log_S_obs), sd(log_S_obs)),
#                rho_R = dist_wrap("pexp", 0, 0.85, 20),
#                sigma_R = dist_normal(0,5),
#                mu_p = dist_beta(1,2),
#                sigma_p = dist_normal(0,5),
#                rho_p = dist_uniform(-1,1),
#                tau = dist_wrap("pexp", 1, 0.85, 30))
# 
# # true parameter values
# true1pop <- sim1pop$pars_out %>% replace("sigma_R", .$sigma_year_R) %>% 
#   c(`log(alpha)` = log(.$alpha), `log(Rmax)` = log(.$Rmax), rho_p = .$R_p[2,1]) %>% 
#   .[par_names] %>% unlist() %>% setNames(gsub("(\\d)", "[\\1]", names(.)))
# 
# # extract and transform draws using posterior package
# post1pop <- as_draws_rvars(fit1pop) %>% 
#   mutate_variables(`log(alpha)` = log(alpha), `log(Rmax)` = log(Rmax),
#                    mu_p = as.vector(mu_p), sigma_p = as.vector(sigma_p), 
#                    rho_p = R_p[1,2,1]) %>% 
#   as_draws_matrix(.[par_names])
# 
# # plot
# par(mfrow = c(3,4), mar = c(5,1,0,1))
# for(j in names(true1pop)) {
#   hist(post1pop[,j], 20, prob = TRUE, col = alpha("slategray4", 0.5), border = "white",
#        xlab = j, ylab = "", yaxt = "n", main = "", cex.axis = 1.2, cex.lab = 1.5)
#   curve(density(prior1pop[gsub("\\[.\\]", "", j)], at = x)[[1]], lwd = 0.5, add = TRUE)
#   abline(v = true1pop[[j]], lwd = 2, lty = 3)
# }
# legend("right", c("true","prior","posterior"), cex = 1.5, text.col = "white", 
#        fill = c(NA, NA, alpha("slategray4", 0.5)), border = NA,
#        inset = c(-1,0), xpd = NA, bty = "n")
# legend("right", c("true","prior","posterior"), cex = 1.5, 
#        lty = c(3,1,NA), lwd = c(2,1,NA), col = c(rep("black",2), "slategray4"),
#        inset = c(-1,0), xpd = NA, bty = "n")
# ## @knitr
# 
# #-----------------------------------------------------
# # Plot true S-R curve, obs, states and fitted draws
# #-----------------------------------------------------
# 
# ## @knitr singlepop_SR_base
# SR <- as_draws_rvars(as.matrix(fit1pop, c("S","R")))
# alphaRmax <- as.data.frame(fit1pop, c("alpha", "Rmax")) %>% 
#   rename(alpha = `alpha[1]`, Rmax = `Rmax[1]`)
# RR <- run_recon(sim1pop$sim_dat)
# SRdat <- cbind(RR, S_true = sim1pop$pars_out$S, R_true = sim1pop$pars_out$R,
#                S = SR$S, R = SR$R)
# 
# curve(SR(SR_fun = "BH", alpha = sim1pop$pars_out$alpha,
#          Rmax = sim1pop$pars_out$Rmax, S = x),
#       from = 0, to = max(SRdat$S_true, SRdat$S_obs, quantile(SR$S, 0.975)), 
#       ylim = range(0, SRdat$R_true, SRdat$R_obs, quantile(SR$R, 0.975), na.rm=TRUE)*1.02,
#       xaxs = "i", yaxs = "i", lty = 3, lwd = 3, xlab = "Spawners", ylab = "Recruits", 
#       las = 1, cex.axis = 1.2, cex.lab = 1.5)
# for(i in sample(4000, 100))
#   curve(SR(SR_fun = "BH", alpha = alphaRmax$alpha[i], Rmax = alphaRmax$Rmax[i], S = x),
#         col = alpha("slategray4", 0.3), from = par("usr")[1], to = par("usr")[2], 
#         add = TRUE)
# segments(x0 = SRdat$S_true, x1 = SRdat$S_obs, y0 = SRdat$R_true, y1 = SRdat$R_obs,
#          col = alpha("black", 0.3))
# segments(x0 = SRdat$S_true, x1 = median(SRdat$S), y0 = SRdat$R_true, y1 = median(SRdat$R),
#          col = alpha("black", 0.3))
# points(R_true ~ S_true, data = SRdat, pch = 21, bg = "white", cex = 1.2)
# points(R_obs ~ S_obs, data = SRdat, pch = 16, cex = 1.2)
# points(median(R) ~ median(S), data = SRdat, pch = 16, cex = 1.2, col = "slategray4")
# segments(x0 = quantile(SRdat$S, 0.025), x1 = quantile(SRdat$S, 0.975),
#          y0 = median(SRdat$R), col = "slategray4")
# segments(x0 = median(SRdat$S), y0 = quantile(SRdat$R, 0.025), 
#          y1 = quantile(SRdat$R, 0.975), col = "slategray4")
# legend("topleft", c("true","obs","states","fit"), cex = 1.2, bty = "n",
#        pch = c(21,16,16,NA), pt.cex = 1.2, pt.bg = c("white",NA,NA,NA), 
#        pt.lwd = 1, lty = c(3,NA,1,1), lwd = c(3,NA,1,1),
#        col = c("black", "black", "slategray4", alpha("slategray4", 0.3)))
# 
# 
# ## @knitr singlepop_SR_ggplot
# SR <- as_draws_rvars(as.matrix(fit1pop, c("S","R")))
# alphaRmax <- as.data.frame(fit1pop, c("alpha", "Rmax")) %>% 
#   rename(alpha = `alpha[1]`, Rmax = `Rmax[1]`)
# RR <- run_recon(sim1pop$sim_dat)
# 
# cbind(RR, S_true = sim1pop$pars_out$S, R_true = sim1pop$pars_out$R,
#       S = draws1pop$S, R = draws1pop$R) %>% 
#   ggplot(aes(x = median(S), y = median(R))) +
#   geom_function(fun = ~ SR(SR_fun = "BH", alpha = sim1pop$pars_out$alpha,
#                            Rmax = sim1pop$pars_out$Rmax, S = .x),
#                 aes(lty = "true", col = "true"), lwd = 1) +
#   lapply(sample(4000, 100), function(i) {
#     geom_function(fun = ~ SR(SR_fun = "BH", alpha = alphaRmax$alpha[i],
#                              Rmax = alphaRmax$Rmax[i], S = .x),
#                   aes(lty = "fit", col = "fit"))
#   }) +
#   geom_segment(aes(x = S_true, xend = S_obs, y = R_true, yend = R_obs),
#                col = "slategray4", alpha = 0.3) +
#   geom_segment(aes(x = S_true, xend = median(S), y = R_true, yend = median(R)),
#                col = "slategray4", alpha = 0.3) +
#   geom_point(aes(x = S_true, y = R_true, pch = "true", col = "true"), 
#              fill = "white", size = 2.5) +
#   geom_point(aes(x = S_obs, y = R_obs, pch = "obs", col = "obs"), size = 2.5) +
#   geom_point(aes(pch = "states", col = "states"), size = 2.5) +
#   geom_segment(aes(x = quantile(S, 0.025), xend = quantile(S, 0.975),
#                    y = median(R), yend = median(R), 
#                    lty = "states", col = "states")) +
#   geom_segment(aes(x = median(S), xend = median(S),
#                    y = quantile(R, 0.025), yend = quantile(R, 0.975),
#                    lty = "states", col = "states")) +
#   scale_x_continuous(limits = c(0,NA), expand = c(0,1.05)) +
#   scale_y_continuous(limits = c(0,NA), expand = c(0,1.05)) +
#   scale_shape_manual(values = c(true = 21, obs = 16, states = 16, fit = NA)) +
#   scale_linetype_manual(values = c(true = "dotted", obs = NA, 
#                                    states = "solid", fit = "solid")) +
#   scale_color_manual(values = c(true = "black", obs = "black", 
#                                 states = "slategray4",
#                                 fit = alpha("slategray4", 0.3))) +
#   labs(x = "Spawners", y = "Recruits", shape = "", linetype = "", color = "") +
#   theme(legend.position = c(0.1,0.93))
# ## @knitr
# 
# 
# #-----------------------------------------------------
# # Spawner time series plots
# #-----------------------------------------------------
# 
# ## @knitr singlepop_spawners_ggplot
# draws1pop <- as_draws_rvars(fit1pop) %>% 
#   mutate_variables(S_ppd = rvar_rng(rlnorm, N, log(S), tau))
# 
# sim1pop$sim_dat %>% cbind(S = draws1pop$S, S_ppd = draws1pop$S_ppd) %>% 
#   ggplot(aes(x = year, y = S_obs)) + 
#   geom_line(aes(y = sim1pop$pars_out$S, color = "true", lty = "true"), lwd = 1) +
#   geom_ribbon(aes(ymin = quantile(S, 0.025), ymax = quantile(S, 0.975), fill = "states")) +
#   geom_ribbon(aes(ymin = quantile(S_ppd, 0.025), ymax = quantile(S_ppd, 0.975), fill = "PPD")) +
#   geom_line(aes(y = median(draws1pop$S), lty = "states", col = "states"), lwd = 1) +
#   geom_point(aes(col = "obs", pch = "obs"), size = 3) + 
#   scale_y_continuous(limits = c(0,NA), expand = c(0,0)) +
#   scale_shape_manual(values = c(true = NA, obs = 16, states = NA, PPD = NA)) +
#   scale_color_manual(values = c(true = "black", obs = "black",
#                                 states = "slategray4", PPD = "white")) +
#   scale_fill_manual(values = c(true = "white", obs = "white", 
#                                states = alpha("slategray4", 0.3),
#                                PPD = alpha("slategray4", 0.2))) +
#   scale_linetype_manual(values = c(true = "dotted", obs = NA,  
#                                    states = "solid", PPD = NA)) +
#   labs(x = "Year", y = "Spawners", shape = "", color = "", fill = "", linetype = "") +
#   theme(legend.position = c(0.9,0.93))
# 
# 
# ## @knitr singlepop_spawners_base
# draws1pop <- as_draws_rvars(fit1pop) %>% 
#   mutate_variables(S_ppd = rvar_rng(rlnorm, N, log(S), tau))
# ppd1pop <- sim1pop$sim_dat %>% cbind(S = draws1pop$S, S_ppd = draws1pop$S_ppd) 
# 
# plot(ppd1pop$year, sim1pop$pars_out$S, type = "l", lty = 3, lwd = 3,
#      ylim = c(0, max(quantile(ppd1pop$S_ppd, 0.975))), yaxs = "i",
#      xlab = "Year", ylab = "Spawners", cex.axis = 1.2, cex.lab = 1.5, las = 1)
# polygon(c(ppd1pop$year, rev(ppd1pop$year)),
#         c(quantile(ppd1pop$S, 0.025), rev(quantile(ppd1pop$S, 0.975))),
#         col = alpha("slategray4", 0.3), border = NA)
# polygon(c(ppd1pop$year, rev(ppd1pop$year)),
#         c(quantile(ppd1pop$S_ppd, 0.025), rev(quantile(ppd1pop$S_ppd, 0.975))),
#         col = alpha("slategray4", 0.2), border = NA)
# lines(median(S) ~ year, data = ppd1pop, col = "slategray4", lwd = 3)
# points(S_obs ~ year, data = ppd1pop, pch = 16, cex = 1.5)
# legend("topright", c("true","obs","states","PPD"), cex = 1.2, y.intersp = 1.2,
#        pch = c(NA,16,NA,NA), pt.cex = 1.5, lty = c(3,NA,1,1), lwd = c(3,NA,15,15), 
#        col = c("black", "black", alpha("slategray4", c(0.5,0.2))), bty = "n")
# ## @knitr
# 






