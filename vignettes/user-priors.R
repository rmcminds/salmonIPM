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
library(matrixStats)
library(Hmisc)           # binomial CI function
library(shinystan)       # interactive exploration of posterior
library(distributional)  # plotting priors
library(posterior)       # working with posterior samples
library(ggplot2)         # alpha function
library(viridis)         # plot colors
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
# SPAWNER-RECRUIT PARAMETERS
#===========================================================================

#------------------------------
# Data dimensions
#------------------------------

## @knitr SR_data_setup
set.seed(23148)
N <- 30
N_age <- 3
max_age <- 5
## @knitr

#------------------------------
# True parameter values 
#------------------------------

## @knitr SR_pars
pars1 <- list(mu_alpha = 1.7, sigma_alpha = 0, mu_Rmax = 7, sigma_Rmax = 0, 
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

## @knitr SR_data_struct
df1 <- data.frame(pop = 1, year = 1:N + 2020 - N,
                  A = 1, p_HOS = 0, B_rate = 0,
                  F_rate = rbeta(N, 1, 4),
                  n_age_obs = runif(N, 10, 100),
                  n_HW_obs = 0)

#------------------------------
# Simulate data 
#------------------------------

## @knitr SR_data
sim1 <- simIPM(life_cycle = "SS", SR_fun = "BH", pars = pars1, 
               fish_data = df1, N_age = N_age, max_age = max_age)
sim1$pars_out[c("alpha","Rmax")]
## @knitr

#-----------------------------------------------------
# Fit IPM with default priors
#-----------------------------------------------------

## @knitr SR_default_fit
fit1a <- salmonIPM(life_cycle = "SS", pool_pops = FALSE, SR_fun = "BH", 
                  fish_data = sim1$sim_dat, 
                  seed = 123)

print(fit1a, pars = c("p","R_p","p_HOS","S","R","q","LL"), 
      include = FALSE, prob = c(0.025, 0.5, 0.975))
## @knitr

#-----------------------------------------------------
# Plot posteriors, default priors, and true values
#-----------------------------------------------------

## @knitr SR_default_posteriors
# specify priors using distributional package
log_R_obs <- log(sim1$sim_dat$S_obs / (1 - sim1$sim_dat$F_rate))
prior1a <- c(`log(alpha)` = dist_normal(2,2),
             `log(Rmax)` = dist_normal(quantile(log_R_obs, 0.8, names = FALSE), 
                                       sd(log_R_obs)))

# true parameter values
true1 <- with(sim1$pars_out, c(`log(alpha)` = log(alpha), `log(Rmax)` = log(Rmax)))

# extract and transform draws using posterior package
post1a <- as.data.frame(fit1a, pars = c("alpha","Rmax")) %>% 
  transmute(`log(alpha)` = log(`alpha[1]`), `log(Rmax)` = log(`Rmax[1]`))

# plot
par(mfcol = c(2,2), mar = c(3,3,2,1))
for(i in names(true1)) 
  for(j in names(true1)) {
    if(i == j) {
      hist(post1a[,i], 20, prob = TRUE, col = alpha("slategray4", 0.5), 
           border = "white", xlab = "", ylab = "", yaxt = "n", main = i,
           cex.axis = 1.2, cex.lab = 1.5, cex.main = 1.5, xpd = NA)
      curve(density(prior1a[i], at = x)[[1]], lwd = 0.5, add = TRUE)
      abline(v = true1[i], lwd = 2, lty = 3)
    } else {
      plot(post1a[,i], post1a[,j], pch = 16, col = alpha("slategray4", 0.3),
           xlab = "", ylab = "", las = 1, cex.axis = 1.2, cex.lab = 1.5, xpd = NA)
    }
  }
## @knitr

#-----------------------------------------------------
# Plot true S-R curve, obs, states and fitted draws
#-----------------------------------------------------

## @knitr SR_default_curve
SR <- as_draws_rvars(as.matrix(fit1a, c("S","R")))
RR <- run_recon(sim1$sim_dat)
SRdat <- cbind(RR, S_true = sim1$pars_out$S, R_true = sim1$pars_out$R,
               S = SR$S, R = SR$R)
alphaRmax <- as.data.frame(fit1a, c("alpha", "Rmax")) %>% 
  rename(alpha = `alpha[1]`, Rmax = `Rmax[1]`)

curve(SR(SR_fun = "BH", alpha = sim1$pars_out$alpha,
         Rmax = sim1$pars_out$Rmax, S = x),
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
## @knitr

#-----------------------------------------------------
# Fit IPM with modified prior on alpha
#-----------------------------------------------------

## @knitr SR_user_fit
fit1b <- salmonIPM(life_cycle = "SS", pool_pops = FALSE, SR_fun = "BH", 
                   prior = list(alpha ~ lognormal(2, 0.5)),
                   fish_data = sim1$sim_dat, 
                   seed = 123)

print(fit1b, pars = c("p","R_p","p_HOS","S","R","q","LL"), 
      include = FALSE, prob = c(0.025, 0.5, 0.975))
## @knitr

#-----------------------------------------------------
# Plot posteriors, modified priors, and true values
#-----------------------------------------------------

## @knitr SR_user_posteriors
log_R_obs <- log(sim1$sim_dat$S_obs / (1 - sim1$sim_dat$F_rate))
prior1b <- c(`log(alpha)` = dist_normal(2, 0.5),
             `log(Rmax)` = dist_normal(quantile(log_R_obs, 0.8, names = FALSE), 
                                       sd(log_R_obs)))

# extract and transform draws using posterior package
post1b <- as.data.frame(fit1b, pars = c("alpha","Rmax")) %>% 
  transmute(`log(alpha)` = log(`alpha[1]`), `log(Rmax)` = log(`Rmax[1]`))

# plot
par(mfcol = c(2,2), mar = c(3,3,2,1))
for(i in names(true1)) 
  for(j in names(true1)) {
    if(i == j) {
      hist(post1b[,i], 20, prob = TRUE, col = alpha("slategray4", 0.5), 
           border = "white", xlab = "", ylab = "", yaxt = "n", main = i,
           cex.axis = 1.2, cex.lab = 1.5, cex.main = 1.5, xpd = NA)
      curve(density(prior1b[i], at = x)[[1]], lwd = 0.5, add = TRUE)
      abline(v = true1[i], lwd = 2, lty = 3)
    } else {
      plot(post1b[,i], post1b[,j], pch = 16, col = alpha("slategray4", 0.3),
           xlab = "", ylab = "", las = 1, cex.axis = 1.2, cex.lab = 1.5, xpd = NA)
    }
  }
## @knitr

#-----------------------------------------------------
# Plot true S-R curve, obs, states and fitted draws
#-----------------------------------------------------

## @knitr SR_user_curve
SR <- as_draws_rvars(as.matrix(fit1b, c("S","R")))
RR <- run_recon(sim1$sim_dat)
SRdat <- cbind(RR, S_true = sim1$pars_out$S, R_true = sim1$pars_out$R,
               S = SR$S, R = SR$R)
alphaRmax <- as.data.frame(fit1b, c("alpha", "Rmax")) %>% 
  rename(alpha = `alpha[1]`, Rmax = `Rmax[1]`)

curve(SR(SR_fun = "BH", alpha = sim1$pars_out$alpha,
         Rmax = sim1$pars_out$Rmax, S = x),
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
## @knitr

#===========================================================================
# AGE-STRUCTURE PARAMETERS
#===========================================================================

#------------------------------
# Data dimensions
#------------------------------

## @knitr age_data_setup
set.seed(666)
N <- 30
N_age <- 5
max_age <- 7
## @knitr

#------------------------------
# True parameter values 
#------------------------------

## @knitr age_pars
pars2 <- list(mu_alpha = 1.7, sigma_alpha = 0, mu_Rmax = 7, sigma_Rmax = 0, 
              rho_alphaRmax = 0, rho_R = 0.6, sigma_year_R = 0.2, sigma_R = 0,
              mu_p = c(0.05, 0.4, 0.4, 0.1, 0.05), 
              sigma_pop_p = rep(0,4), R_pop_p = diag(4), 
              sigma_p = c(0.1, 0.2, 0.2, 0.1), R_p = 1 - 0.7*(1 - diag(4)), 
              mu_SS = 0.1, rho_SS = 0.6, sigma_year_SS = 0.2, sigma_SS = 0,
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

## @knitr age_data_struct
df2 <- data.frame(pop = 1, year = 1:N + 2020 - N,
                  A = 1, p_HOS = 0, B_rate = 0,
                  F_rate = rbeta(N, 3, 2),
                  n_age_obs = 0, #runif(N, 10, 100),
                  n_HW_obs = 0)

#------------------------------
# Simulate data 
#------------------------------

## @knitr age_data
sim2 <- simIPM(life_cycle = "SSiter", SR_fun = "BH", pars = pars2, 
               fish_data = df2, N_age = N_age, max_age = max_age)
## @knitr

#-----------------------------------------------------
# Fit IPM with default priors
#-----------------------------------------------------

## @knitr age_default_fit
fit2a <- salmonIPM(life_cycle = "SSiter", pool_pops = FALSE, SR_fun = "BH", 
                   fish_data = sim2$sim_dat, 
                   seed = 123)

print(fit2a, pars = c("p","R_p","s_SS","p_HOS","S","R","q","LL"), 
      include = FALSE, prob = c(0.025, 0.5, 0.975))
## @knitr

#-----------------------------------------------------
# Plot posteriors, modified priors, and true values
#-----------------------------------------------------

## @knitr age_default_posteriors
prior2a <- c(mu_p = dist_beta(1, N_age - 1),
             sigma_p = dist_normal(0,5),
             rho_p = 2*dist_beta((N_age - 1)/2, (N_age - 1)/2) - 1, # LKJ  
             mu_SS = dist_uniform(0,1),
             rho_SS = dist_wrap("gnorm", 0, 0.85, 20),
             sigma_SS = dist_normal(0,3))

# true parameter values
true2 <- sim2$pars_out %>% 
  replace("sigma_SS", .$sigma_year_SS) %>%
  c(rho_p = list(.$R_p[lower.tri(diag(N_age - 1))])) %>%
  .[names(prior2a)] %>% unlist() %>% setNames(gsub("(\\d)", "[\\1]", names(.)))

# extract and transform draws using posterior package
post2a <- as_draws_rvars(fit2a) %>%
  mutate_variables(mu_p = as.vector(mu_p), sigma_p = as.vector(sigma_p),
                   rho_p = R_p[1,,][lower.tri(diag(N_age - 1))]) %>%
  .[names(prior2a)] %>% as_draws_matrix()

# plot
par(mfrow = c(5,4), mar = c(5,1,1,1))
for(j in names(true2)) {
  hist(post2a[,j], 20, prob = TRUE, col = alpha("slategray4", 0.5), border = "white",
       xlab = j, ylab = "", yaxt = "n", main = "", cex.axis = 1.2, cex.lab = 1.4)
  curve(density(prior2a[gsub("\\[.\\]", "", j)], at = x)[[1]], lwd = 0.5, add = TRUE)
  abline(v = true2[[j]], lwd = 2, lty = 3)
}
legend("right", c("true","prior","posterior"), cex = 1.4, text.col = "white",
       fill = c(NA, NA, alpha("slategray4", 0.5)), border = NA,
       inset = c(-1,0), xpd = NA, bty = "n")
legend("right", c("true","prior","posterior"), cex = 1.4,
       lty = c(3,1,NA), lwd = c(2,1,NA), col = c(rep("black",2), "slategray4"),
       inset = c(-1,0), xpd = NA, bty = "n")
## @knitr

#-----------------------------------------------------------
# Spawner age structure time series
#-----------------------------------------------------------

## @knitr age_structure_default
q <- extract1(fit2a, "q")

gg <- sim2$sim_dat %>% 
  select(year, starts_with("n_age")) %>% 
  mutate(total = rowSums(across(starts_with("n_age"))),
         across(starts_with("n_age"), ~ binconf(.x, total, alpha = 0.1))) %>% 
  do.call(data.frame, .) %>% # unpack cols of nested data frames
  pivot_longer(cols = -c(year, total), names_to = c("age", "MK", ".value"),
               names_pattern = "n_age(.)(.)_obs.(.*)") %>% 
  mutate(MK = ifelse(MK == "M", "Maiden", "Repeat")) %>%
  cbind(array(aperm(sapply(1:10, function(k) colQuantiles(q[,,k], probs = c(0.05, 0.5, 0.95)), 
                           simplify = "array"), c(3,1,2)), dim = c(nrow(.), 3), 
              dimnames = list(NULL, paste0("q_age_", c("lo","med","up"))))) %>%
  ggplot(aes(x = year, group = age, color = age, fill = age)) +
  geom_line(aes(y = q_age_med), lwd = 1, alpha = 0.8) +
  geom_ribbon(aes(ymin = q_age_lo, ymax = q_age_up), color = NA, alpha = 0.3) +
  geom_point(aes(y = PointEst), pch = 16, size = 2.5, alpha = 0.8) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0, alpha = 0.8) +
  scale_color_manual(values = viridis(N_age + 1, end = 0.8, direction = -1)) +
  scale_fill_manual(values = viridis(N_age + 1, end = 0.8, direction = -1)) +
  labs(x = "Year", y = "Proportion at age") + 
  facet_wrap(vars(MK), nrow = 2, scales = "free_y") + theme_bw(base_size = 16) + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
        strip.background = element_rect(fill = NA),
        strip.text = element_text(margin = margin(b = 3, t = 3)), 
        legend.box.margin = margin(0,-10,0,-15))

show(gg)  
## @knitr


#-----------------------------------------------------
# Fit IPM with modified prior on mu_p
#-----------------------------------------------------

## @knitr age_default_fit
fit2b <- salmonIPM(life_cycle = "SSiter", pool_pops = FALSE, SR_fun = "BH", 
                   prior = list(mu_p ~ dirichlet(sim2$pars_out$mu_p * 100),
                                mu_SS ~ beta(10,90)),
                   fish_data = sim2$sim_dat, 
                   seed = 123)

print(fit2b, pars = c("p","R_p","s_SS","p_HOS","S","R","q","LL"), 
      include = FALSE, prob = c(0.025, 0.5, 0.975))
## @knitr

#-----------------------------------------------------------
# Spawner age structure time series
#-----------------------------------------------------------

## @knitr age_structure_default
q <- extract1(fit2b, "q")

gg <- sim2$sim_dat %>% 
  select(year, starts_with("n_age")) %>% 
  mutate(total = rowSums(across(starts_with("n_age"))),
         across(starts_with("n_age"), ~ binconf(.x, total, alpha = 0.1))) %>% 
  do.call(data.frame, .) %>% # unpack cols of nested data frames
  pivot_longer(cols = -c(year, total), names_to = c("age", "MK", ".value"),
               names_pattern = "n_age(.)(.)_obs.(.*)") %>% 
  mutate(MK = ifelse(MK == "M", "Maiden", "Repeat")) %>%
  cbind(array(aperm(sapply(1:10, function(k) colQuantiles(q[,,k], probs = c(0.05, 0.5, 0.95)), 
                           simplify = "array"), c(3,1,2)), dim = c(nrow(.), 3), 
              dimnames = list(NULL, paste0("q_age_", c("lo","med","up"))))) %>%
  ggplot(aes(x = year, group = age, color = age, fill = age)) +
  geom_line(aes(y = q_age_med), lwd = 1, alpha = 0.8) +
  geom_ribbon(aes(ymin = q_age_lo, ymax = q_age_up), color = NA, alpha = 0.3) +
  geom_point(aes(y = PointEst), pch = 16, size = 2.5, alpha = 0.8) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0, alpha = 0.8) +
  scale_color_manual(values = viridis(N_age + 1, end = 0.8, direction = -1)) +
  scale_fill_manual(values = viridis(N_age + 1, end = 0.8, direction = -1)) +
  labs(x = "Year", y = "Proportion at age") + 
  facet_wrap(vars(MK), nrow = 2, scales = "free_y") + theme_bw(base_size = 16) + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
        strip.background = element_rect(fill = NA),
        strip.text = element_text(margin = margin(b = 3, t = 3)), 
        legend.box.margin = margin(0,-10,0,-15))

show(gg)  
## @knitr





####################
par_names <- c("log(alpha)","log(Rmax)","rho_R","sigma_R",
               "mu_p","sigma_p","rho_p", "mu_SS","rho_SS","sigma_SS","tau")

# specify priors using distributional package
log_S_obs <- log(sim2$sim_dat$S_obs)
concentration <- sim2$pars_out$mu_p * 100
prior2 <- c(`log(alpha)` = dist_normal(2,2),
            `log(Rmax)` = dist_normal(max(log_S_obs), sd(log_S_obs)),
            rho_R = dist_wrap("gnorm", 0, 0.85, 20),
            sigma_R = dist_normal(0,5),
            # mu_p = dist_beta(1,4),
            setNames(dist_beta(concentration, 100 - concentration), paste0("mu_p[", 1:5, "]")),
            # sigma_p = dist_normal(0,5),
            setNames(rep(dist_normal(0,5), 4), paste0("sigma_p[", 1:4, "]")),
            rho_p = 2*dist_beta((N_age - 1)/2, (N_age - 1)/2) - 1, # LKJ  
            mu_SS = dist_beta(10,90), #dist_uniform(0,1),
            rho_SS = dist_wrap("gnorm", 0, 0.85, 20),
            sigma_SS = dist_normal(0,3),
            tau = dist_wrap("gnorm", 1, 0.85, 30))

# true parameter values
true2 <- sim2$pars_out %>% 
  replace(c("sigma_R", "sigma_SS"), c(.$sigma_year_R, .$sigma_year_SS)) %>%
  c(`log(alpha)` = log(.$alpha), `log(Rmax)` = log(.$Rmax), rho_p = .$R_p[2,1]) %>%
  .[par_names] %>% unlist() %>% setNames(gsub("(\\d)", "[\\1]", names(.)))

# extract and transform draws using posterior package
post2 <- as_draws_rvars(fit2b) %>%
  mutate_variables(`log(alpha)` = log(alpha), `log(Rmax)` = log(Rmax),
                   mu_p = as.vector(mu_p), sigma_p = as.vector(sigma_p),
                   rho_p = as.vector(R_p[1,2,1])) %>%
  .[par_names] %>% as_draws_matrix()

# plot
par(mfrow = c(5,4), mar = c(5,1,1,1))
for(j in names(true2)) {
  hist(post2[,j], 20, prob = TRUE, col = alpha("slategray4", 0.5), border = "white",
       xlab = j, ylab = "", yaxt = "n", main = "", cex.axis = 1.2, cex.lab = 1.4)
  # curve(density(prior2[gsub("\\[.\\]", "", j)], at = x)[[1]], lwd = 0.5, add = TRUE)
  curve(density(prior2[j], at = x)[[1]], lwd = 0.5, add = TRUE)
  abline(v = true2[[j]], lwd = 2, lty = 3)
}
legend("right", c("true","prior","posterior"), cex = 1.4, text.col = "white",
       fill = c(NA, NA, alpha("slategray4", 0.5)), border = NA,
       inset = c(-1,0), xpd = NA, bty = "n")
legend("right", c("true","prior","posterior"), cex = 1.4,
       lty = c(3,1,NA), lwd = c(2,1,NA), col = c(rep("black",2), "slategray4"),
       inset = c(-1,0), xpd = NA, bty = "n")
####################
