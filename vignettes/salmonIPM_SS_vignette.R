# Test spawner-spawner IPM with simulated data

library(salmonIPM)
options(device=windows)


#===========================================================================
# SIMULATE DATA
#===========================================================================

# Simulate data
N_pop <- 20
N_year <- 30
N <- N_pop * N_year

sim_out <- IPM_sim(pars = list(mu_alpha = 2, sigma_alpha = 0.5, mu_Rmax = 5, sigma_Rmax = 0.5,
                               rho_alphaRmax = 0.3, beta_phi = 0, rho_phi = 0.7, sigma_phi = 0.5, 
                               sigma = 0.3, tau = 0.5, 
                               mu_p = c(0.05, 0.55, 0.4), sigma_gamma = c(0.1, 0.2), 
                               R_gamma = diag(2), sigma_p = c(0.5, 0.5), R_p = diag(2),
                               S_init_K = 0.7),
                   fish_data = data.frame(pop = rep(1:N_pop, each = N_year),
                                          year = rep(1:N_year, N_pop),
                                          A = 1, p_HOS = 0, F_rate = rbeta(N,7,3), B_rate = 0,
                                          n_age_obs = 50, n_HW_obs = 0),
                   N_age = 3, max_age = 5)

#===========================================================================
# CALL STAN TO FIT MODELS
#===========================================================================

# No pooling across populations
fit_np <- salmonIPM(fish_data = sim_out$sim_dat, stan_model = "IPM_SS_np",
                    chains = 3, iter = 200, warmup = 100, thin = 1, cores = 3,
                    control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 13))

# "old priors" version of partial pooling
fit_pp1 <- salmonIPM(fish_data = sim_out$sim_dat, stan_model = "IPM_SSpa_pp",
                     age_S_obs = rep(1,3), age_S_eff = rep(1,3),
                     chains = 3, iter = 200, warmup = 100, thin = 1, cores = 3,
                     control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 13))

# Partial pooling across populations
fit_pp <- salmonIPM(fish_data = sim_out$sim_dat, stan_model = "IPM_SS_pp",
                    chains = 3, iter = 200, warmup = 100, thin = 1, cores = 3,
                    control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 13))


#===========================================================================
# PLOTS
#===========================================================================

##############################################
# temp: test new priors on initial states
##############################################
library(matrixStats)
library(yarrr)

S_true <- rowSums(sim_out$pars_out$S_W_a) / (1 - sim_out$pars_out$p_HOS)

graphics.off()
windows(width = 14, height = 10)
par(mfrow=c(4,5), mar = c(4,4,2,1))

for(i in 1:20) {
  plot(sim_out$sim_dat$year[sim_out$sim_dat$pop==i], sim_out$sim_dat$S_obs[sim_out$sim_dat$pop==i],
       pch = "", xlab = "year", ylab = "S", main = i, 
       ylim = range(colQuantiles(extract1(fit_pp,"S")[,sim_out$sim_dat$pop==i], probs = c(0.05,0.95))))
  lines(sim_out$sim_dat$year[sim_out$sim_dat$pop==i], S_true[sim_out$sim_dat$pop==i])
  lines(sim_out$sim_dat$year[sim_out$sim_dat$pop==i], colMedians(extract1(fit_pp1,"S")[,sim_out$sim_dat$pop==i]),
        col = "orange", lwd = 2)
  polygon(c(sim_out$sim_dat$year[sim_out$sim_dat$pop==i], rev(sim_out$sim_dat$year[sim_out$sim_dat$pop==i])),
          c(colQuantiles(extract1(fit_pp1,"S")[,sim_out$sim_dat$pop==i], probs = 0.05),
            rev(colQuantiles(extract1(fit_pp1,"S")[,sim_out$sim_dat$pop==i], probs = 0.95))), 
          col = transparent("orange", 0.8), border = NA)
  lines(sim_out$sim_dat$year[sim_out$sim_dat$pop==i], colMedians(extract1(fit_pp,"S")[,sim_out$sim_dat$pop==i]),
        col = "blue", lwd = 2)
  polygon(c(sim_out$sim_dat$year[sim_out$sim_dat$pop==i], rev(sim_out$sim_dat$year[sim_out$sim_dat$pop==i])),
          c(colQuantiles(extract1(fit_pp,"S")[,sim_out$sim_dat$pop==i], probs = 0.05),
            rev(colQuantiles(extract1(fit_pp,"S")[,sim_out$sim_dat$pop==i], probs = 0.95))), 
          col = transparent("blue", 0.8), border = NA)
  points(sim_out$sim_dat$year[sim_out$sim_dat$pop==i], sim_out$sim_dat$S_obs[sim_out$sim_dat$pop==i],
         pch = 16)
}


dev.new(width = 10, height = 7)
par(mfrow = c(1,2))

for(i in c(6,15)) {
  hist(extract1(fit_pp1,"S")[,sim_out$sim_dat$pop==i & sim_out$sim_dat$year==1],
       prob = TRUE, col = transparent("orange",0.8),
       xlim = range(extract1(fit_pp1,"S")[,sim_out$sim_dat$pop==i & sim_out$sim_dat$year==1],
                    extract1(fit_pp,"S")[,sim_out$sim_dat$pop==i & sim_out$sim_dat$year==1]),
       xlab = "S", main = i)
  curve(dlnorm(x,0,10), n = 500, col = "orange", add = TRUE)
  hist(extract1(fit_pp,"S")[,sim_out$sim_dat$pop==i & sim_out$sim_dat$year==1],
       prob = TRUE, col = transparent("blue",0.8), add = TRUE)
  curve(dgamma(x, 1/3*3, 0.001), n = 500, col = "blue", add = TRUE)
  curve(dlnorm(x, log(S_true[sim_out$sim_dat$pop==i & sim_out$sim_dat$year==1]), stan_mean(fit_pp1,"tau")), add = TRUE)
}

#####
#####

# Time series of observed and estimated S_tot under unpooled and partially pooled models, 
# panel for each pop
# png(filename="S_tot_simdata.png", width=16*0.6, height=10*0.6, units="in", res=200, type="cairo-png")
dev.new(width=16,height=10)
par(mfrow=c(3,4), mar=c(1,2,2,1), oma=c(4.1,4.1,0,0))
S_tot <- extract1(stan_BH,"S_tot")
S_tot_fixedpop <- extract1(stan_BH_fixedpop,"S_tot")
c1 <- "darkgray"
c1t <- col2rgb(c1)
c1t <- rgb(c1t[1], c1t[2], c1t[3], maxColorValue = 255, alpha = 255*0.6)
for(i in unique(pop))
{
  plot(year[pop==i], sim_dat$sim_dat$S_tot_obs[pop==i], pch="", cex=1.2, cex.axis=1.2, las=1,
       ylim=c(min(apply(S_tot[,pop==i], 2, quantile, 0.025)),
              max(apply(S_tot[,pop==i], 2, quantile, 0.975))), 
       xlab="", ylab="", log = "y", yaxt = "n")
  at <- maglab(10^par("usr")[3:4], log = T)
  axis(2, at$labat, cex.axis=1.2, las=1,
       labels = sapply(log10(at$labat), function(i) as.expression(bquote(10^ .(i)))))
  lines(year[pop==i], apply(S_tot[,pop==i], 2, quantile, 0.5), lwd=2)
  lines(year[pop==i], apply(S_tot[,pop==i], 2, quantile, 0.025), lwd=1)
  lines(year[pop==i], apply(S_tot[,pop==i], 2, quantile, 0.975), lwd=1)
  lines(year[pop==i], apply(S_tot_fixedpop[,pop==i], 2, quantile, 0.5), col=c1, lwd=2)
  lines(year[pop==i], apply(S_tot_fixedpop[,pop==i], 2, quantile, 0.025), col=c1, lwd=1)
  lines(year[pop==i], apply(S_tot_fixedpop[,pop==i], 2, quantile, 0.975), col=c1, lwd=1)
  # polygon(c(year[pop==i], rev(year[pop==i])),
  #         c(apply(S_tot[,pop==i], 2, quantile, 0.025), rev(apply(S_tot[,pop==i], 2, quantile, 0.975))),
  #         col = c1t, border = NA)
  points(year[pop==i],  sim_dat$sim_dat$S_tot_obs[pop==i], pch = 16, type = "p")
}
mtext("Year", outer = T, side=1, line=2, cex=2*par("cex"))
mtext("Spawners", outer = T, side=2, line=1.1, cex=2*par("cex"))
rm(list=c("S_tot","S_tot_fixedpop","c1","c1t"))
# dev.off()


# Posterior densities under partially pooled model vs. true values
# png(filename="posterior_densities_simdata.png", width=16*0.6, height=10*0.6, units="in", res=200, type="cairo-png")
dev.new(width=16, height=10)
par(mfrow=c(3,4), mar=c(4.1,2,2,1), oma=c(0,4.1,0,0))
mod <- stan_BH
plot(density(extract1(mod,"mu_alpha")), xlab = bquote(mu[a]), ylab="", main="", cex.axis=1.2, cex.lab=2, las=1)
abline(v = sim_dat$pars_out$mu_alpha, lwd = 3)
plot(density(extract1(mod,"sigma_alpha")), xlab = bquote(sigma[a]), ylab="", main="", cex.axis=1.2, cex.lab=2, las=1)
abline(v = sim_dat$pars_out$sigma_alpha, lwd = 3)
plot(density(extract1(mod,"mu_Rmax")), xlab = bquote(mu[b]), ylab="", main="", cex.axis=1.2, cex.lab=2, las=1)
abline(v = sim_dat$pars_out$mu_Rmax, lwd = 3)
plot(density(extract1(mod,"sigma_Rmax")), xlab = bquote(sigma[b]), ylab="", main="", cex.axis=1.2, cex.lab=2, las=1)
abline(v = sim_dat$pars_out$sigma_Rmax, lwd = 3)
plot(density(extract1(mod,"beta_phi")), xlab = bquote(beta[phi]), ylab="", main="", cex.axis=1.2, cex.lab=2, las=1)
abline(v = sim_dat$pars_out$beta_phi, lwd = 3)
plot(density(extract1(mod,"rho_phi")), xlab = bquote(rho[phi]), ylab="", main="", cex.axis=1.2, cex.lab=2, las=1)
abline(v = plogis(sim_dat$pars_out$logit_rho_log_phi), lwd = 3)
plot(density(extract1(mod,"sigma_phi")), xlab = bquote(sigma[phi]), ylab="", main="", cex.axis=1.2, cex.lab=2, las=1)
abline(v = sim_dat$pars_out$sigma_phi, lwd = 3)
plot(density(extract1(mod,"mu_p")[,1]), xlab = bquote(mu[p]), ylab="", main="", cex.axis=1.2, cex.lab=2, las=1, 
     xlim = range(extract1(mod,"mu_p")), col = "lightgray")
abline(v = sim_dat$pars_out$mu_p[1], lwd = 3, col = "lightgray")
dd <- density(extract1(mod,"mu_p")[,2])
lines(dd$x, dd$y, col = "darkgray")
abline(v = sim_dat$pars_out$mu_p[2], lwd = 3, col = "darkgray")
dd <- density(extract1(mod,"mu_p")[,3])
lines(dd$x, dd$y)
abline(v = sim_dat$pars_out$mu_p[3], lwd = 3)
plot(density(extract1(mod,"sigma_p")[,1]), xlab = bquote(sigma[p]), ylab="", main="", 
     cex.axis=1.2, cex.lab=2, las=1, col = "darkgray")
abline(v = sim_dat$pars_out$sigma_p[1], lwd = 3, col = "darkgray")
dd <- density(extract1(mod,"sigma_p")[,2])
lines(dd$x, dd$y)
abline(v = sim_dat$pars_out$sigma_p[2], lwd = 3)
plot(density(extract1(mod,"sigma")), xlab = bquote(mu[sigma[proc]]), ylab="", main="", cex.axis=1.2, cex.lab=2, las=1)
abline(v = sim_dat$pars_out$sigma, lwd = 3)
plot(density(extract1(mod,"sigma_log_sigma_proc")), xlab = bquote(sigma[sigma[proc]]), ylab="", main="", cex.axis=1.2, cex.lab=2, las=1)
abline(v = sim_dat$pars_out$sigma_log_sigma_proc, lwd = 3)
plot(density(extract1(mod,"tau")), xlab = bquote(sigma[obs]), ylab="", main="", cex.axis=1.2, cex.lab=2, las=1)
abline(v = sim_dat$pars_out$tau, lwd = 3)
mtext("Posterior density", side = 2, line = 2, outer = T, cex = 2*par("cex"))
rm(list=c("mod","dd"))
# dev.off()


# Posterior densities of pop-level random effects under unpooled and partially pooled models
# png(filename="random_effects_simdata.png", width=10*0.8, height=10*0.8, units="in", res=200, type="cairo-png")
dev.new(width=10, height=10)
par(mfcol=c(2,1), mar=c(2,5,1,1), oma=c(3.1,0,0,0))
log_a <- log(extract1(stan_BH,"a"))
log_a_fixedpop <- log(extract1(stan_BH_fixedpop,"a"))
log_b <- log(extract1(stan_BH,"b"))
log_b_fixedpop <- log(extract1(stan_BH_fixedpop,"b"))
plot(c(1,ncol(log_a)), range(apply(cbind(log_a, log_a_fixedpop), 2, quantile, c(0.01,0.99))), 
     pch="", ylab = bquote(log(italic(a))), xlab = "", cex.lab=1.5, cex.axis=1.2, las = 1)
for(i in 1:ncol(log_a))
{
  dd1 <- density(log_a_fixedpop[,i])
  # lines(i - 0.4*dd1$y/max(dd1$y), dd1$x, col = "darkgray")
  polygon(c(i - 0.4*dd1$y/max(dd1$y), i, i), c(dd1$x, range(dd1$x)),
          col = "darkgray", border = NA)
  dd2 <- density(log_a[,i])
  # lines(i + 0.4*dd2$y/max(dd2$y), dd2$x)
  polygon(c(i + 0.4*dd2$y/max(dd2$y), i, i), c(dd2$x, range(dd2$x)),
          col = "black", border = NA)
  points(i, log(sim_dat$pars_out$a[i]), pch = 16, col = "red")
}
legend("topleft", c("single-population","multi-population","true value"), pch=c(15,15,16),
       col=c("darkgray","black","red"))
plot(c(1,ncol(log_b)), range(apply(cbind(log_b, log_b_fixedpop), 2, quantile, c(0.01,0.99))), 
     pch="", ylab = bquote(log(italic(b))), xlab = "", cex.lab=1.5, cex.axis=1.2, las = 1)
for(i in 1:ncol(log_b))
{
  dd1 <- density(log_b_fixedpop[,i])
  # lines(i - 0.4*dd1$y/max(dd1$y), dd1$x, col = "darkgray")
  polygon(c(i - 0.4*dd1$y/max(dd1$y), i, i), c(dd1$x, range(dd1$x)),
          col = "darkgray", border = NA)
  dd2 <- density(log_b[,i])
  # lines(i + 0.4*dd2$y/max(dd2$y), dd2$x)
  polygon(c(i + 0.4*dd2$y/max(dd2$y), i, i), c(dd2$x, range(dd2$x)),
          col = "black", border = NA)
  points(i, log(sim_dat$pars_out$b[i]), pch = 16, col = "red")
}
mtext("Population", side=1, outer=T, line=1, cex=1.5*par("cex"))
rm(list=c("log_a","log_a_fixedpop","log_b","log_b_fixedpop","dd1","dd2"))
# dev.off()











