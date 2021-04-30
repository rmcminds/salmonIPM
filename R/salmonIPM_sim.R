#' Simulate data under an integrated population model
#'
#' Generate pseudo-data, group-level parameters, and states from a specified
#' **salmonIPM** integrated population model given values for the hyperparameters.
#' This may be useful, e.g., in simulation-based calibration and model sensitivity
#' checking.
#'
#' @param pars Named list of (hyper)parameters to be used for simulations:
#'   * `mu_alpha`  Hyper-mean of log intrinsic productivity.
#'   * `sigma_alpha`  Hyper-SD of log intrinsic productivity.
#'   * `mu_Rmax`  Hyper-mean of log asymptotic recruitment.
#'   * `sigma_Rmax`  Hyper-SD of log asymptotic recruitment.
#'   * `rho_alphaRmax`  Correlation between log(alpha) and log(Rmax).
#'   * `beta_R`  If `life_cycle=="SS"`, vector of regression
#'   coefficients for log productivity anomalies.
#'   * `rho_R`  If `life_cycle=="SS"`, AR(1) coefficient for 
#'   log productivity anomalies.
#'   * `sigma_year_`  If `life_cycle=="SS"`, hyper-SD of brood year 
#'   log productivity anomalies.
#'   * `sigma_R`  If `life_cycle=="SS"`, unique recruitment process error SD.
#'   * `beta_M`  If `life_cycle=="SMS"`, vector of regression coefficients for
#'   spawner-smolt log productivity anomalies.
#'   * `rho_year_M`  If `life_cycle=="SMS"`, AR(1) coefficient of spawner-smolt log
#'   productivity anomalies.
#'   * `sigma_year_M`  If `life_cycle=="SMS"`, process error SD 
#'   of spawner-smolt log productivity anomalies.
#'   * `sigma_M`  If `life_cycle=="SMS"`,
#'   SD of unique spawner-smolt productivity process errors.
#'   * `mu_MS`  If `life_cycle=="SMS"`, mean SAR.
#'   * `beta_MS`  If `life_cycle=="SMS"`, vector of regression
#'   coefficients for logit SAR anomalies.
#'   * `rho_MS`  If `life_cycle=="SMS"`, AR(1) coefficient for  logit SAR anomalies.
#'   * `sigma_year_MS`  If `life_cycle=="SMS"`, process error SD of logit SAR anomalies.
#'   * `sigma_MS`  If `life_cycle=="SMS"`, SD of unique SAR process errors.
#'   * `mu_p`  Among-population mean simplex of age distributions.
#'   * `sigma_pop_p`  Vector of among-population SDs of mean log-ratio age distributions.
#'   * `R_pop_p`  Among-population correlation matrix of mean log-ratio age distributions. 
#'   * `sigma_pop_p`  Vector of SDs of log-ratio cohort age distributions.
#'   * `R_pop_p`  Correlation matrix of cohort log-ratio age distributions. 
#'   * `sigma_p`  Vector of SDs of log-ratio cohort age distributions.
#'   * `R_p`  Correlation matrix of cohort log-ratio age distributions. 
#'   * `tau_M`  If `life_cycle=="SMS"`, smolt observation error SD.
#'   * `tau_S`  Spawner observation error SD.
#'   * `S_init_K`  Mean of initial spawning population size as a fraction of carrying capacity. 
#' @param fish_data Data frame with columns:
#'   * `pop`  Population ID.
#'   * `year`  Calendar year.
#'   * `A`  Spawning habitat area.
#'   * `p_HOS`  True fraction of hatchery-origin spawners.
#'   * `F_rate`  Fishing mortality rate.
#'   * `B_rate`  Hatchery broodstock removal rate.
#'   * `n_age_obs`  Number of spawners sampled for age composition.
#'   * `n_HW_obs`  Number of spawners sampled for hatchery/wild origin.
#' @param env_data Optional data frame or named list of data frames whose
#'   variables are time-varying environmental covariates, sequentially ordered
#'   with each row corresponding to a unique year in fish_data. If a named list,
#'   element names correspond to stage- or transition-specific covariate
#'   matrices defined in the simulation model being used. (This is required if
#'   `life_cycle != "SS"`.)
#' @param life_cycle Character string indicating which life-cycle model to
#'   simulate Currently available options are spawner-to-spawner (`"SS"`,
#'   the default) or spawner-smolt-spawner (`"SMS"`).
#' @param N_age The number of adult age classes.
#' @param max_age Oldest adult age class.
#' @param ages If `life_cycle != "SS"`, a named list giving the fixed ages
#'   in years of all subadult life stages.
#' @param SR_fun One of `"exp"`, `"BH"` (the default), or
#'   `"Ricker"`, indicating which spawner-recruit function to simulate.
#'
#' @return A named list with elements
#' 
#' * `sim_dat`  A data frame containing simulated data in the structure of `fish_data` 
#' (see [stan_data()]) appropriate for the specified `life_cycle`, ready to be passed 
#' to [salmonIPM()].
#' * `pars_out`  A named list of hyperparameters, group-level parameters, and states
#' used in generating the pseudo-data.
#' 
#'
#' @importFrom stats rbinom rlnorm rmultinom rnorm runif
#'
#' @export

salmonIPM_sim <- function(pars, fish_data, env_data = NULL, life_cycle = "SS", 
                          N_age, max_age, ages = NULL, SR_fun = "BH")
{
  # spawner-recruit functions
  SR <- function(SR_fun, alpha, Rmax, S, A) 
  {
    R <- switch(SR_fun,
                exp = alpha*S/A,
                BH = alpha*S/(A + alpha*S/Rmax),
                Ricker = alpha*(S/A)*exp(-alpha*S/(A*exp(1)*Rmax)))
    return(R)
  }
  
  for(i in 1:length(pars)) assign(names(pars)[i], pars[[i]])
  for(i in 1:ncol(fish_data)) assign(colnames(fish_data)[i], fish_data[,i])
  N <- nrow(fish_data)                  
  N_pop <- max(fish_data$pop)             
  A_pop <- tapply(A, pop, mean)
  pop <- as.numeric(factor(pop))
  year <- as.numeric(factor(year))
  adult_ages <- (max_age - N_age + 1):max_age
  if(life_cycle == "SMS") 
  {
    smolt_age <- ages$M
    ocean_ages <- adult_ages - smolt_age
  }
  if(is.null(env_data)) env_data <- switch(life_cycle,
                                           SS = matrix(0, max(fish_data$year), 1),
                                           SMS = list(M = matrix(0, max(fish_data$year), 1),
                                                      MS = matrix(0, max(fish_data$year), 1)))
  # parameters
  Sigma_alphaRmax <- diag(c(sigma_alpha, sigma_Rmax)^2)
  Sigma_alphaRmax[1,2] <- rho_alphaRmax*sigma_alpha*sigma_Rmax
  Sigma_alphaRmax[2,1] <- Sigma_alphaRmax[1,2]
  alphaRmax <- matrix(exp(mvrnorm(N_pop, c(mu_alpha, mu_Rmax), Sigma_alphaRmax)), ncol = 2)
  alpha <- alphaRmax[,1]
  Rmax <- alphaRmax[,2]
  K <- A_pop*switch(life_cycle, 
                    SS = (alpha - 1)*Rmax/alpha,
                    SMS = mu_MS*(alpha - 1)*Rmax/alpha) # assumes BH
  if(life_cycle == "SS")
  {
    eta_year_R <- rep(NA, max(year))
    eta_year_R[1] <- rnorm(1, 0, sigma_year_R/sqrt(1 - rho_R^2))
    for(i in 2:length(eta_year_R))
      eta_year_R[i] <- rnorm(1, rho_R*eta_year_R[i-1], sigma_year_R)
    eta_year_R <- eta_year_R + env_data %*% beta_R
  }
  if(life_cycle == "SMS")
  {
    eta_year_M <- rep(NA, max(year))
    eta_year_MS <- rep(NA, max(year))
    eta_year_M[1] <- rnorm(1, 0, sigma_year_M/sqrt(1 - rho_M^2))
    eta_year_MS[1] <- rnorm(1, 0, sigma_year_MS/sqrt(1 - rho_MS^2))
    for(i in 2:max(year))
    {
      eta_year_M[i] <- rnorm(1, rho_M*eta_year_M[i-1], sigma_year_M)
      eta_year_MS[i] <- rnorm(1, rho_MS*eta_year_MS[i-1], sigma_year_MS) 
    }
    eta_year_M <- eta_year_M + env_data$M %*% beta_M
    eta_year_MS <- eta_year_MS + env_data$MS %*% beta_MS
    s_MS <- plogis(rnorm(N, qlogis(mu_MS) + eta_year_MS[year], sigma_MS))
  }
  mu_alr_p <- log(mu_p[1:(N_age-1)]) - log(mu_p[N_age])
  Sigma_pop_p <-  (sigma_pop_p %*% t(sigma_pop_p)) * R_pop_p
  mu_pop_alr_p <- matrix(mvrnorm(N_pop, mu_alr_p, Sigma_pop_p), ncol = N_age - 1)
  Sigma_alr_p <- (sigma_p %*% t(sigma_p)) * R_p
  alr_p <- t(apply(mu_pop_alr_p[pop,], 1, function(x) mvrnorm(1, x, Sigma_alr_p)))
  e_alr_p <- exp(cbind(alr_p, 0))
  p <- sweep(e_alr_p, 1, rowSums(e_alr_p), "/")
  
  dat_init <- data.frame(pop = rep(1:N_pop, each = max_age), year = NA, S = NA)
  N_init <- nrow(dat_init)
  for(i in 1:N_pop)
    dat_init$year[dat_init$pop==i] <- min(year[pop==i]) - (max_age:1)
  dat_init$S <- rlnorm(nrow(dat_init), log(S_init_K*K[dat_init$pop]), 
                       ifelse(life_cycle=="SS", tau, tau_S))
  alr_p_init <- t(apply(mu_pop_alr_p[dat_init$pop,], 1, function(x) mvrnorm(1, x, Sigma_alr_p)))
  e_alr_p_init <- exp(cbind(alr_p_init, 0))
  p_init <- sweep(e_alr_p_init, 1, rowSums(e_alr_p_init), "/")
  dat_init$M0 <- switch(life_cycle, 
                        SS = NULL,
                        SMS = A[dat_init$pop]*SR(SR_fun, alpha[dat_init$pop], Rmax[dat_init$pop], 
                                                 dat_init$S, A[dat_init$pop])*rlnorm(N_init, 0, sigma_M))
  dat_init$R <- switch(life_cycle,
                       SS = A[dat_init$pop]*SR(SR_fun, alpha[dat_init$pop], Rmax[dat_init$pop], 
                                               dat_init$S, A[dat_init$pop])*rlnorm(N_init, 0, sigma_R),
                       SMS = dat_init$M0*plogis(qlogis(mu_MS) + rnorm(N_init, 0, sigma_MS)))
  
  # Simulate recruits and calculate total spawners and spawner age distributions
  S_W_a <- matrix(NA, N, N_age)  # true wild spawners by age  
  S_W <- vector("numeric",N)     # true total wild spawners
  S_H <- vector("numeric",N)     # true total hatchery spawners
  S <- vector("numeric",N)       # true total spawners
  if(life_cycle=="SMS")
  {
    M_hat <- vector("numeric",N) # expected smolts
    M0 <- vector("numeric",N)    # true smolts by brood year
    M <- vector("numeric",N)     # true smolts by calendar year
  }
  R_hat <- vector("numeric",N)   # expected recruits
  R <- vector("numeric",N)       # true recruits
  B_take <- vector("numeric",N)  # adult broodstock removals
  
  for(i in 1:N)
  {
    if(life_cycle=="SS")
    {
      for(a in 1:N_age)
      {
        if(year[i] - adult_ages[a] < min(year[pop==pop[i]])) # initialize spawners in yrs 1:max_age
        {
          ii <- dat_init$pop == pop[i] & dat_init$year == year[i] - adult_ages[a]
          S_W_a[i,a] <- dat_init$R[ii]*p_init[ii,a]
        } else {
          S_W_a[i,a] <- R[i - adult_ages[a]]*p[i - adult_ages[a],a]
        }
      }
    }
    
    if(life_cycle=="SMS")
    {
      if(year[i] - smolt_age < min(year[pop==pop[i]])) # initialize smolts in yrs 1:smolt_age
      {
        M[i] <- dat_init$M0[dat_init$pop==pop[i] & dat_init$year==year[i] - smolt_age]
      } else {
        M[i] <- M0[i - smolt_age]
      }
      
      for(a in 1:N_age)
      {
        if(year[i] - ocean_ages[a] < min(year[pop==pop[i]])) # initialize spawners in yrs 1:max_ocean_age
        {
          ii <- dat_init$pop == pop[i] & dat_init$year == year[i] - ocean_ages[a]
          S_W_a[i,a] <- dat_init$R[ii]*p_init[ii,a]
        } else {
          S_W_a[i,a] <- M[i - ocean_ages[a]]*s_MS[i - ocean_ages[a]]*p[i - ocean_ages[a],a]
        }
      }
    }
    
    S_W_a[i,-1] <- S_W_a[i,-1]*(1 - F_rate[i]) # catch (assumes no take of age 1)
    B_take[i] <- B_rate[i]*sum(S_W_a[i,-1])
    S_W_a[i,-1] <- S_W_a[i,-1]*(1 - B_rate[i]) # broodstock removal (assumes no take of age 1)
    S_W[i] <- sum(S_W_a[i,])
    S_H[i] <- S_W[i]*p_HOS[i]/(1 - p_HOS[i])
    S[i] <- S_W[i] + S_H[i]
    
    if(life_cycle=="SS")
    {
      R_hat[i] <- A[i]*SR(SR_fun, alpha[pop[i]], Rmax[pop[i]], S[i], A[i])
      R[i] <- R_hat[i]*rlnorm(1, eta_year_R[year[i]], sigma_R)
    }
    if(life_cycle=="SMS")
    {
      M_hat[i] <- A[i]*SR(SR_fun, alpha[pop[i]], Rmax[pop[i]], S[i], A[i])
      M0[i] <- M_hat[i]*rlnorm(1, eta_year_M[year[i]], sigma_M)
    }
  }
  
  if(life_cycle=="SS") S_obs <- rlnorm(N, log(S), tau)  # obs total spawners
  if(life_cycle=="SMS")
  {
    M_obs <- rlnorm(N, log(M), tau_M)              # obs smolts
    S_obs <- rlnorm(N, log(S), tau_S)              # obs total spawners
  }
  q <- sweep(S_W_a, 1, S_W, "/")                   # true spawner age distn 
  n_age_obs <- pmax(round(pmin(n_age_obs, S)), 0)  # cap age samples at pop size
  n_HW_obs <- pmax(round(pmin(n_HW_obs, S)), 0)    # cap H/W samples at pop size
  n_age_obs <- t(sapply(1:N, function(i) rmultinom(1, n_age_obs[i], q[i,]))) # obs wild age frequencies
  dimnames(n_age_obs)[[2]] <- paste0("n_age", adult_ages, "_obs")
  n_H_obs <- rbinom(N, n_HW_obs, p_HOS)            # obs count of hatchery spawners
  n_W_obs <- n_HW_obs - n_H_obs                    # obs count of wild spawners
  
  return(list(sim_dat = data.frame(pop = pop, A = A, year = year, fit_p_HOS = p_HOS > 0,
                                   S_obs = S_obs, M_obs = switch(life_cycle, SS = NA, SMS = M_obs),
                                   n_age_obs, n_H_obs = n_H_obs, n_W_obs = n_W_obs, 
                                   B_take_obs = B_take, F_rate = F_rate),
              pars_out = c(pars, 
                           list(S_W_a = S_W_a, alpha = alpha, Rmax = Rmax, 
                                eta_year_R = switch(life_cycle, SS = eta_year_R, SMS = NULL),
                                eta_year_M = switch(life_cycle, SS = NULL, SMS = eta_year_M),
                                eta_year_MS = switch(life_cycle, SS = NULL, SMS = eta_year_MS),
                                s_MS = switch(life_cycle, SS = NULL, SMS = s_MS),
                                mu_pop_alr_p = mu_pop_alr_p, alr_p = alr_p, p = p, p_HOS = p_HOS, 
                                R_hat = switch(life_cycle, SS = R_hat, SMS = NULL), 
                                R = switch(life_cycle, SS = R, SMS = NULL),
                                M_hat = switch(life_cycle, SS = NULL, SMS = M_hat),
                                M = switch(life_cycle, SS = NULL, SMS = M)))))
  
}
