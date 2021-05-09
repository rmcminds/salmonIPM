#' Simulate data from a salmonid integrated population model
#'
#' Generate pseudo-data, group-level parameters, and states from a specified
#' **salmonIPM** integrated population model given values for the hyperparameters.
#' This may be useful, e.g., in simulation-based calibration and model sensitivity
#' checking.
#'
#' @param life_cycle Character string indicating which life-cycle model to
#'   simulate. Currently available options are spawner-to-spawner (`"SS"`,
#'   the default) or spawner-smolt-spawner (`"SMS"`).
#' @param SR_fun One of `"exp"`, `"BH"` (the default), or `"Ricker"`, 
#' indicating which spawner-recruit function to simulate.
#' @param pars Named list of (hyper)parameters to be used for simulations:
#'   * `mu_alpha`  Hyper-mean of log intrinsic productivity.
#'   * `beta_alpha`  Vector of regression coefficients for log intrinsic productivity.
#'   * `sigma_alpha`  Hyper-SD of log intrinsic productivity.
#'   * `mu_Rmax`  If `life_cycle == "SS"`, hyper-mean of log maximum recruitment.
#'   * `beta_Rmax`  If `life_cycle == "SS"`, vector of regression coefficients for 
#'   log maximum recruitment.
#'   * `sigma_Rmax`  If `life_cycle == "SS"`, hyper-SD of log maximum recruitment.
#'   * `rho_alphaRmax`  If `life_cycle == "SS"`, correlation between log(alpha) and log(Rmax).
#'   * `beta_R`  If `life_cycle == "SS"`, vector of regression coefficients for 
#'   log recruitment.
#'   * `rho_R`  If `life_cycle == "SS"`, AR(1) coefficient of brood-year 
#'   log productivity anomalies.
#'   * `sigma_year_R`  If `life_cycle == "SS"`, hyper-SD of brood-year 
#'   log productivity anomalies.
#'   * `sigma_R`  If `life_cycle == "SS"`, unique recruitment process error SD.
#'   * `mu_Mmax`  If `life_cycle == "SMS"`, hyper-mean of log maximum smolt recruitment.
#'   * `beta_Mmax`  If `life_cycle == "SMS"`, vector of regression coefficients for 
#'   log maximum smolt recruitment.
#'   * `sigma_Mmax`  If `life_cycle == "SMS"`, hyper-SD of log maximum smolt recruitment.
#'   * `rho_alphaMmax`  If `life_cycle == "SMS"`, correlation between log(alpha) and log(Mmax).
#'   * `beta_M`  If `life_cycle == "SMS"`, vector of regression coefficients for
#'   log smolt recruitment.
#'   * `rho_M`  If `life_cycle == "SMS"`, AR(1) coefficient of spawner-smolt log
#'   productivity anomalies.
#'   * `sigma_year_M`  If `life_cycle == "SMS"`, process error SD 
#'   of spawner-smolt log productivity anomalies.
#'   * `sigma_M`  If `life_cycle == "SMS"`,
#'   SD of unique spawner-smolt productivity process errors.
#'   * `mu_MS`  If `life_cycle == "SMS"`, mean SAR.
#'   * `beta_MS`  If `life_cycle == "SMS"`, vector of regression coefficients for 
#'   logit SAR anomalies.
#'   * `rho_MS`  If `life_cycle == "SMS"`, AR(1) coefficient for  logit SAR anomalies.
#'   * `sigma_year_MS`  If `life_cycle == "SMS"`, process error SD of logit SAR anomalies.
#'   * `sigma_MS`  If `life_cycle == "SMS"`, SD of unique SAR process errors.
#'   * `mu_p`  Among-population mean simplex of age distributions.
#'   * `sigma_pop_p`  Vector of among-population SDs of mean log-ratio age distributions.
#'   * `R_pop_p`  Among-population correlation matrix of mean log-ratio age distributions. 
#'   * `sigma_pop_p`  Vector of SDs of log-ratio cohort age distributions.
#'   * `R_pop_p`  Correlation matrix of cohort log-ratio age distributions. 
#'   * `sigma_p`  Vector of SDs of log-ratio cohort age distributions.
#'   * `R_p`  Correlation matrix of cohort log-ratio age distributions. 
#'   * `tau_M`  If `life_cycle == "SMS"`, smolt observation error SD.
#'   * `tau_S`  Spawner observation error SD.
#'   * `S_init_K`  Mean of initial spawning population size as a fraction of carrying capacity. 
#' @param par_models  Optional list of two-sided formulas of the form 
#' `theta ~ t1 + ... + tK`, where `theta` is a parameter or state in `pars` 
#'  with corresponding regression coefficients `beta_theta` and `t1 ... tK` 
#'  are terms involving variables in `fish_data`. Standard formula syntax 
#'  such as `:` and `*` may be used; see [stats::formula()].
#' @param scale  Logical indicating whether the model matrices constructed from
#' `fish_data` using the formulas in `par_models` should be scaled to have 
#' column SDs of 1 in addition to being centered (`TRUE`) or centered only (`FALSE`). 
#' @param N_age Number of adult age classes.
#' @param max_age Oldest adult age class.
#' @param ages If `life_cycle != "SS"`, a named list giving the fixed ages
#'   in years of all subadult life stages.
#' @param fish_data Data frame with columns:
#'   * `pop`  Population ID.
#'   * `year`  Calendar year.
#'   * `A`  Spawning habitat area.
#'   * `p_HOS`  True fraction of hatchery-origin spawners.
#'   * `F_rate`  Fishing mortality rate.
#'   * `B_rate`  Hatchery broodstock removal rate.
#'   * `n_age_obs`  Number of spawners sampled for age composition.
#'   * `n_HW_obs`  Number of spawners sampled for hatchery/wild origin.
#'   * `...`  Additional variables to be used as covariates.  
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

sim_salmonIPM <- function(life_cycle = "SS", SR_fun = "BH", pars, par_models = NULL, 
                          scale = TRUE, N_age, max_age, ages = NULL, fish_data)
{
  # Spawner-recruit functions
  SR <- function(SR_fun, alpha, Rmax, S, A) 
  {
    R <- switch(SR_fun,
                exp = alpha*S,
                BH = alpha*S/(1 + alpha*S/(A*Rmax)),
                Ricker = alpha*S*exp(-alpha*S/(A*exp(1)*Rmax)))
    return(R)
  }
  
  # Function to simulate correlated pop-specific intrinsic productivity and max recruitment
  alphaRmax_mvn <- function(N_pop, mu_alpha, sigma_alpha, mu_Rmax, sigma_Rmax, rho_alphaRmax) 
  {
    Sigma_alphaRmax <- diag(c(sigma_alpha, sigma_Rmax)^2)
    Sigma_alphaRmax[1,2] <- rho_alphaRmax*sigma_alpha*sigma_Rmax
    Sigma_alphaRmax[2,1] <- Sigma_alphaRmax[1,2]
    alphaRmax <- matrix(exp(mvrnorm(N_pop, c(mu_alpha, mu_Rmax), Sigma_alphaRmax)), ncol = 2)
    return(list(alpha = alphaRmax[,1], Rmax = alphaRmax[,2]))
  }
  
  # Assign objects
  for(i in 1:length(pars)) assign(names(pars)[i], pars[[i]])
  for(i in 1:ncol(fish_data)) assign(colnames(fish_data)[i], fish_data[,i])
  N <- nrow(fish_data)                  
  N_pop <- max(fish_data$pop)             
  A_pop <- tapply(A, pop, mean)
  pop <- as.numeric(factor(pop))
  year <- as.numeric(factor(year))
  adult_ages <- (max_age - N_age + 1):max_age
  if(life_cycle == "SMS") {
    smolt_age <- ages$M
    ocean_ages <- adult_ages - smolt_age
  }
  
  # Simulate correlated pop-specific intrinsic productivity and max recruitment
  alphaRmax <- alphaRmax_mvn(N_pop, mu_alpha, sigma_alpha, 
                             switch(life_cycle, SS = mu_Rmax, SMS = mu_Mmax),
                             switch(life_cycle, SS = sigma_Rmax, SMS = sigma_Mmax),
                             switch(life_cycle, SS = rho_alphaRmax, SMS = rho_alphaMmax))
  alpha <- alphaRmax$alpha
  if(life_cycle == "SS") Rmax <- alphaRmax$Rmax else Mmax <- alphaRmax$Rmax

  # Carrying capacity (for scaling S_init)
  # assumes BH
  K <- switch(life_cycle, 
              SS = (alpha - 1)*A_pop*Rmax/alpha,
              SMS = mu_MS*(alpha - 1)*A_pop*Mmax/alpha) 
  
  # Covariate model matrices
  X <- par_model_matrix(par_models = par_models, scale = scale, fish_data = fish_data)
  for(i in c("alpha","Rmax","R","Mmax","M","s_MS"))
    assign(paste("X", tail(unlist(strsplit(i, "_")), 1), sep = "_"), 
           if(is.null(X[[i]])) matrix(0,N,1) else X[[i]])

  # Annual recruitment and SAR anomalies
  if(life_cycle == "SS")
  {
    eta_year_R <- rep(NA, max(year))
    eta_year_R[1] <- rnorm(1, 0, sigma_year_R/sqrt(1 - rho_R^2))
    for(i in 2:length(eta_year_R))
      eta_year_R[i] <- rnorm(1, rho_R*eta_year_R[i-1], sigma_year_R)
  }
  if(life_cycle == "SMS")
  {
    eta_year_M <- rep(NA, max(year))
    eta_year_M[1] <- rnorm(1, 0, sigma_year_M/sqrt(1 - rho_M^2))
    eta_year_MS <- rep(NA, max(year))
    eta_year_MS[1] <- rnorm(1, 0, sigma_year_MS/sqrt(1 - rho_MS^2))
    for(i in 2:max(year))
    {
      eta_year_M[i] <- rnorm(1, rho_M*eta_year_M[i-1], sigma_year_M)
      eta_year_MS[i] <- rnorm(1, rho_MS*eta_year_MS[i-1], sigma_year_MS) 
    }
    s_MS <- plogis(rnorm(N, qlogis(mu_MS) + X_MS %*% beta_MS + eta_year_MS[year], sigma_MS))
  }
  
  # Conditional age-at-return
  mu_alr_p <- log(head(mu_p, N_age-1)) - log(tail(mu_p,1))
  Sigma_pop_p <-  tcrossprod(sigma_pop_p) * R_pop_p
  mu_pop_alr_p <- matrix(mvrnorm(N_pop, mu_alr_p, Sigma_pop_p), ncol = N_age - 1)
  Sigma_alr_p <- tcrossprod(sigma_p) * R_p
  alr_p <- t(apply(mu_pop_alr_p[pop,], 1, function(x) mvrnorm(1, x, Sigma_alr_p)))
  exp_alr_p <- cbind(exp(alr_p), 1)
  p <- sweep(exp_alr_p, 1, rowSums(exp_alr_p), "/")
  
  # Initial states
  dat_init <- data.frame(pop = rep(1:N_pop, each = max_age), year = NA, S = NA)
  N_init <- nrow(dat_init)
  for(i in 1:N_pop)
    dat_init$year[dat_init$pop==i] <- min(year[pop==i]) - (max_age:1)
  dat_init$S <- rlnorm(nrow(dat_init), log(S_init_K*K[dat_init$pop]), 
                       ifelse(life_cycle == "SS", tau, tau_S))
  alr_p_init <- t(apply(mu_pop_alr_p[dat_init$pop,], 1, function(x) mvrnorm(1, x, Sigma_alr_p)))
  exp_alr_p_init <- cbind(exp(alr_p_init), 1)
  p_init <- sweep(exp_alr_p_init, 1, rowSums(exp_alr_p_init), "/")
  dat_init$M0 <- switch(life_cycle, 
                        SS = NULL,
                        SMS = SR(SR_fun, alpha[dat_init$pop], Mmax[dat_init$pop], 
                                 dat_init$S, A[dat_init$pop])*rlnorm(N_init, 0, sigma_M))
  dat_init$R <- switch(life_cycle,
                       SS = SR(SR_fun, alpha[dat_init$pop], Rmax[dat_init$pop], 
                               dat_init$S, A[dat_init$pop])*rlnorm(N_init, 0, sigma_R),
                       SMS = dat_init$M0 * plogis(qlogis(mu_MS) + rnorm(N_init, 0, sigma_MS)))
  
  ## Simulate recruits and calculate total spawners and spawner age distributions
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
  
  # calculate true total wild and hatchery spawners and spawner age distribution
  # and predict recruitment from brood year i
  for(i in 1:N)
  {
    # pre-removal spawners by age
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
    
    # smolts and pre-removal spawners by age
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
    
    # catch and broodstock removal (assumes no take of age 1)
    S_W_a[i,-1] <- S_W_a[i,-1]*(1 - F_rate[i]) 
    B_take[i] <- B_rate[i]*sum(S_W_a[i,-1])
    S_W_a[i,-1] <- S_W_a[i,-1]*(1 - B_rate[i]) 
    S_W[i] <- sum(S_W_a[i,])
    S_H[i] <- S_W[i]*p_HOS[i]/(1 - p_HOS[i])
    S[i] <- S_W[i] + S_H[i]
    
    # recruitment
    if(life_cycle=="SS")
    {
      R_hat[i] <- SR(SR_fun, 
                     alpha[pop[i]] * exp(sum(X_alpha[i,]*beta_alpha)), 
                     Rmax[pop[i]] * exp(sum(X_Rmax[i,]*beta_Rmax)), 
                     S[i], A[i])
      R[i] <- R_hat[i] * rlnorm(1, sum(X_R[i,]*beta_R) + eta_year_R[year[i]], sigma_R)
    }
    if(life_cycle=="SMS")
    {
      M_hat[i] <- SR(SR_fun, 
                     alpha[pop[i]] * exp(sum(X_alpha[i,]*beta_alpha)), 
                     Mmax[pop[i]] * exp(sum(X_Mmax[i,]*beta_Mmax)), 
                     S[i], A[i])
      M0[i] <- M_hat[i] * rlnorm(1, sum(X_M[i,]*beta_M) + eta_year_M[year[i]], sigma_M)
    }
  }  # end loop over rows of fish_data
  
  # Observation model
  if(life_cycle == "SS") S_obs <- rlnorm(N, log(S), tau)  # obs total spawners
  if(life_cycle == "SMS")
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
  
  # Return results
  list(
    sim_dat = data.frame(pop = pop, A = A, year = year, fit_p_HOS = p_HOS > 0,
                         S_obs = S_obs, M_obs = switch(life_cycle, SS = NA, SMS = M_obs),
                         n_age_obs, n_H_obs = n_H_obs, n_W_obs = n_W_obs, 
                         B_take_obs = B_take, F_rate = F_rate,
                         fish_data[, names(fish_data) %in% unlist(lapply(par_models, all.vars)), drop = FALSE]),
    pars_out = c(pars, 
                 list(S_W_a = S_W_a, alpha = alpha, 
                      Rmax = switch(life_cycle, SS = Rmax, SMS = NULL), 
                      eta_year_R = switch(life_cycle, SS = eta_year_R, SMS = NULL),
                      R_hat = switch(life_cycle, SS = R_hat, SMS = NULL), 
                      R = switch(life_cycle, SS = R, SMS = NULL),
                      Mmax = switch(life_cycle, SS = NULL, SMS = Mmax), 
                      eta_year_M = switch(life_cycle, SS = NULL, SMS = eta_year_M), 
                      M_hat = switch(life_cycle, SS = NULL, SMS = M_hat),
                      M = switch(life_cycle, SS = NULL, SMS = M),
                      eta_year_MS = switch(life_cycle, SS = NULL, SMS = eta_year_MS),
                      s_MS = switch(life_cycle, SS = NULL, SMS = s_MS),
                      mu_pop_alr_p = mu_pop_alr_p, alr_p = alr_p, p = p, p_HOS = p_HOS))
  )
}
