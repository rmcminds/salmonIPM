#' Simulate data from a salmonid integrated population model
#'
#' Generate pseudo-data, group-level parameters, and states from a specified
#' **salmonIPM** integrated population model given values for the hyperparameters.
#' This may be useful, e.g., in simulation-based calibration and model sensitivity
#' checking.
#'
#' @param life_cycle Character string indicating which life-cycle model to
#'   simulate. Currently available options are spawner-to-spawner (`"SS"`,
#'   the default), iteroparous spawner-to-spawner (`"SSiter"`), or spawner-smolt-spawner (`"SMS"`).
#' @param pars Named list of (hyper)parameters to be used for simulations:
#'   * `mu_alpha`  Hyper-mean of log intrinsic productivity.
#'   * `beta_alpha`  Vector of regression coefficients for log intrinsic productivity.
#'   Must be specified if `par_models` includes `alpha ~ ...`; otherwise will be ignored
#'   and may be omitted.
#'   * `sigma_alpha`  Hyper-SD of log intrinsic productivity.
#'   * `mu_Rmax`  If `life_cycle == "SS"`, hyper-mean of log maximum recruitment.
#'   * `beta_Rmax`  If `life_cycle == "SS"`, vector of regression coefficients for 
#'   log maximum recruitment. Must be specified if `par_models` includes `Rmax ~ ...`; 
#'   otherwise will be ignored and may be omitted.
#'   * `sigma_Rmax`  If `life_cycle == "SS"`, hyper-SD of log maximum recruitment.
#'   * `rho_alphaRmax`  If `life_cycle == "SS"`, correlation between log(alpha) and log(Rmax).
#'   * `beta_R`  If `life_cycle == "SS"`, vector of regression coefficients for 
#'   log recruitment. Must be specified if `par_models` includes `R ~ ...`; 
#'   otherwise will be ignored and may be omitted.
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
#'   log smolt recruitment. Must be specified if `par_models` includes `M ~ ...`; 
#'   otherwise will be ignored and may be omitted.
#'   * `rho_M`  If `life_cycle == "SMS"`, AR(1) coefficient of spawner-smolt log
#'   productivity anomalies.
#'   * `sigma_year_M`  If `life_cycle == "SMS"`, process error SD 
#'   of spawner-smolt log productivity anomalies.
#'   * `sigma_M`  If `life_cycle == "SMS"`,
#'   SD of unique spawner-smolt productivity process errors.
#'   * `mu_MS`  If `life_cycle == "SMS"`, mean SAR.
#'   * `beta_MS`  If `life_cycle == "SMS"`, vector of regression coefficients for 
#'   logit SAR anomalies. Must be specified if `par_models` includes `s_MS ~ ...`; 
#'   otherwise will be ignored and may be omitted.
#'   * `rho_MS`  If `life_cycle == "SMS"`, AR(1) coefficient for logit SAR anomalies.
#'   * `sigma_year_MS`  If `life_cycle == "SMS"`, process error SD of logit SAR anomalies.
#'   * `sigma_MS`  If `life_cycle == "SMS"`, SD of unique logit SAR process errors.
#'   * `mu_p`  Among-population mean simplex of age distributions.
#'   * `sigma_pop_p`  Vector of among-population SDs of mean log-ratio age distributions.
#'   * `R_pop_p`  Among-population correlation matrix of mean log-ratio age distributions. 
#'   * `sigma_pop_p`  Vector of SDs of log-ratio cohort age distributions.
#'   * `R_pop_p`  Correlation matrix of cohort log-ratio age distributions. 
#'   * `sigma_p`  Vector of SDs of log-ratio cohort age distributions.
#'   * `R_p`  Correlation matrix of cohort log-ratio age distributions. 
#'   * `mu_SS`  If `life_cycle == "SSiter"`, mean kelt survival.
#'   * `beta_SS`  If `life_cycle == "SSiter"`, vector of regression coefficients for 
#'   logit kelt survival. Must be specified if `par_models` includes `s_SS ~ ...`; 
#'   otherwise will be ignored and may be omitted.
#'   * `rho_SS`  If `life_cycle == "SSiter"`, AR(1) coefficient for annual logit kelt survival anomalies.
#'   * `sigma_year_SS`  If `life_cycle == "SSiter"`, process error SD of annual logit kelt survival anomalies.
#'   * `sigma_SS`  If `life_cycle == "SSiter"`, SD of unique logit kelt survival process errors.
#'   * `tau`  If `life_cycle != "SMS"`, spawner observation error SD.
#'   * `tau_M`  If `life_cycle == "SMS"`, smolt observation error SD.
#'   * `tau_S`  If `life_cycle == "SMS"`, pawner observation error SD.
#'   * `S_init_K`  Mean of initial spawning population size as a fraction of carrying capacity. 
#' @param N_age Number of (maiden) adult age classes.
#' @param max_age Oldest (maiden) adult age class.
#' @param ages For multi-stage models, a named list giving the fixed ages
#'   in years of all subadult life stages. For example, if `life_cycle == "SMS"`, 
#'   `list(M = a)` where `a` is integer smolt age.
#' @param fish_data Data frame with columns:
#'   * `pop`  Numeric, character or factor population ID.  
#'   * `year`  Numeric or integer giving the year the fish spawned (i.e., the brood year).
#'   * `A`  Spawning habitat size (either stream length or area). Will often be 
#'   time-invariant within a population, but need not be.   
#'   * `p_HOS`  True fraction of hatchery-origin spawners.
#'   * `F_rate`  Total harvest rate (proportion) of natural-origin fish.   
#'   * `B_rate`  Adult hatchery broodstock removal rate.
#'   * `n_age_obs`  Number of spawners sampled for age composition.
#'   * `n_HW_obs`  Number of spawners sampled for hatchery / wild origin.
#'   * `...`  Additional variables to be used as covariates.  
#' @inheritParams salmonIPM
#'
#' @return A named list with elements
#' 
#' * `sim_dat`  A data frame containing simulated data in the structure of 
#' (see [salmonIPM()]) appropriate for the specified `life_cycle`, ready to be passed 
#' to [salmonIPM()].
#' * `pars_out`  A named list of hyperparameters, group-level parameters, and states
#' used in generating the pseudo-data.
#' 
#' @seealso [salmonIPM()] for fitting models
#'
#' @export

simIPM <- function(life_cycle = "SS", SR_fun = "BH", pars, par_models = NULL, 
                   scale = TRUE, N_age, max_age, ages = NULL, fish_data)
{
  # Function to simulate correlated pop-specific intrinsic productivity and max recruitment
  alphaRmax_mvn <- function(N_pop, mu_alpha, sigma_alpha, mu_Rmax, sigma_Rmax, rho_alphaRmax) 
  {
    Sigma_alphaRmax <- diag(c(sigma_alpha, sigma_Rmax)^2)
    Sigma_alphaRmax[1,2] <- rho_alphaRmax*sigma_alpha*sigma_Rmax
    Sigma_alphaRmax[2,1] <- Sigma_alphaRmax[1,2]
    alphaRmax <- matrix(exp(mvrnorm(N_pop, c(mu_alpha, mu_Rmax), Sigma_alphaRmax)), N_pop, 2)
    return(list(alpha = alphaRmax[,1], Rmax = alphaRmax[,2]))
  }
  
  # Assign objects
  if(SR_fun %in% c("B-H","bh","b-h")) SR_fun <- "BH"
  if(SR_fun == "ricker") SR_fun <- "Ricker"
  for(i in 1:length(pars)) assign(names(pars)[i], pars[[i]])
  for(i in 1:ncol(fish_data)) assign(colnames(fish_data)[i], fish_data[,i])
  N <- nrow(fish_data)                  
  N_pop <- length(unique(fish_data$pop))             
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
                             switch(life_cycle, SMS = mu_Mmax, mu_Rmax),
                             switch(life_cycle, SMS = sigma_Mmax, sigma_Rmax),
                             switch(life_cycle, SMS = rho_alphaMmax, rho_alphaRmax))
  alpha <- alphaRmax$alpha
  if(life_cycle == "SMS") Mmax <- alphaRmax$Rmax else Rmax <- alphaRmax$Rmax
  
  # Carrying capacity (for scaling S_init)
  # assumes BH
  K <- switch(life_cycle,
              SMS = mu_MS*(alpha - 1)*A_pop*Mmax/alpha, 
              (alpha - 1)*A_pop*Rmax/alpha)
  
  # Covariate model matrices
  X <- par_model_matrix(par_models = par_models, scale = scale, fish_data = fish_data)
  for(i in c("alpha","Rmax","R","Mmax","M","s_MS","s_SS")) 
  {
    X_i <- paste("X", tail(unlist(strsplit(i, "_")), 1), sep = "_")
    beta_i <- paste("beta", tail(unlist(strsplit(i, "_")), 1), sep = "_")
    if(is.null(X[[i]])) 
    {
      assign(X_i, matrix(0,N,1))
      if(exists(beta_i)) message("Value of ", beta_i, " is being ignored.")
      assign(beta_i, 0)
    } else {
      if(!exists(beta_i)) 
        stop("par_models specifies covariates of ", i, " but pars$", beta_i, " is missing.")
      if(length(get(beta_i)) != ncol(X[[i]]))
        stop("length(", beta_i, ") (", length(get(beta_i)), ") does not equal ncol(", X_i, ") (", ncol(get(X_i)), ").")
      assign(X_i, X[[i]])
    }
  }
  
  # Annual recruitment and SAR anomalies
  if(life_cycle %in% c("SS","SSiter")) 
  {
    eta_year_R <- rep(NA, max(year))
    eta_year_R[1] <- rnorm(1, 0, sigma_year_R/sqrt(1 - rho_R^2))
    for(i in 2:max(year))
      eta_year_R[i] <- rnorm(1, rho_R*eta_year_R[i-1], sigma_year_R)
  }
  
  if(life_cycle == "SSiter") {
    eta_year_SS <- rep(NA, max(year))
    eta_year_SS[1] <- rnorm(1, 0, sigma_year_SS/sqrt(1 - rho_SS^2))
    for(i in 2:max(year))
      eta_year_SS[i] <- rnorm(1, rho_SS*eta_year_SS[i-1], sigma_year_SS)
    s_SS <- plogis(rnorm(N, qlogis(mu_SS) + X_SS %*% beta_SS + eta_year_SS[year], sigma_SS))
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
  
  # Conditional age-at-(maiden)-return
  mu_alr_p <- log(head(mu_p, N_age-1)) - log(tail(mu_p,1))
  Sigma_pop_p <-  tcrossprod(sigma_pop_p) * R_pop_p
  mu_pop_alr_p <- matrix(mvrnorm(N_pop, mu_alr_p, Sigma_pop_p), ncol = N_age - 1)
  Sigma_alr_p <- tcrossprod(sigma_p) * R_p
  alr_p <- t(matrix(apply(mu_pop_alr_p[pop,,drop = FALSE], 1, 
                          function(x) mvrnorm(1, x, Sigma_alr_p)), nrow = N_age - 1))
  exp_alr_p <- cbind(exp(alr_p), 1)
  p <- sweep(exp_alr_p, 1, rowSums(exp_alr_p), "/")
  
  # Initial states
  dat_init <- data.frame(pop = rep(1:N_pop, each = max_age), year = NA, S = NA)
  N_init <- nrow(dat_init)
  for(i in 1:N_pop)
    dat_init$year[dat_init$pop==i] <- min(year[pop==i]) - (max_age:1)
  dat_init$S <- rlnorm(nrow(dat_init), log(S_init_K*K[dat_init$pop]), 
                       ifelse(life_cycle == "SMS", tau_S, tau))
  if(life_cycle == "SMS")
    dat_init$M0 <- SR(SR_fun, alpha[dat_init$pop], Mmax[dat_init$pop], 
                      dat_init$S, A[dat_init$pop]) * rlnorm(N_init, 0, sigma_M)
  dat_init$R <- switch(life_cycle,
                       SMS = dat_init$M0 * plogis(qlogis(mu_MS) + rnorm(N_init, 0, sigma_MS)),
                       SR(SR_fun, alpha[dat_init$pop], Rmax[dat_init$pop], 
                          dat_init$S, A[dat_init$pop]) * rlnorm(N_init, 0, sigma_R))
  alr_p_init <- t(matrix(apply(mu_pop_alr_p[dat_init$pop,,drop = FALSE], 1, 
                               function(x) mvrnorm(1, x, Sigma_alr_p)), nrow = N_age - 1))
  exp_alr_p_init <- cbind(exp(alr_p_init), 1)
  p_init <- sweep(exp_alr_p_init, 1, rowSums(exp_alr_p_init), "/")
  if(life_cycle == "SSiter")
    dat_init$s_SS <- plogis(qlogis(mu_SS) + rnorm(N_init, 0, sigma_SS))
  
  ## Simulate recruits and calculate total spawners and spawner age distributions
  if(life_cycle == "SSiter") {
    S_M_a <- matrix(NA, N, N_age)     # true wild maiden spawners by age  
    S_K_a <- matrix(NA, N, N_age)     # true wild repeat spawners by age  
    S_W_a <- matrix(NA, N, N_age + 1) # true wild total spawners by age
  } else {
    S_W_a <- matrix(NA, N, N_age)  # true wild spawners by age  
  }
  S_W <- vector("numeric",N)     # true total wild spawners
  S_H <- vector("numeric",N)     # true total hatchery spawners
  S <- vector("numeric",N)       # true total spawners
  if(life_cycle=="SMS") {
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
    if(life_cycle == "SS")
    {
      for(a in 1:N_age)
      {
        if(year[i] - adult_ages[a] < min(year[pop == pop[i]])) # use initial states in yrs 1:max_age
        {
          ii <- dat_init$pop == pop[i] & dat_init$year == year[i] - adult_ages[a]
          S_W_a[i,a] <- dat_init$R[ii]*p_init[ii,a]
        } else {  # use process model
          S_W_a[i,a] <- R[i - adult_ages[a]]*p[i - adult_ages[a],a]
        }
      }
    }
    
    # pre-removal maiden and repeat spawners by age
    if(life_cycle == "SSiter")
    {
      for(a in 1:N_age)
      {
        # maiden spawners
        if(year[i] - adult_ages[a] < min(year[pop == pop[i]])) # use initial states in yrs 1:max_age
        {
          ii <- dat_init$pop == pop[i] & dat_init$year == year[i] - adult_ages[a]
          S_M_a[i,a] <- dat_init$R[ii]*p_init[ii,a]
        } else {  # use process model
          S_M_a[i,a] <- R[i - adult_ages[a]]*p[i - adult_ages[a],a]
        }
      }
      
      # repeat spawners
      if(year[i] == min(year[pop == pop[i]])) # ise initial states in year 1
      { 
        ii <- dat_init$pop == pop[i] & dat_init$year == year[i] - 1
        S_K_a[i,] <- dat_init$S[ii]*p_init[ii,]*dat_init$s_SS[ii]  # kludge age structure
      } else {  # use process model
        S_W_plus <- c(head(S_W_a[i-1,], N_age - 1), sum(tail(S_W_a[i-1,], 2)))
        S_K_a[i,] <- S_W_plus*s_SS[i-1]
      }
    }
    
    # smolts and pre-removal spawners by age
    if(life_cycle == "SMS")
    {
      # smolt recruitment
      if(year[i] - smolt_age < min(year[pop == pop[i]])) # use initial states in yrs 1:smolt_age
      {
        M[i] <- dat_init$M0[dat_init$pop == pop[i] & dat_init$year == year[i] - smolt_age]
      } else {
        M[i] <- M0[i - smolt_age]
      }
      
      # adult recruitment
      for(a in 1:N_age)
      {
        if(year[i] - ocean_ages[a] < min(year[pop == pop[i]])) # use initial states in yrs 1:max_ocean_age
        {
          ii <- dat_init$pop == pop[i] & dat_init$year == year[i] - ocean_ages[a]
          S_W_a[i,a] <- dat_init$R[ii]*p_init[ii,a]
        } else {  # use process model
          S_W_a[i,a] <- M[i - ocean_ages[a]]*s_MS[i - ocean_ages[a]]*p[i - ocean_ages[a],a]
        }
      }
    }
    
    # catch and broodstock removal (assumes no take of age 1)
    if(life_cycle == "SSiter")
    {      
      S_M_a[i,-1] <- S_M_a[i,-1]*(1 - F_rate[i]) 
      S_K_a[i,] <- S_K_a[i,]*(1 - F_rate[i]) 
      B_take[i] <- B_rate[i]*(sum(S_M_a[i,-1]) + sum(S_K_a[i,]))
      S_M_a[i,-1] <- S_M_a[i,-1]*(1 - B_rate[i]) 
      S_K_a[i,] <- S_K_a[i,]*(1 - B_rate[i])
      S_W_a[i,] <- c(S_M_a[i,], 0) + c(0, S_K_a[i,])
      S_W[i] <- sum(S_W_a[i,])
    } else {
      S_W_a[i,-1] <- S_W_a[i,-1]*(1 - F_rate[i]) 
      B_take[i] <- B_rate[i]*sum(S_W_a[i,-1])
      S_W_a[i,-1] <- S_W_a[i,-1]*(1 - B_rate[i]) 
      S_W[i] <- sum(S_W_a[i,])
    }
    S_H[i] <- S_W[i]*p_HOS[i]/(1 - p_HOS[i])
    S[i] <- S_W[i] + S_H[i]
    
    # recruitment
    if(life_cycle %in% c("SS","SSiter"))
    {
      R_hat[i] <- SR(SR_fun, 
                     alpha[pop[i]] * exp(sum(X_alpha[i,]*beta_alpha)), 
                     Rmax[pop[i]] * exp(sum(X_Rmax[i,]*beta_Rmax)), 
                     S[i], A[i])
      R[i] <- R_hat[i] * rlnorm(1, sum(X_R[i,]*beta_R) + eta_year_R[year[i]], sigma_R)
    }
    
    if(life_cycle == "SMS")
    {
      M_hat[i] <- SR(SR_fun, 
                     alpha[pop[i]] * exp(sum(X_alpha[i,]*beta_alpha)), 
                     Mmax[pop[i]] * exp(sum(X_Mmax[i,]*beta_Mmax)), 
                     S[i], A[i])
      M0[i] <- M_hat[i] * rlnorm(1, sum(X_M[i,]*beta_M) + eta_year_M[year[i]], sigma_M)
    }
  }  # end loop over rows of fish_data
  
  # Observation model
  if(life_cycle %in% c("SS","SSiter")) 
    S_obs <- rlnorm(N, log(S), tau)    # obs total spawners
  if(life_cycle == "SMS")
  {
    M_obs <- rlnorm(N, log(M), tau_M)  # obs smolts
    S_obs <- rlnorm(N, log(S), tau_S)  # obs total spawners
  }
  n_age_obs <- pmax(round(pmin(n_age_obs, S)), 0)   # cap age samples at pop size
  if(life_cycle == "SSiter")
  {
    q_MK <- sweep(cbind(S_M_a, S_K_a), 1, S_W, "/") # true [maiden | kelt] spawner age distn
    n_MKage_obs <- t(sapply(1:N, function(i) rmultinom(1, n_age_obs[i], q_MK[i,]))) # obs wild age frequencies
    colnames(n_MKage_obs) <- c(paste0("n_age", adult_ages, "M_obs"),
                               paste0("n_age", adult_ages + 1, "K_obs"))
  } else {
    q <- sweep(S_W_a, 1, S_W, "/")                  # true spawner age distn 
    n_age_obs <- t(sapply(1:N, function(i) rmultinom(1, n_age_obs[i], q[i,]))) # obs wild age frequencies
    colnames(n_age_obs) <- paste0("n_age", adult_ages, "_obs")
  }
  n_HW_obs <- pmax(round(pmin(n_HW_obs, S)), 0)     # cap H/W samples at pop size
  n_H_obs <- rbinom(N, n_HW_obs, p_HOS)             # obs count of hatchery spawners
  n_W_obs <- n_HW_obs - n_H_obs                     # obs count of wild spawners
  
  # Return results
  list(
    sim_dat = data.frame(
      pop = fish_data$pop, A = A, year = fish_data$year, 
      cbind(S_obs = S_obs, M_obs = switch(life_cycle, SMS = M_obs, NULL),
            n_age_obs = switch(life_cycle, SSiter = NULL, n_age_obs),
            n_MKage_obs = switch(life_cycle, SSiter = n_MKage_obs, NULL)),
      n_H_obs = n_H_obs, n_W_obs = n_W_obs, 
      fit_p_HOS = p_HOS > 0, B_take_obs = B_take, F_rate = F_rate,
      fish_data[, names(fish_data) %in% unlist(lapply(par_models, all.vars)), drop = FALSE]
    ),
    pars_out = c(pars, 
                 list(M = switch(life_cycle, SMS = M, NULL), 
                      S = S, S_W_a = S_W_a, 
                      q = switch(life_cycle, SSiter = NULL, q),
                      q_MK = switch(life_cycle, SSiter = q_MK, NULL),
                      alpha = alpha, 
                      Rmax = switch(life_cycle, SMS = NULL, Rmax), 
                      eta_year_R = switch(life_cycle, SMS = NULL, eta_year_R),
                      R_hat = switch(life_cycle, SMS = NULL, R_hat), 
                      R = switch(life_cycle, SMS = NULL, R),
                      Mmax = switch(life_cycle, SMS = Mmax, NULL), 
                      eta_year_M = switch(life_cycle, SMS = eta_year_M, NULL), 
                      M_hat = switch(life_cycle, SMS = M_hat, NULL),
                      M = switch(life_cycle, SMS = M, NULL),
                      eta_year_MS = switch(life_cycle, SMS = eta_year_MS, NULL),
                      s_MS = switch(life_cycle, SMS = s_MS, NULL),
                      mu_pop_alr_p = mu_pop_alr_p, alr_p = alr_p, p = p, 
                      eta_year_SS = switch(life_cycle, SSiter = eta_year_SS, NULL),
                      s_SS = switch(life_cycle, SSiter = s_SS, NULL),
                      B_rate = B_rate, p_HOS = p_HOS))
  )
}
