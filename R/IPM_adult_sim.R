#' Simulate data under adult-to-adult integrated population model
#'
#' \code{IPM_adult_sim} simulates initial values for parameters and states in Stan.
#'
#' @param pars Model parameters to be used for simulations.
#' @param pop Population ID.
#' @param year Calendar year.
#' @param X Covariate(s).
#' @param N_age The number of adult age classes.
#' @param max_age Oldest adult age class.
#' @param A Area of spawning habitat.
#' @param F_rate Harvest rate.
#' @param B_rate Broodstock take rate.
#' @param SR_func Character code for the type of stocl-recruit model to fit. At present, only 'BH' for "Beverton-Holt' is allowed.
#' @param n_age_obs Number of adults for age comp.
#' @param n_HW_obs Number of adults for H vs W origin.
#' 
#' @return A list with initial starting values for all of the parameters and states in the Stan model.
#' 
#' @importFrom stats rbinom rlnorm rmultinom rnorm runif
#' 
#' @export
IPM_adult_sim <- function(pars, pop, year, X = NULL, N_age, max_age, A, 
                          F_rate, B_rate, SR_func = "BH", n_age_obs, n_HW_obs)
{
  # spawner-recruit functions
  BH <- function(alpha, Rmax, S, A) 
  {
    R <- alpha*S/(A + alpha*S/Rmax)
    return(R)
  }
  
  N <- length(pop)                        # number of observations 
  N_pop <- max(pop)                       # number of populations
  ages <- (max_age - N_age + 1):max_age   # adult ages
  if(is.null(X)) X <- matrix(0, nrow = max(year), ncol = 1)
  A_pop <- tapply(A, pop, mean)
  
  with(pars, {
    # parameters
    Sigma_alphaRmax <- diag(c(sigma_alpha, sigma_Rmax)^2)
    Sigma_alphaRmax[1,2] <- rho_alphaRmax*sigma_alpha*sigma_Rmax
    Sigma_alphaRmax[2,1] <- Sigma_alphaRmax[1,2]
    aRmax <- exp(mvrnorm(N_pop, c(mu_alpha, mu_Rmax), Sigma_alphaRmax))
    alpha <- aRmax[,1]
    Rmax <- aRmax[,2]
    K <- (alpha - 1)*Rmax/alpha
    phi <- rep(NA, max(year))
    phi[1] <- rnorm(1, 0, sigma_phi/sqrt(1 - rho_phi^2))
    for(i in 2:length(phi))
      phi[i] <- rnorm(1, rho_phi*phi[i-1], sigma_phi)
    phi <- phi + X %*% beta_phi
    mu_alr_p <- log(mu_p[1:(N_age-1)]) - log(mu_p[N_age])
    Sigma_gamma <-  diag(sigma_gamma^2) * (L_gamma %*% t(L_gamma))
    gamma <- mvrnorm(N_pop, mu_alr_p, Sigma_gamma)
    Sigma_alr_p <- diag(sigma_alr_p^2) * (L_alr_p %*% t(L_alr_p))
    alr_p <- t(apply(gamma[pop,], 1, function(x) mvrnorm(1, x, Sigma_alr_p)))
    e_alr_p <- exp(cbind(alr_p, 0))
    p <- sweep(e_alr_p, 1, rowSums(e_alr_p), "/")
    R_init <- data.frame(pop = rep(1:N_pop, each = max_age), year = NA, R = NA)
    for(i in 1:N_pop)
      R_init$year[R_init$pop==i] <- min(year[pop==i]) - (max_age:1)
    R_init$R <- rlnorm(nrow(R_init), log((A_pop*K)[R_init$pop]), sigma_proc)
    alr_p_init <- t(apply(gamma[R_init$pop,], 1, function(x) mvrnorm(1, x, Sigma_alr_p)))
    e_alr_p_init <- exp(cbind(alr_p_init, 0))
    p_init <- sweep(e_alr_p_init, 1, rowSums(e_alr_p_init), "/")
    R_init <- sweep(p_init, 1, R_init$R, "*")
    
    # Simulate recruits and calculate total spawners
    # and spawner age distributions
    S_W_a <- matrix(NA, N, N_age)  # true wild spawners by age  
    S_W <- vector("numeric",N)     # true total wild spawners
    S_H <- vector("numeric",N)     # true total hatchery spawners
    S <- vector("numeric",N)       # true total spawners
    R_hat <- vector("numeric",N)   # expected recruits
    R <- vector("numeric",N)       # true recruits
    B_take <- vector("numeric",N)  # adult broodstock removals
    
    for(i in 1:N)
    {
      for(j in 1:N_age)
        if(year[i] - ages[j] < min(year[pop==pop[i]])) # initialize years 1:max_age
        {
          S_W_a[i,j] <- R_init[R_init$pop==pop[i] & R_init$year==year[i]-ages[j], j]
        } else
        {
          S_W_a[i,j] <- R[i-ages[j]]*p[i-ages[j],j]
        }
      S_W_a[i,-1] <- S_W_a[i,-1]*(1 - F_rate[i])     # catch (assumes no take of age 1)
      B_take[i] <- B_rate[i]*sum(S_W_a[i,-1])
      S_W_a[i,-1] <- S_W_a[i,-1]*(1 - B_rate[i])     # broodstock removal (assumes no take of age 1)
      S_W[i] <- sum(S_W_a[i,])
      S_H[i] <- S_W[i]*p_HOS[i]/(1 - p_HOS[i])
      S[i] <- S_W[i] + S_H[i]
      R_hat[i] <- switch(SR_func,
                             BH = A[i]*BH(alpha[pop[i]], Rmax[pop[i]], S[i], A[i]))
      R[i] <- rlnorm(1, log(R_hat[i]) + phi[year[i]], sigma_proc)
    }
    
    S_obs <- rlnorm(N, log(S), sigma_obs)            # obs total spawners
    q <- sweep(S_W_a, 1, S_W, "/")                   # true spawner age distn 
    n_age_obs <- pmax(round(pmin(n_age_obs, S)), 1)  # cap age samples at pop size
    n_HW_obs <- pmax(round(pmin(n_HW_obs, S)), 1)    # cap H/W samples at pop size
    n_age_obs <- t(sapply(1:N, function(i) rmultinom(1, n_age_obs[i], q[i,]))) # obs wild age frequencies
    dimnames(n_age_obs)[[2]] <- paste0("n_age", ages, "_obs")
    n_H_obs <- rbinom(N, n_HW_obs, p_HOS)            # obs count of hatchery spawners
    n_W_obs <- n_HW_obs - n_H_obs                    # obs count of wild spawners
    
    return(list(sim_dat = data.frame(pop = pop, A = A, year = year, fit_p_HOS = p_HOS > 0,
                                     S_obs = S_obs, n_age_obs, 
                                     n_H_obs = n_H_obs, n_W_obs = n_W_obs, 
                                     B_take_obs = B_take, F_rate = F_rate),
                pars_out = c(pars, list(S_W_a = S_W_a, alpha = alpha, Rmax = Rmax, phi = phi, gamma = gamma, alr_p = alr_p,
                                        p_HOS = p_HOS, p = p, R_hat = R_hat, R = R))))
  })
}
