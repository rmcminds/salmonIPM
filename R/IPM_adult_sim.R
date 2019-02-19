#' Simulate data under an integrated population model
#'
#' \code{IPM_adult_sim} simulates initial values for parameters and states in
#' Stan.
#'
#' @param pars Named list of (hyper-) parameters to be used for simulations:
#'   \describe{\item{\code{mu_alpha}}{Hyper-mean of log intrinsic productivity}
#'   \item{\code{sigma_alpha}}{Hyper-SD of log intrinsic productivity}
#'   \item{\code{mu_Rmax}}{Hyper-mean of log asymptotic recruitment}
#'   \item{\code{sigma_Rmax}}{Hyper-SD of log asymptotic recruitment}
#'   \item{\code{rho_alphaRmax}}{Correlation between log(alpha) and log(Rmax)}
#'   \item{\code{beta_phi}}{Vector of regression coefficients for log
#'   productivity anomalies} \item{\code{rho_phi}}{AR(1) coefficient for log
#'   productivity anomalies} \item{\code{sigma_phi}}{Hyper-SD of brood year log
#'   productivity anomalies} \item{\code{sigma}}{Unique recruitment process
#'   error SD}  \item{\code{beta_M}}{If \code{life_cycle=="SMS"}, vector of
#'   regression coefficients for spawner-smolt productivity}
#'   \item{\code{rho_M}}{If \code{life_cycle=="SMS"}, AR(1) coefficient of
#'   spawner-smolt productivity} \item{\code{sigma_M}}{If
#'   \code{life_cycle=="SMS"}, spawner-smolt process error SD}
#'   \item{\code{mu_MS}}{If \code{life_cycle=="SMS"}, mean SAR}
#'   \item{\code{beta_MS}}{If \code{life_cycle=="SMS"}, vector of regression
#'   coefficients for SAR} \item{\code{rho_MS}}{If \code{life_cycle=="SMS"},
#'   AR(1) coefficient for logit SAR} \item{\code{sigma_MS}}{If
#'   \code{life_cycle=="SMS"}, process error SD of SAR}
#'   \item{\code{mu_p}}{Among-population mean simplex of age distributions}
#'   \item{\code{sigma_gamma}}{Vector of among-population SDs of mean log-ratio
#'   age distributions} \item{\code{R_gamma}}{Among-population correlation
#'   matrix of mean log-ratio age distributions}
#'   \item{\code{sigma_gamma}}{Vector of SDs of log-ratio cohort age
#'   distributions} \item{\code{R_gamma}}{Correlation matrix of cohort log-ratio
#'   age distributions} \item{\code{tau_M}}{If \code{life_cycle=="SMS"}, smolt
#'   observation error SD} \item{\code{tau_S}}{Spawner observation error SD}
#'   \item{\code{S_init_K}}{Mean of initial spawning population size as a
#'   fraction of carrying capacity} }
#' @param fish_data Data frame with columns
#'   \describe{\item{\code{pop}}{Population ID} \item{\code{year}}{Calendar
#'   year} \item{\code{A}}{Spawning habitat area}  \item{\code{F_rate}}{Fishing
#'   mortality rate} \item{\code{B_rate}}{Hatchery broodstock removal rate}
#'   \item{\code{n_age_obs}}{Number of spawners sampled for age composition}
#'   \item{\code{n_HW_obs}}{Number of spawners sampled for hatchery/wild
#'   origin}}
#' @param env_data Optional data frame or named list of data frames whose
#'   variables are time-varying environmental covariates, sequentially ordered
#'   with each row corresponding to a unique year in fish_data. If a named list,
#'   element names correspond to stage- or transition-specific covariate
#'   matrices defined in the simulation model being used. (This is required if
#'   \code{life_cycle != "SS"}.)
#' @param life_cycle Character string indicating which life-cycle model to
#'   simulate Currently available options are spawner-to-spawner (\code{"SS"},
#'   the default) or spawner-smolt-spawner (\code{"SMS"}).
#' @param N_age The number of adult age classes.
#' @param max_age Oldest adult age class.
#' @param ages If \code{life_cycle != "SS"}, a named list giving the fixed ages
#'   in years of all subadult life stages.
#' @param SR_fun One of \code{"exp"}, \code{"BH"} (the default), or
#'   \code{"Ricker"}, indicating which spawner-recruit function to simulate.
#'
#' @return A list with initial starting values for all of the parameters and
#'   states in the Stan model.
#'
#' @importFrom stats rbinom rlnorm rmultinom rnorm runif
#'
#' @export

IPM_adult_sim <- function(pars, fish_data, env_data = NULL, life_cycle = "SS", 
                          N_age, max_age, ages = NULL, SR_fun = "BH")
{
  # spawner-recruit functions
  SR <- function(alpha, Rmax, S, A, SR_fun) 
  {
    R <- switch(SR_fun,
                exp = alpha*S/A,
                BH = alpha*S/(A + alpha*S/Rmax),
                Ricker = alpha*(S/A)*exp(-alpha*S/(A*e()*Rmax)))
    return(R)
  }
  
  for(i in 1:length(pars)) assign(names(pars)[i], pars[i])
  for(j in 1:ncol(fish_data)) assign(colnames(fish_data[,i]), fish_data[,i])
  N <- nrow(fish_data)                  
  N_pop <- max(fish_data$pop)             
  A_pop <- tapply(A, pop, mean)
  pop <- as.numeric(factor(pop))
  year <- as.numeric(factor(year))
  adult_ages <- (max_age - N_age + 1):max_age   
  if(is.null(env_data)) env_data <- switch(life_cycle,
                                           SS = matrix(0, max(fish_data$year), 1),
                                           SMS = list(M = matrix(0, max(fish_data$year), 1),
                                                      MS = matrix(0, max(fish_data$year), 1)))
  # parameters
  Sigma_alphaRmax <- diag(c(sigma_alpha, sigma_Rmax)^2)
  Sigma_alphaRmax[1,2] <- rho_alphaRmax*sigma_alpha*sigma_Rmax
  Sigma_alphaRmax[2,1] <- Sigma_alphaRmax[1,2]
  aRmax <- exp(mvrnorm(N_pop, c(mu_alpha, mu_Rmax), Sigma_alphaRmax))
  alpha <- aRmax[,1]
  Rmax <- aRmax[,2]
  K <- switch(life_cycle, 
              SS = (alpha - 1)*Rmax/alpha,
              SMS = exp(mu_MS)*(alpha - 1)*Rmax/alpha)
  if(life_cycle == "SS")
  {
    phi <- rep(NA, max(year))
    phi[1] <- rnorm(1, 0, sigma_phi/sqrt(1 - rho_phi^2))
    for(i in 2:length(phi))
      phi[i] <- rnorm(1, rho_phi*phi[i-1], sigma_phi)
    phi <- phi + X %*% beta_phi
  }
  if(life_cycle = "SMS")  # only works for single pop at the moment
  {
    epsilon_M <- rep(NA, max(year))
    epsilon_MS <- rep(NA, max(year))
    epsilon_M[1] <- rnorm(1, 0, sigma_M/sqrt(1 - rho_M^2))
    epsilon_MS[1] <- rnorm(1, 0, sigma_MS/sqrt(1 - rho_MS^2))
    for(i in 2:max(year))
    {
      epsilon_M[i] <- rnorm(1, rho_M*epsilon_M[i-1], sigma_M)
      epsilon_MS[i] <- rnorm(1, rho_MS*epsilon_MS[i-1], sigma_MS) 
    }
  }
  mu_alr_p <- log(mu_p[1:(N_age-1)]) - log(mu_p[N_age])
  Sigma_gamma <-  diag(sigma_gamma^2) * (L_gamma %*% t(L_gamma))
  gamma <- mvrnorm(N_pop, mu_alr_p, Sigma_gamma)
  Sigma_alr_p <- diag(sigma_alr_p^2) * (L_alr_p %*% t(L_alr_p))
  alr_p <- t(apply(gamma[pop,], 1, function(x) mvrnorm(1, x, Sigma_alr_p)))
  e_alr_p <- exp(cbind(alr_p, 0))
  p <- sweep(e_alr_p, 1, rowSums(e_alr_p), "/")
  
  ###
  R_init <- data.frame(pop = rep(1:N_pop, each = max_age), year = NA, R = NA)
  for(i in 1:N_pop)
    R_init$year[R_init$pop==i] <- min(year[pop==i]) - (max_age:1)
  R_init$R <- rlnorm(nrow(R_init), log((A_pop*K)[R_init$pop]), sigma_proc)
  alr_p_init <- t(apply(gamma[R_init$pop,], 1, function(x) mvrnorm(1, x, Sigma_alr_p)))
  e_alr_p_init <- exp(cbind(alr_p_init, 0))
  p_init <- sweep(e_alr_p_init, 1, rowSums(e_alr_p_init), "/")
  R_init <- sweep(p_init, 1, R_init$R, "*")
  ###
  
  
  
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
      if(year[i] - adult_ages[j] < min(year[pop==pop[i]])) # initialize years 1:max_age
      {
        S_W_a[i,j] <- R_init[R_init$pop==pop[i] & R_init$year==year[i]-adult_ages[j], j]
      } else
      {
        S_W_a[i,j] <- R[i-adult_ages[j]]*p[i-adult_ages[j],j]
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
  dimnames(n_age_obs)[[2]] <- paste0("n_age", adult_ages, "_obs")
  n_H_obs <- rbinom(N, n_HW_obs, p_HOS)            # obs count of hatchery spawners
  n_W_obs <- n_HW_obs - n_H_obs                    # obs count of wild spawners
  
  return(list(sim_dat = data.frame(pop = pop, A = A, year = year, fit_p_HOS = p_HOS > 0,
                                   S_obs = S_obs, n_age_obs, 
                                   n_H_obs = n_H_obs, n_W_obs = n_W_obs, 
                                   B_take_obs = B_take, F_rate = F_rate),
              pars_out = c(pars, list(S_W_a = S_W_a, alpha = alpha, Rmax = Rmax, phi = phi, gamma = gamma, alr_p = alr_p,
                                      p_HOS = p_HOS, p = p, R_hat = R_hat, R = R))))
  
}
