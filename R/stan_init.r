#' Generate initial values for fitting either integrated or run-reconstruction
#' spawner-recruit models in Stan.
#'
#' @param data Named list of input data for fitting either an integrated or
#'   run-reconstruction spawner-recruit model in \code{stan}, as returned by
#'   \code{stan_data}.
#' @param chains A positive integer specifying the number of Markov chains.
#' @param stan_model Character string giving the name of the Stan model being
#'   fit (".stan" filetype extension is not included).
#'
#' @importFrom stats aggregate na.omit
#'
#' @return A list with initial starting values for all of the parameters and
#'   states in the Stan model.
#'
#' @export
stan_init <- function(data, stan_model, chains) 
{
  if(stan_model %in% c("IPM_SS_np","IPM_SS_pp","IPM_SS_F_pp","IPM_SMS_np"))
  {
    with(data, {
      S_obs_noNA <- S_obs
      if(N_S_obs < N)
        S_obs[-which_S_obs] <- NA
      p_HOS_obs <- pmin(pmax(n_H_obs/(n_H_obs + n_W_obs), 0.01), 0.99)
      p_HOS_obs[n_H_obs + n_W_obs == 0] <- 0.5
      p_HOS_all <- rep(0,N)
      if(N_H > 0)
        p_HOS_all[which_H] <- p_HOS_obs
      min_age <- max_age - N_age + 1
      adult_ages <- min_age:max_age
      q_obs <- sweep(n_age_obs, 1, rowSums(n_age_obs), "/")
      q_obs_NA <- apply(is.na(q_obs), 1, any)
      q_obs[q_obs_NA,] <- rep(colMeans(na.omit(q_obs)), each = sum(q_obs_NA))
      R_a <- matrix(NA, N, N_age)
      S_W_obs <- S_obs*(1 - p_HOS_all)
      B_rate_all <- rep(0,N)
      B_rate <- pmin(pmax(B_take_obs/(S_W_obs[which_B]*(1 - q_obs[which_B,1]) + B_take_obs), 0.01), 0.99)
      B_rate[is.na(B_rate)] <- 0.1
      B_rate_all[which_B] <- B_rate
      year <- as.numeric(factor(year))
      if(stan_model == "IPM_SS_F_pp")
        F_rate_obs <- pmin(pmax(catch_data$C_obs/catch_data$R_F_obs, 0.01), 0.99)
      
      # Maybe figure out a way to do this with a call to run_recon?
      for(i in 1:N)
        for(j in 1:N_age)
        {
          if(year[i] + adult_ages[j] <= max(year[pop==pop[i]]))
          {
            b <- ifelse(j==1, 0, B_rate_all[i+adult_ages[j]])
            f <- ifelse(stan_model == "IPM_SS_F_pp",
                        F_rate_obs[year[i] + adult_ages[j]],
                        ifelse(j==1, 0, F_rate[i + adult_ages[j]]))
            R_a[i,j] <- S_W_obs[i + adult_ages[j]]*q_obs[i + adult_ages[j],j]/((1 - b)*(1 - f))
          }
        }
      
      R_a <- pmax(R_a, min(1, R_a[R_a > 0], na.rm = T))
      R <- rowSums(R_a)
      R[is.na(R)] <- max(R, na.rm = T)
      p <- sweep(R_a, 1, R, "/")
      p_NA <- apply(is.na(p), 1, any)
      p[p_NA, ] <- rep(colMeans(na.omit(p)), each = sum(p_NA))
      alr_p <- sweep(log(p[, 1:(N_age-1), drop = FALSE]), 1, log(p[,N_age]), "-")
      epsilon_p_z <- apply(alr_p, 2, scale)
      mu_p <- aggregate(p, list(pop), mean)
      epsilon_gamma_z <- aggregate(alr_p, list(pop), mean)[,-1, drop = FALSE]
      epsilon_gamma_z <- apply(epsilon_gamma_z, 2, scale)
      
      if(stan_model == "IPM_SMS_np") 
      {
        if(N_M_obs < N)
          M_obs[-which_M_obs] <- median(M_obs[which_M_obs])
        M0 <- c(M_obs[-(1:smolt_age)], rep(median(M_obs), smolt_age))
        s_MS <- R/M0
      }
      
      if(stan_model == "IPM_SS_np") {
        return(lapply(1:chains, function(i)
          list(alpha = array(exp(runif(max(pop), 1, 3)), dim = max(pop)),
               Rmax = array(rlnorm(max(pop), log(tapply(R/A, pop, quantile, 0.9)), 0.5), dim = max(pop)),
               beta = matrix(rnorm(N_X*max(pop), 0,1), max(pop), N_X),
               rho = array(runif(max(pop), 0.1, 0.7), dim = max(pop)),
               sigma = array(runif(max(pop), 0.05, 2), dim = max(pop)), 
               tau = array(runif(max(pop), 0.5, 1), dim = max(pop)),
               mu_p = mu_p,
               sigma_p = matrix(runif(max(pop)*(N_age-1), 0.5, 1), max(pop), N_age-1),
               epsilon_p_z = epsilon_p_z,
               S_init = rep(median(S_obs_noNA), max_age*max(pop)),
               q_init = matrix(colMeans(q_obs), max_age*max(pop), N_age, byrow = T),
               p_HOS = p_HOS_obs,
               epsilon_R_z = as.vector(scale(log(R)))*0.1,
               B_rate = B_rate)))
      } else if(stan_model == "IPM_SS_pp") {
        return(lapply(1:chains, function(i)
          list(mu_alpha = runif(1, 1, 3),
               sigma_alpha = runif(1, 0.1, 0.5),
               epsilon_alpha_z = array(runif(max(pop), -1, 1), dim = max(pop)),
               mu_Rmax = rnorm(1, log(quantile(R/A,0.9)), 0.5),
               sigma_Rmax = runif(1, 0.1, 0.5),
               epsilon_Rmax_z = array(runif(max(pop),-1,1), dim = max(pop)),
               rho_alphaRmax = runif(1, -0.5, 0.5),
               beta_phi = array(rnorm(N_X, 0, 1), dim = N_X),
               rho_phi = runif(1, 0.1, 0.7),
               sigma_phi = runif(1, 0.1, 0.5),
               epsilon_phi_z = array(rnorm(max(year), 0, 0.1), dim = max(year,year_fwd)),
               sigma = runif(1, 0.5, 1),
               tau = runif(1, 0.5, 1),
               mu_p = colMeans(p), sigma_gamma = array(runif(N_age-1, 0.5, 1), dim = N_age-1),
               epsilon_gamma_z = epsilon_gamma_z,
               sigma_p = array(runif(N_age-1, 0.5, 1), dim = N_age-1),
               epsilon_p_z = epsilon_p_z,
               S_init = rep(median(S_obs_noNA), max_age*max(pop)),
               q_init = matrix(colMeans(q_obs), max_age*max(pop), N_age, byrow = T),
               p_HOS = p_HOS_obs,
               epsilon_R_z = as.vector(scale(log(R)))*0.1,
               B_rate = B_rate)))
      } else if(stan_model == "IPM_SS_F_pp") {
        return(lapply(1:chains, function(i)
          list(mu_alpha = runif(1, 1, 3),
               sigma_alpha = runif(1, 0.1, 0.5),
               epsilon_alpha_z = array(runif(max(pop), -1, 1), dim = max(pop)),
               mu_Rmax = rnorm(1, log(quantile(R/A,0.9)), 0.5),
               sigma_Rmax = runif(1, 0.1, 0.5),
               epsilon_Rmax_z = array(runif(max(pop), -1, 1), dim = max(pop)),
               rho_alphaRmax = runif(1, -0.5, 0.5),
               beta_phi = array(rnorm(N_X, 0, 1), dim = N_X),
               rho_phi = runif(1, 0.1, 0.7),
               sigma_phi = runif(1, 0.1, 0.5),
               epsilon_phi_z = array(rnorm(max(year), 0, 0.1), dim = max(year)),
               sigma = runif(1, 0.5, 1),
               tau = runif(1, 0.5, 1),
               mu_p = colMeans(p), sigma_gamma = array(runif(N_age-1, 0.5, 1), dim = N_age-1),
               epsilon_gamma_z = epsilon_gamma_z,
               sigma_p = array(runif(N_age-1, 0.5, 1), dim = N_age-1),
               epsilon_p_z = epsilon_p_z,
               S_init = rep(median(S_obs_noNA), max_age*max(pop)),
               q_init = matrix(colMeans(q_obs), max_age*max(pop), N_age, byrow = T),
               p_HOS = p_HOS_obs,
               epsilon_R_z = as.vector(scale(log(R)))*0.1,
               c1 = rnorm(1, 0, 0.5), c2 = 0,
               sigma_log_C = runif(1, 0.1, 0.5),
               B_rate = B_rate)))
      } else if(stan_model == "IPM_SMS_np") {
        return(lapply(1:chains, function(i)
          list(alpha = array(exp(runif(max(pop),1,3)), dim = max(pop)),
               Rmax = array(rlnorm(max(pop), log(tapply(R/A, pop, quantile, 0.9)), 0.5), dim = max(pop)),
               beta_M = matrix(rnorm(N_X_M*max(pop),0,1), max(pop), N_X_M),
               rho_M = array(runif(max(pop), 0.1, 0.7), dim = max(pop)),
               sigma_M = array(runif(max(pop), 0.05, 2), dim = max(pop)), 
               mu_MS = plogis(rnorm(1, mean(qlogis(s_MS)), 0.5)),
               beta_MS = matrix(rnorm(N_X_MS*max(pop),0,1), max(pop), N_X_MS),
               rho_MS = array(runif(max(pop), 0.1, 0.7), dim = max(pop)),
               sigma_MS = array(runif(max(pop), 0.05, 2), dim = max(pop)), 
               tau_M = array(runif(max(pop), 0.5, 1), dim = max(pop)),
               tau_S = array(runif(max(pop), 0.5, 1), dim = max(pop)),
               mu_p = mu_p,
               sigma_p = matrix(runif(max(pop)*(N_age-1),0.5,1), max(pop), N_age-1),
               epsilon_p_z = epsilon_p_z,
               M_init = rep(median(M_obs), smolt_age*max(pop)),
               S_init = rep(median(S_obs_noNA), max_age*max(pop)),
               q_init = matrix(colMeans(q_obs), max_age*max(pop), N_age, byrow = T),
               p_HOS = p_HOS_obs,
               epsilon_M_z = as.vector(scale(log(M_obs)))*0.1,
               epsilon_MS_z = as.vector(scale(qlogis(s_MS))),
               B_rate = B_rate)))
      }
    })
  } else if(stan_model %in% c("RR_SS_np","RR_SS_pp")) {
    with(data, {
      if(pool_pops)
      {
        # This is currently not based on the input data
        return(lapply(1:chains, function(i)
          list(mu_alpha = runif(1, 3, 6), 
               sigma_alpha = runif(1, 0.1, 0.5),
               epsilon_alpha_z = array(runif(max(pop), -1, 1), dim = max(pop)), 
               mu_Rmax = rnorm(1, log(quantile(S/A, 0.9, na.rm = T)), 0.5),
               sigma_Rmax = runif(1, 0.1, 0.5),
               epsilon_Rmax_z = array(runif(max(pop), -1, 1), dim = max(pop)), 
               rho_alphaRmax = runif(1, -0.5, 0.5),
               rho_phi = runif(1, 0.1, 0.7),
               sigma_phi = runif(1, 0.1, 0.5), 
               epsilon_phi_z = array(rnorm(max(year), 0, 0.1), dim = max(year)),
               sigma = runif(1, 0.1, 2))))
      } else {
        return(lapply(1:chains, function(i)
          list(alpha = array(exp(runif(max(pop), 1, 3)), dim = max(pop)),
               Rmax = array(exp(runif(max(pop), -1, 0)), dim = max(pop)),
               rho = array(runif(max(pop), 0.1, 0.7), dim = max(pop)),
               sigma = array(runif(max(pop), 0.5, 1), dim = max(pop)))))
      }
    })
  }
}
