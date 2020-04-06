#' Generate initial values for fitting either integrated or run-reconstruction
#' spawner-recruit models in Stan.
#'
#' @param data Named list of input data for fitting either an integrated or
#'   run-reconstruction spawner-recruit model in Stan, as returned by
#'   \code{stan_data}.
#' @param chains A positive integer specifying the number of Markov chains.
#' @param stan_model Character string giving the name of the Stan model being
#'   fit (".stan" filetype extension is not included).
#'
#' @importFrom stats aggregate na.omit
#'
#' @return A list with initial starting values for the parameters and states in
#'   the Stan model.
#'
#' @export
stan_init <- function(data, stan_model, chains) 
{
  if(stan_model %in% c("IPM_SS_np","IPM_SS_pp","IPM_SSpa_pp","IPM_SMS_np","IPM_ICchinook_pp"))
  {
    with(data, {
      N_pop <- max(pop)
      N_year <- max(year)
      S_obs_noNA <- S_obs
      S_obs[-which_S_obs] <- NA
      p_HOS_obs <- pmin(pmax(n_H_obs/(n_H_obs + n_W_obs), 0.01), 0.99)
      p_HOS_obs[n_H_obs + n_W_obs == 0] <- 0.5
      p_HOS_all <- rep(0,N)
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
      N_R_a_init <- N_pop*N_age*(min_age + (N_age - 1)/2) # total no. initial orphan recruit age classes
      
      # Maybe figure out a way to do this with a call to run_recon?
      for(i in 1:N)
        for(j in 1:N_age)
        {
          if(year[i] + adult_ages[j] <= max(year[pop==pop[i]]))
          {
            b <- ifelse(j==1, 0, B_rate_all[i+adult_ages[j]])
            f <- ifelse(j==1, 0, F_rate[i + adult_ages[j]])
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
      zeta_p <- apply(alr_p, 2, scale)
      mu_p <- aggregate(p, list(pop), mean)
      zeta_gamma <- aggregate(alr_p, list(pop), mean)[,-1, drop = FALSE]
      zeta_gamma <- apply(zeta_gamma, 2, scale)
      
      if(stan_model == "IPM_SMS_np") 
      {
        if(N_M_obs < N)
          M_obs[-which_M_obs] <- median(M_obs[which_M_obs])
        M0 <- c(M_obs[-(1:smolt_age)], rep(median(M_obs), smolt_age))
        s_MS <- pmin(S_obs_noNA/M0, 0.99)
      }
      
      if(stan_model == "IPM_SS_np") {
        return(lapply(1:chains, function(i)
          list(
            # recruitment
            alpha = array(exp(runif(N_pop, 1, 3)), dim = N_pop),
            Rmax = array(rlnorm(N_pop, log(tapply(R/A, pop, quantile, 0.9)), 0.5), dim = N_pop),
            beta = matrix(rnorm(N_X*N_pop, 0,1), N_pop, N_X),
            rho = array(runif(N_pop, 0.1, 0.7), dim = N_pop),
            sigma = array(runif(N_pop, 0.05, 2), dim = N_pop), 
            zeta_R = as.vector(scale(log(R)))*0.1,
            # spawner age structure
            mu_p = mu_p,
            sigma_p = matrix(runif(N_pop*(N_age-1), 0.5, 1), N_pop, N_age-1),
            zeta_p = zeta_p,
            # H/W composition, removals
            p_HOS = p_HOS_obs,
            B_rate = B_rate,
            # initial spawners, observation error
            S_init = rep(median(S_obs_noNA), max_age*N_pop),
            q_init = matrix(colMeans(q_obs), max_age*N_pop, N_age, byrow = T),
            tau = array(runif(N_pop, 0.5, 1), dim = N_pop)
          )
        ))
      } else if(stan_model %in% c("IPM_SS_pp","IPM_SSpa_pp")) {
        return(lapply(1:chains, function(i)
          list(
            # recruitment
            mu_alpha = runif(1, 1, 3),
            sigma_alpha = runif(1, 0.1, 0.5),
            zeta_alpha = array(runif(N_pop, -1, 1), dim = N_pop),
            mu_Rmax = rnorm(1, log(quantile(R/A,0.9)), 0.5),
            sigma_Rmax = runif(1, 0.1, 0.5),
            zeta_Rmax = array(runif(N_pop,-1,1), dim = N_pop),
            rho_alphaRmax = runif(1, -0.5, 0.5),
            beta_phi = array(rnorm(N_X, 0, 1), dim = N_X),
            rho_phi = runif(1, 0.1, 0.7),
            sigma_phi = runif(1, 0.1, 0.5),
            zeta_phi = array(rnorm(max(year,year_fwd), 0, 0.1), dim = max(year,year_fwd)),
            sigma = runif(1, 0.5, 1),
            zeta_R = as.vector(scale(log(R)))*0.1,
            # spawner age structure
            mu_p = colMeans(p), sigma_gamma = array(runif(N_age-1, 0.5, 1), dim = N_age-1),
            sigma_gamma = array(runif(N_age - 1, 0.5, 1), dim = N_age - 1),
            zeta_gamma = zeta_gamma,
            sigma_p = array(runif(N_age-1, 0.5, 1), dim = N_age-1),
            zeta_p = zeta_p,
            # H/W composition, removals
            p_HOS = p_HOS_obs,
            B_rate = B_rate,
            # initial spawners, observation error
            R_a_init = rep(median(S_obs_noNA/((1 - B_rate_all)*(1 - F_rate))), N_R_a_init),
            tau = runif(1, 0.5, 1)
          )
        ))
      } else if(stan_model == "IPM_SMS_np") {
        return(lapply(1:chains, function(i)
          list(
            # smolt recruitment
            alpha = array(exp(runif(N_pop,1,3)), dim = N_pop),
            Rmax = array(rlnorm(N_pop, log(tapply(R/A, pop, quantile, 0.9)), 0.5), dim = N_pop),
            beta_M = matrix(rnorm(N_X_M*N_pop,0,1), N_pop, N_X_M),
            rho_M = array(runif(N_pop, 0.1, 0.7), dim = N_pop),
            sigma_M = array(runif(N_pop, 0.05, 2), dim = N_pop), 
            zeta_M = as.vector(scale(log(M_obs)))*0.1,
            # SAR
            mu_MS = array(plogis(rnorm(N_pop, mean(qlogis(s_MS)), 0.5)), dim = N_pop),
            beta_MS = matrix(rnorm(N_X_MS*N_pop,0,1), N_pop, N_X_MS),
            rho_MS = array(runif(N_pop, 0.1, 0.7), dim = N_pop),
            sigma_MS = array(runif(N_pop, 0.05, 2), dim = N_pop), 
            zeta_MS = as.vector(scale(qlogis(s_MS))),
            # spawner age structure
            tau_M = array(runif(N_pop, 0.5, 1), dim = N_pop),
            tau_S = array(runif(N_pop, 0.5, 1), dim = N_pop),
            mu_p = mu_p,
            sigma_p = matrix(runif(N_pop*(N_age-1),0.5,1), N_pop, N_age-1),
            zeta_p = zeta_p,
            # H/W composition, removals
            p_HOS = p_HOS_obs,
            B_rate = B_rate,
            # initial states, observation error
            M_init = array(rep(median(M_obs), smolt_age*N_pop), dim = smolt_age*N_pop),
            S_init = rep(median(S_obs_noNA), max_age*N_pop),
            q_init = matrix(colMeans(q_obs), max_age*N_pop, N_age, byrow = T)
          )
        ))
      } else if(stan_model == "IPM_ICchinook_pp") {
        return(lapply(1:chains, function(i)
          list(
            # smolt recruitment
            mu_alpha = runif(1, 1, 3),
            sigma_alpha = runif(1, 0.1, 0.5),
            zeta_alpha = array(runif(N_pop, -1, 1), dim = N_pop),
            mu_Rmax = rnorm(1, log(quantile(R/A,0.9)), 0.5),
            sigma_Rmax = runif(1, 0.1, 0.5),
            zeta_Rmax = array(runif(N_pop, -1, 1), dim = N_pop),
            rho_alphaRmax = runif(1, -0.5, 0.5),
            beta_M = array(rnorm(N_X_M, 0, 1), dim = N_X_M),
            rho_M = runif(1, 0.1, 0.7),
            sigma_M = runif(1, 0.05, 2), 
            zeta_M = as.vector(scale(log(R)))*0.01,
            M_init = array(rep(median(S_obs_noNA)*100, smolt_age*N_pop), dim = smolt_age*N_pop),
            # downstream, SAR, upstream survival
            mu_D = qlogis(0.8),
            beta_D = array(rnorm(N_X_D, 0, 1), dim = N_X_D),
            rho_D = runif(1, 0.1, 0.7),
            sigma_D = runif(1, 0.05, 2),
            zeta_D = array(rnorm(max(year,year_fwd), 0, 0.1), dim = max(year,year_fwd)),
            mu_SAR = qlogis(0.01),
            beta_SAR = array(rnorm(N_X_SAR, 0, 1), dim = N_X_SAR),
            rho_SAR = runif(1, 0.1, 0.7),
            sigma_SAR = runif(1, 0.05, 2),
            zeta_SAR = array(rnorm(max(year,year_fwd), 0, 0.1), dim = max(year,year_fwd)),
            mu_U = qlogis(0.8),
            beta_U = array(rnorm(N_X_U, 0, 1), dim = N_X_U),
            rho_U = runif(1, 0.1, 0.7),
            sigma_U = runif(1, 0.05, 2),
            zeta_U = array(rnorm(max(year,year_fwd), 0, 0.1), dim = max(year,year_fwd)),
            # spawner age structure
            mu_p = colMeans(p), sigma_gamma = array(runif(N_age-1, 0.5, 1), dim = N_age-1),
            sigma_gamma = array(runif(N_age - 1, 0.5, 1), dim = N_age - 1),
            zeta_gamma = zeta_gamma,
            sigma_p = array(runif(N_age-1, 0.5, 1), dim = N_age-1),
            zeta_p = zeta_p,
            # H/W composition, removals
            p_HOS = p_HOS_obs,
            B_rate = B_rate,
            # initial spawners, observation error
            S_init = rep(median(S_obs_noNA), max_age*N_pop),
            q_init = matrix(colMeans(q_obs), max_age*N_pop, N_age, byrow = T),
            tau_S = runif(1, 0.5, 1)
          )
        ))
      }
    })
  } else if(stan_model == "IPM_SMaS_np") {
    with(data, {
      # This is a bit of a hack to avoid tedious age-structured run reconstuction
      N_pop <- max(pop)
      N_year <- max(year)
      N_GRage <- N_Mage*N_MSage
      max_age <- max_Mage + max_MSage
      q_M_obs <- sweep(n_Mage_obs, 1, rowSums(n_Mage_obs), "/")
      s_MS <- mean(pmin(S_obs/M_obs, 0.99), na.rm = TRUE)
      q_MS_obs <- sweep(n_MSage_obs, 1, rowSums(n_MSage_obs), "/")
      q_MS_obs_NA <- apply(is.na(q_MS_obs), 1, any)
      q_MS_obs[q_MS_obs_NA,] <- rep(colMeans(na.omit(q_MS_obs)), each = sum(q_MS_obs_NA))
      q_MS_pop <- as.matrix(aggregate(q_MS_obs, list(pop), mean))[,-1,drop=FALSE]
      mu_q_MS <- array(0, c(N_pop,N_Mage,N_MSage))
      for(i in 1:N_pop)
        for(j in 1:N_Mage)
          mu_q_MS[i,j,] <- as.vector(q_MS_pop[i,])
      q_GR_obs <- sweep(n_GRage_obs, 1, rowSums(n_GRage_obs), "/")
      p_HOS_obs <- pmin(pmax(n_H_obs/(n_H_obs + n_W_obs), 0.01), 0.99)
      p_HOS_obs[n_H_obs + n_W_obs == 0] <- 0.5
      p_HOS_all <- rep(0,N)
      p_HOS_all[which_H] <- p_HOS_obs
      S_W_obs <- S_obs*(1 - p_HOS_all)
      B_rate <- pmin(pmax(B_take_obs/(S_W_obs[which_B] + B_take_obs), 0.01), 0.99)
      B_rate[is.na(B_rate)] <- 0.1
      
      return(lapply(1:chains, function(i)
        list(
          # smolt recruitment
          alpha = array(rlnorm(N_pop, max(log(M_obs/S_obs), na.rm = TRUE), 1), dim = N_pop),
          Rmax = array(rlnorm(N_pop, log(tapply(M_obs/A, pop, quantile, 0.9, na.rm = TRUE)), 0.5), dim = N_pop),
          beta_M = matrix(rnorm(N_X_M*N_pop,0,1), N_pop, N_X_M),
          rho_M = array(runif(N_pop, 0.1, 0.7), dim = N_pop),
          sigma_M = array(runif(N_pop, 0.05, 2), dim = N_pop), 
          zeta_M = rnorm(N,0,0.1), 
          # smolt age structure
          mu_p_M = aggregate(q_M_obs, list(pop), mean, na.rm = TRUE),
          sigma_p_M = matrix(runif(N_pop*(N_Mage - 1), 0.05, 2), N_pop, N_Mage - 1),
          zeta_p_M = matrix(rnorm(N*(N_Mage - 1), 0, 0.1), N, N_Mage - 1),
          # SAR
          mu_MS = matrix(plogis(rnorm(N_pop*N_Mage, qlogis(s_MS), 0.5)), N_pop, N_Mage),
          beta_MS = matrix(rnorm(N_X_MS*N_pop,0,1), N_pop, N_X_MS),
          rho_MS = matrix(runif(N_pop, 0.1, 0.7), N_pop, N_Mage),
          sigma_MS = matrix(runif(N_pop, 0.05, 2), N_pop, N_Mage), 
          zeta_MS = matrix(rnorm(N*N_Mage, 0, 0.1), N, N_Mage),
          # ocean age structure
          mu_p_MS = mu_q_MS,
          sigma_p_MS = array(runif(N_pop*N_Mage*(N_MSage - 1), 0.05, 2), 
                             c(N_pop, N_Mage, N_MSage - 1)), 
          zeta_p_MS = matrix(rnorm(N*N_Mage*(N_MSage - 1), 0, 0.5), N, N_Mage*(N_MSage - 1)),
          # H/W composition, removals
          p_HOS = p_HOS_obs,
          B_rate = B_rate,
          # initial states, observation error
          M_init = rep(median(M_obs), max_Mage*N_pop),
          q_M_init = matrix(colMeans(q_M_obs, na.rm = TRUE), max_Mage*N_pop, N_Mage, byrow = T),
          S_init = rep(median(S_obs, na.rm = TRUE), N_pop*max_age),
          q_GR_init = matrix(colMeans(q_GR_obs, na.rm = TRUE), max_age*N_pop, N_GRage, byrow = T),
          tau_S = array(runif(N_pop, 0.01, 0.05), dim = N_pop),
          tau_M = array(runif(N_pop, 0.01, 0.05), dim = N_pop)
        )
      ))
    })
  } else if(stan_model %in% c("RR_SS_np","RR_SS_pp")) {
    with(data, {
      N_pop <- max(pop)
      N_year <- max(year)
      
      if(pool_pops)
      {
        # This is currently not based on the input data
        return(lapply(1:chains, function(i)
          list(mu_alpha = runif(1, 3, 6), 
               sigma_alpha = runif(1, 0.1, 0.5),
               zeta_alpha = array(runif(N_pop, -1, 1), dim = N_pop), 
               mu_Rmax = rnorm(1, log(quantile(S/A, 0.9, na.rm = T)), 0.5),
               sigma_Rmax = runif(1, 0.1, 0.5),
               zeta_Rmax = array(runif(N_pop, -1, 1), dim = N_pop), 
               rho_alphaRmax = runif(1, -0.5, 0.5),
               rho_phi = runif(1, 0.1, 0.7),
               sigma_phi = runif(1, 0.1, 0.5), 
               zeta_phi = array(rnorm(N_year, 0, 0.1), dim = N_year),
               sigma = runif(1, 0.1, 2))))
      } else {
        return(lapply(1:chains, function(i)
          list(alpha = array(exp(runif(N_pop, 1, 3)), dim = N_pop),
               Rmax = array(exp(runif(N_pop, -1, 0)), dim = N_pop),
               rho = array(runif(N_pop, 0.1, 0.7), dim = N_pop),
               sigma = array(runif(N_pop, 0.5, 1), dim = N_pop))))
      }
    })
  }
}
