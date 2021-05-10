#' Assemble input data for integrated or run-reconstruction spawner-recruit
#' models.
#'
#' @inheritParams salmonIPM
#'
#' @return A named list that can be passed to `stan` as the `data`
#'   argument.
#'
#' @export

stan_data <- function(stan_model, SR_fun = "BH", par_models = NULL, scale = TRUE, 
                      ages = NULL, age_S_obs = NULL, age_S_eff = NULL, conditionGRonMS = FALSE,
                      fish_data, fish_data_fwd = NULL, fecundity_data = NULL, prior_data = NULL)
{
  if(!stan_model %in% c("IPM_SS_np","IPM_SS_pp","IPM_SMS_np","IPM_SMS_pp",
                        "IPM_SMaS_np","IPM_LCRchum_pp","IPM_ICchinook_pp",
                        "RR_SS_np","RR_SS_pp"))
    stop("Stan model ", stan_model, " does not exist")
  
  fish_data <- as.data.frame(fish_data)
  life_cycle <- strsplit(stan_model, "_")[[1]][2]
  model_life_cycle <- paste(strsplit(stan_model, "_")[[1]][1], life_cycle, sep = "_")
  for(i in names(fish_data)) assign(i, fish_data[[i]])
  N <- nrow(fish_data)
  pop <- as.numeric(factor(pop))
  year <- as.numeric(factor(year))
  n_age_obs <- as.matrix(fish_data[, grep("n_age", names(fish_data))])
  n_age_obs <- replace(n_age_obs, is.na(n_age_obs), 0)
  fit_p_HOS <- as.logical(fit_p_HOS)
  n_W_obs = as.vector(replace(n_W_obs[fit_p_HOS], is.na(n_W_obs[fit_p_HOS]), 0))
  n_H_obs = as.vector(replace(n_H_obs[fit_p_HOS], is.na(n_H_obs[fit_p_HOS]), 0))
  
  if(life_cycle == "SMaS") {
    if(!is.logical(conditionGRonMS))
      stop("conditionGRonMS must be TRUE or FALSE for model 'IPM_SMaS_np'")
    max_Mage <- max(as.numeric(substring(names(fish_data)[grep("n_Mage", names(fish_data))], 7, 7)))
    max_MSage <- max(as.numeric(substring(names(fish_data)[grep("n_MSage", names(fish_data))], 8, 8))) 
    max_age <- max_Mage + max_MSage
  } else {
    max_age <- max(as.numeric(substring(names(fish_data)[grep("n_age", names(fish_data))], 6, 6)))
  }
  
  if(any(unlist(tapply(year, pop, diff)) != 1))
    stop("Non-consecutive years not allowed in fish_data")
  
  for(i in c("pop","year","A","fit_p_HOS","B_take_obs"))
    if(any(is.na(fish_data[,i])))
      stop("Missing values not allowed in fish_data$", i)
  
  if(any(is.na(fish_data_fwd)))
    stop("Missing values not allowed in fish_data_fwd")
  
  if(stan_model == "IPM_SS_pp" & is.null(age_S_obs))
    age_S_obs <- rep(1, sum(grepl("n_age", names(fish_data))))
  age_S_obs <- as.numeric(age_S_obs)
  
  if(stan_model == "IPM_SS_pp" & is.null(age_S_eff))
    age_S_eff <- rep(1, sum(grepl("n_age", names(fish_data))))
  age_S_eff <- as.numeric(age_S_eff)
  
  if(life_cycle != "SS" & (any(is.na(ages)) | is.null(ages)))
    stop("Multi-stage models must specify age in years for all stages.")
  
  age_NA_check <- is.na(fish_data[,grep("n_age", names(fish_data))])
  if(any(!rowSums(age_NA_check) %in% c(0, nrow(age_NA_check))))
    stop("Conflicting NAs in age frequency data in rows ", 
         which(!rowSums(age_NA_check) %in% c(0, nrow(age_NA_check))))
  
  if(any(is.na(n_W_obs) != is.na(n_H_obs)))
    stop("Conflicting NAs in n_W_obs and n_H_obs in rows ", 
         which(is.na(n_W_obs) != is.na(n_H_obs)))
  
  if(stan_model == "IPM_LCRchum_pp") {
    if(any(is.na(n_M_obs) != is.na(n_F_obs)))
      stop("Conflicting NAs in n_M_obs and n_F_obs in rows ", 
           which(is.na(n_M_obs) != is.na(n_F_obs)))
    if(is.null(fecundity_data)) {
      stop("Fecundity data must be provided for Lower Columbia chum IPM") 
    } else if(any(is.na(fecundity_data))) {
      stop("Missing values are not allowed in fecundity data")
    }
  }
  
  F_rate_check <- tapply(F_rate, pop, function(x) any(is.na(x[-c(1:max_age)])))
  if(any(F_rate_check))
    stop("Missing values not allowed in fish_data$F_rate except in years 1:max_age", i)
  
  X <- par_model_matrix(par_models = par_models, scale = scale, fish_data = fish_data)
  
  if(is.null(fish_data_fwd)) {
    N_fwd <- 0
    fish_data_fwd <- data.frame(pop = 1, year = 1, A = 0, F_rate = 0, B_rate = 0, p_HOS = 0)
    fish_data_fwd <- fish_data_fwd[c(),]
  } else {
    if(stan_model != "IPM_SS_pp")
      warning("Argument fish_data_fwd is ignored unless stan_model == 'IPM_SS_pp'")
    N_fwd <- nrow(fish_data_fwd)
    fish_data_fwd <- as.data.frame(fish_data_fwd)
    fish_data_fwd$pop <- factor(fish_data_fwd$pop, levels = levels(factor(fish_data$pop)))
    if(any(!fish_data_fwd$pop %in% fish_data$pop))
      stop("All populations in fish_data_fwd must appear in fish_data")
    year_check <- tapply(fish_data$year, fish_data$pop, max)
    year_check <- year_check[names(year_check) %in% fish_data_fwd$pop]
    year_fwd_check <- tapply(fish_data_fwd$year, fish_data_fwd$pop, min)
    year_fwd_check <- year_fwd_check[names(year_fwd_check) %in% fish_data_fwd$pop]
    if(any(year_fwd_check != year_check + 1))
      stop("First year in fish_data_fwd must equal 1 + last year in fish_data for each population")
    fish_data_fwd$pop <- as.numeric(fish_data_fwd$pop)
    fish_data_fwd$year <- as.numeric(factor(fish_data_fwd$year, 
                                            levels = levels(factor(c(fish_data$year, fish_data_fwd$year)))))
  }
  
  if(stan_model == "IPM_ICchinook_pp") {
    if(is.null(prior_data)) {
      stop("Priors for survival must be specified for model IPM_ICchinook_pp")
    } else {
      s <- prior_data$s[prior_data$s$year %in% fish_data$year,]
      
      prior_D <- na.omit(s[,c("year","mu_prior_D","sigma_prior_D")])
      N_prior_D <- nrow(prior_D)
      which_prior_D <- match(prior_D$year, sort(unique(fish_data$year)))
      mu_prior_D <- prior_D$mu_prior_D
      sigma_prior_D <- prior_D$sigma_prior_D
      
      prior_SAR <- na.omit(s[,c("year","mu_prior_SAR","sigma_prior_SAR")])
      N_prior_SAR <- nrow(prior_SAR)
      which_prior_SAR <- match(prior_SAR$year, sort(unique(fish_data$year)))
      mu_prior_SAR <- prior_SAR$mu_prior_SAR
      sigma_prior_SAR <- prior_SAR$sigma_prior_SAR
      
      prior_U <- na.omit(s[,c("year","mu_prior_U","sigma_prior_U")])
      N_prior_U <- nrow(prior_U)
      which_prior_U <- match(prior_U$year, sort(unique(fish_data$year)))
      mu_prior_U <- prior_U$mu_prior_U
      sigma_prior_U <- prior_U$sigma_prior_U
    }
  }
  
  if(model_life_cycle == "RR_SS") {
    rr_dat <- run_recon(fish_data)
    which_fit <- which(!is.na(rr_dat$R_obs) & !is.na(rr_dat$S_obs))
    N_fit <- length(which_fit)
    R_obs <- rr_dat$R_obs
    p_pop_obs <- aggregate(rr_dat[,grep("p_age", names(rr_dat))], list(pop = rr_dat$pop), mean, na.rm = TRUE)[,-1]
  }
  
  out <- switch(model_life_cycle,
                IPM_SS = list(
                    # info for observed data
                    N = N,
                    pop = pop, 
                    year = year,
                    # info for forward simulations
                    N_fwd = N_fwd,
                    pop_fwd = as.vector(fish_data_fwd$pop),
                    year_fwd = as.vector(fish_data_fwd$year),
                    A_fwd = as.vector(fish_data_fwd$A),
                    B_rate_fwd = as.vector(fish_data_fwd$B_rate),
                    F_rate_fwd = as.vector(fish_data_fwd$F_rate),
                    p_HOS_fwd = as.vector(fish_data_fwd$p_HOS),
                    # recruitment
                    SR_fun = switch(SR_fun, exp = 1, BH = 2, Ricker = 3),
                    A = A,
                    K_alpha = ifelse(is.null(X[["alpha"]]), 0, ncol(X[["alpha"]])), 
                    X_alpha = if(is.null(X[["alpha"]])) matrix(0,N,0) else X[["alpha"]],
                    K_Rmax = ifelse(is.null(X[["Rmax"]]), 0, ncol(X[["Rmax"]])), 
                    X_Rmax = if(is.null(X[["Rmax"]])) matrix(0,N,0) else X[["Rmax"]],
                    K_R = ifelse(is.null(X[["R"]]), 0, ncol(X[["R"]])), 
                    X_R = if(is.null(X[["R"]])) matrix(0,N,0) else X[["R"]],
                    # fishery and hatchery removals
                    F_rate = replace(F_rate, is.na(F_rate), 0),
                    N_B = sum(B_take_obs > 0),
                    which_B = as.vector(which(B_take_obs > 0)),
                    B_take_obs = B_take_obs[B_take_obs > 0],
                    # spawner abundance
                    N_S_obs = sum(!is.na(S_obs)),
                    which_S_obs = as.vector(which(!is.na(S_obs))),
                    S_obs = replace(S_obs, is.na(S_obs) | S_obs==0, 1),
                    # spawner age structure
                    N_age = sum(grepl("n_age", names(fish_data))), 
                    max_age = max_age,
                    n_age_obs = n_age_obs,
                    age_S_obs = as.vector(age_S_obs),
                    age_S_eff = as.vector(age_S_eff),
                    # H/W composition
                    N_H = sum(fit_p_HOS),
                    which_H = as.vector(which(fit_p_HOS)),
                    n_W_obs = n_W_obs,
                    n_H_obs = n_H_obs
                  ),
                
                IPM_SMS = list(
                    # info for observed data
                    N = N,
                    pop = pop, 
                    year = year,
                    # smolt production
                    SR_fun = switch(SR_fun, exp = 1, BH = 2, Ricker = 3),
                    A = A,
                    K_alpha = ifelse(is.null(X[["alpha"]]), 0, ncol(X[["alpha"]])), 
                    X_alpha = if(is.null(X[["alpha"]])) matrix(0,N,0) else X[["alpha"]],
                    K_Mmax = ifelse(is.null(X[["Mmax"]]), 0, ncol(X[["Mmax"]])), 
                    X_Mmax = if(is.null(X[["Mmax"]])) matrix(0,N,0) else X[["Mmax"]],
                    K_M = ifelse(is.null(X[["M"]]), 0, ncol(X[["M"]])), 
                    X_M = if(is.null(X[["M"]])) matrix(0,N,0) else X[["M"]],
                    # smolt abundance
                    N_M_obs = sum(!is.na(M_obs)),
                    which_M_obs = as.vector(which(!is.na(M_obs))),
                    M_obs = replace(M_obs, is.na(M_obs) | M_obs==0, 1),
                    N_age = sum(grepl("n_age", names(fish_data))), 
                    smolt_age = ages$M,
                    # SAR
                    K_MS = ifelse(is.null(X[["s_MS"]]), 0, ncol(X[["s_MS"]])), 
                    X_MS = if(is.null(X[["s_MS"]])) matrix(0,N,0) else X[["s_MS"]],
                    # fishery and hatchery removals
                    F_rate = replace(F_rate, is.na(F_rate), 0),
                    N_B = sum(B_take_obs > 0),
                    which_B = as.vector(which(B_take_obs > 0)),
                    B_take_obs = B_take_obs[B_take_obs > 0],
                    # spawner abundance
                    N_S_obs = sum(!is.na(S_obs)),
                    which_S_obs = as.vector(which(!is.na(S_obs))),
                    S_obs = replace(S_obs, is.na(S_obs) | S_obs==0, 1),
                    # spawner age structure
                    max_age = max_age,
                    n_age_obs = n_age_obs,
                    # H/W composition
                    N_H = sum(fit_p_HOS),
                    which_H = as.vector(which(fit_p_HOS)),
                    n_W_obs = n_W_obs,
                    n_H_obs = n_H_obs
                ),
                
                IPM_SMaS = list(
                    # info for observed data
                    N = N,
                    pop = pop, 
                    year = year,
                    # smolt production
                    SR_fun = switch(SR_fun, exp = 1, BH = 2, Ricker = 3),
                    A = A,
                    K_M = ifelse(is.null(X[["M"]]), 0, ncol(X[["M"]])), 
                    X_M = if(is.null(X[["M"]])) matrix(0,N,0) else X[["M"]],
                    # smolt abundance
                    N_M_obs = sum(!is.na(M_obs)),
                    which_M_obs = as.vector(which(!is.na(M_obs))),
                    M_obs = replace(M_obs, is.na(M_obs) | M_obs==0, 1),
                    # smolt age structure
                    N_Mage = sum(grepl("n_Mage", names(fish_data))),
                    max_Mage = max_Mage,
                    n_Mage_obs = as.matrix(fish_data[,grep("n_Mage", names(fish_data))]),
                    # SAR
                    K_MS = ifelse(is.null(X[["s_MS"]]), 0, ncol(X[["s_MS"]])), 
                    X_MS = if(is.null(X[["s_MS"]])) matrix(0,N,0) else X[["s_MS"]],
                    # fishery and hatchery removals
                    F_rate = replace(F_rate, is.na(F_rate), 0),
                    N_B = sum(B_take_obs > 0),
                    which_B = as.vector(which(B_take_obs > 0)),
                    B_take_obs = B_take_obs[B_take_obs > 0],
                    # spawner abundance
                    N_S_obs = sum(!is.na(S_obs)),
                    which_S_obs = as.vector(which(!is.na(S_obs))),
                    S_obs = replace(S_obs, is.na(S_obs) | S_obs==0, 1),
                    # spawner ocean and Gilbert-Rich age structure
                    N_MSage = sum(grepl("n_MSage", names(fish_data))), 
                    max_MSage = max_MSage,
                    n_MSage_obs = as.matrix(fish_data[,grep("n_MSage", names(fish_data))]),
                    n_GRage_obs = as.matrix(fish_data[,grep("n_GRage", names(fish_data))]),
                    conditionGRonMS = as.numeric(conditionGRonMS),
                    # H/W composition
                    N_H = sum(fit_p_HOS),
                    which_H = as.vector(which(fit_p_HOS)),
                    n_W_obs = n_W_obs,
                    n_H_obs = n_H_obs
                ),
                
                IPM_LCRchum = list(
                    # info for observed data
                    N = N,
                    pop = pop, 
                    year = year,
                    # egg deposition
                    SR_fun = switch(SR_fun, exp = 1, BH = 2, Ricker = 3),
                    A = A,
                    N_E = nrow(fecundity_data),
                    age_E = fecundity_data$age_E,
                    E_obs = fecundity_data$E_obs,
                    # spawner-smolt productivity
                    smolt_age = ages$M,
                    K_psi = ifelse(is.null(X[["psi"]]), 0, ncol(X[["psi"]])), 
                    X_psi = if(is.null(X[["psi"]])) matrix(0,N,0) else X[["psi"]],
                    K_Mmax = ifelse(is.null(X[["Mmax"]]), 0, ncol(X[["Mmax"]])), 
                    X_Mmax = if(is.null(X[["Mmax"]])) matrix(0,N,0) else X[["Mmax"]],
                    K_M = ifelse(is.null(X[["M"]]), 0, ncol(X[["M"]])), 
                    X_M = if(is.null(X[["M"]])) matrix(0,N,0) else X[["M"]],
                    # smolt abundance and observation error
                    N_M_obs = sum(!is.na(M_obs)),
                    which_M_obs = as.vector(which(!is.na(M_obs))),
                    M_obs = replace(M_obs, is.na(M_obs) | M_obs==0, 1),
                    N_age = sum(grepl("n_age", names(fish_data))), 
                    N_tau_M_obs = sum(!is.na(tau_M_obs)),
                    which_tau_M_obs = as.vector(which(!is.na(tau_M_obs))),
                    tau_M_obs = replace(tau_M_obs, is.na(tau_M_obs), 0),
                    N_upstream = sum(!is.na(downstream_trap)),
                    which_upstream = as.vector(which(!is.na(downstream_trap))),
                    downstream_trap = na.omit(downstream_trap),
                    # SAR
                    K_MS = ifelse(is.null(X[["s_MS"]]), 0, ncol(X[["s_MS"]])), 
                    X_MS = if(is.null(X[["s_MS"]])) matrix(0,N,0) else X[["s_MS"]],
                    # fishery and hatchery removals and translocations
                    F_rate = replace(F_rate, is.na(F_rate), 0),
                    N_B = sum(B_take_obs > 0),
                    which_B = as.vector(which(B_take_obs > 0)),
                    B_take_obs = B_take_obs[B_take_obs > 0],
                    S_add_obs = replace(S_add_obs, is.na(S_add_obs), 0),
                    # spawner abundance and observation error
                    N_S_obs = sum(!is.na(S_obs)),
                    which_S_obs = as.vector(which(!is.na(S_obs))),
                    S_obs = replace(S_obs, is.na(S_obs) | S_obs==0, 1),
                    N_tau_S_obs = sum(!is.na(tau_S_obs)),
                    which_tau_S_obs = as.vector(which(!is.na(tau_S_obs))),
                    tau_S_obs = replace(tau_S_obs, is.na(tau_S_obs), 0),
                    # spawner age structure and sex ratio
                    max_age = max_age,
                    n_age_obs = n_age_obs,
                    n_M_obs = as.vector(n_M_obs),
                    n_F_obs = as.vector(n_F_obs),
                    p_G_obs = as.vector(p_G_obs),
                    # H/W composition
                    N_H = sum(fit_p_HOS),
                    which_H = as.vector(which(fit_p_HOS)),
                    n_W_obs = n_W_obs,
                    n_H_obs = n_H_obs
                ),
                
                IPM_ICchinook = list(
                    # info for observed data
                    N = N,
                    pop = pop, 
                    year = year,
                    # info for forward simulations
                    N_fwd = N_fwd,
                    pop_fwd = as.vector(fish_data_fwd$pop),
                    year_fwd = as.vector(fish_data_fwd$year),
                    A_fwd = as.vector(fish_data_fwd$A),
                    B_rate_fwd = as.vector(fish_data_fwd$B_rate),
                    F_rate_fwd = as.vector(fish_data_fwd$F_rate),
                    p_HOS_fwd = as.vector(fish_data_fwd$p_HOS),
                    # recruitment
                    SR_fun = switch(SR_fun, exp = 1, BH = 2, Ricker = 3),
                    smolt_age = ages$M,
                    A = A,
                    K_alpha = ifelse(is.null(X[["alpha"]]), 0, ncol(X[["alpha"]])), 
                    X_alpha = if(is.null(X[["alpha"]])) matrix(0,N,0) else X[["alpha"]],
                    K_Mmax = ifelse(is.null(X[["Mmax"]]), 0, ncol(X[["Mmax"]])), 
                    X_Mmax = if(is.null(X[["Mmax"]])) matrix(0,N,0) else X[["Mmax"]],
                    K_M = ifelse(is.null(X[["M"]]), 0, ncol(X[["M"]])), 
                    X_M = if(is.null(X[["M"]])) matrix(0,N,0) else X[["M"]],
                    # downstream, SAR, upstream survival
                    # K_D = ncol(env_data$s_D),
                    # X_D = as.matrix(env_data$s_D),
                    N_prior_D = N_prior_D,
                    which_prior_D = as.vector(which_prior_D),
                    mu_prior_D = as.vector(mu_prior_D),
                    sigma_prior_D = as.vector(sigma_prior_D),
                    # K_SAR = ncol(env_data$s_SAR),
                    # X_SAR = as.matrix(env_data$s_SAR),
                    N_prior_SAR = N_prior_SAR,
                    which_prior_SAR = as.vector(which_prior_SAR),
                    mu_prior_SAR = as.vector(mu_prior_SAR),
                    sigma_prior_SAR = as.vector(sigma_prior_SAR),
                    # K_U = ncol(env_data$s_U),
                    # X_U = as.matrix(env_data$s_U),
                    N_prior_U = N_prior_U,
                    which_prior_U = as.vector(which_prior_U),
                    mu_prior_U = as.vector(mu_prior_U),
                    sigma_prior_U = as.vector(sigma_prior_U),
                    # fishery and hatchery removals
                    F_rate = replace(F_rate, is.na(F_rate), 0),
                    N_B = sum(B_take_obs > 0),
                    which_B = as.vector(which(B_take_obs > 0)),
                    B_take_obs = B_take_obs[B_take_obs > 0],
                    # spawner abundance
                    N_S_obs = sum(!is.na(S_obs)),
                    which_S_obs = as.vector(which(!is.na(S_obs))),
                    S_obs = replace(S_obs, is.na(S_obs) | S_obs==0, 1),
                    # spawner age structure
                    N_age = sum(grepl("n_age", names(fish_data))), 
                    max_age = max_age,
                    n_age_obs = n_age_obs,
                    # H/W composition
                    N_H = sum(fit_p_HOS),
                    which_H = as.vector(which(fit_p_HOS)),
                    n_W_obs = n_W_obs,
                    n_H_obs = n_H_obs
                ),
                
                RR_SS = list(
                    SR_fun = switch(SR_fun, exp = 1, BH = 2, Ricker = 3),
                    N = length(S_obs),
                    pop = pop, 
                    year = year,
                    N_fit = N_fit,
                    which_fit = as.vector(which_fit),
                    S_obs = replace(S_obs, S_obs == 0 | is.na(S_obs), 1),
                    R_obs = replace(R_obs, R_obs == 0 | is.na(R_obs), 1),
                    A = A,
                    S_NA = as.vector(as.integer(is.na(S_obs))),
                    R_NA = as.vector(as.integer(is.na(R_obs))),
                    N_age = sum(grepl("p_age", names(rr_dat))),
                    max_age = max(as.numeric(substring(names(rr_dat)[grep("p_age", names(rr_dat))], 6, 6))),
                    p_pop_obs = as.matrix(p_pop_obs)
                  )
  )  # end switch()
  
  return(out)
}
