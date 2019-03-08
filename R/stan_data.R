#' Assemble input data for integrated or run-reconstruction spawner-recruit
#' models.
#'
#' @param fish_data Data frame that includes the following \code{colnames}, in
#'   no particular order except where noted: \describe{
#'   \item{\code{pop}}{Numeric or character population ID.}
#'   \item{\code{year}}{Numeric variable giving the year the fish spawned (i.e.,
#'   the brood year).} \item{\code{A}}{Spawning habitat size (either stream
#'   length or area). Will usually be time-invariant within a population, but
#'   need not be.} \item{\code{M_obs}}{Total number of wild-origin smolts (only
#'   needed for models including smolt stage).}
#'   \item{\code{n_Mage[min_Mage]_obs...n_Mage[max_Mage]_obs}}{If
#'   \code{life_cycle=="SMaS"}, multiple columns of observed smolt age
#'   frequencies (i.e., counts), where \code{[min_Mage]} and \code{[max_Mage]}
#'   are the numeral age in years of the youngest and oldest smolts,
#'   respectively. Note that age is measured in calendar years from the brood
#'   year (i.e., the Gilbert-Rich system).} \item{\code{S_obs}}{Total number
#'   (not density) of wild and hatchery-origin spawners.}
#'   \item{\code{n_age[min_age]_obs...n_age[max_age]_obs}}{Multiple columns of
#'   observed spawner age frequencies (i.e., counts), where \code{[min_age]} and
#'   \code{[max_age]} are the numeral age in years (total, not ocean age) of the
#'   youngest and oldest spawners, respectively.}
#'   \item{\code{n_MSage[min_MSage]_obs...n_MSage[max_MSage]_obs}}{If
#'   \code{life_cycle=="SMaS"}, multiple columns of observed ocean age
#'   frequencies (i.e., counts), where \code{[min_MSage]} and \code{[max_MSage]}
#'   are the youngest and oldest ocean age in years, respectively.}
#'   \item{\code{n_GRage[min_age]_[min_Mage]_obs...n_GRage[max_age]_[max_Mage]_obs}}{If
#'    \code{life_cycle=="SMaS"}, multiple columns of observed Gilbert-Rich age
#'   frequencies, sorted first by smolt age (\code{min_Mage:max_Mage}) and then
#'   by total age \code{min_age:max_age}. For example, a life history with
#'   subyearling or yearling smolts and ocean ages 2:3 would have column names
#'   c("n_GRage_3_1_obs", "n_GRage_4_1_obs", "n_GRage_4_2_obs",
#'   "n_GRage_5_2_obs")} \item{\code{n_W_obs}}{Observed frequency of
#'   natural-origin spawners.} \item{\code{n_H_obs}}{Observed frequency of
#'   hatchery-origin spawners.} \item{\code{fit_p_HOS}}{Logical or 0/1
#'   indicating for each row in fish_data whether the model should estimate
#'   \code{p_HOS > 0}. This is only required if \code{model == "IPM"}.}
#'   \item{\code{F_rate}}{Total harvest rate (proportion) of natural-origin
#'   fish.} \item{\code{B_take_obs}}{Number of adults taken for hatchery
#'   broodstock.} }
#' @param fish_data_fwd Only if model == "IPM", life_cycle == "SS", and
#'   pool_pops == TRUE, optional data frame with the following \code{colnames},
#'   representing "forward" or "future" simulations. Unlike \code{fish_data}, a
#'   given combination of population and year may occur multiple times, perhaps
#'   to facilitate comparisons across scenarios or "branches" with different
#'   inputs (e.g., harvest rate). In this case, all branches are subjected to
#'   the same sequence of process errors in recruitment and age structure.
#'   \describe{ \item{\code{pop}}{Numeric or character population ID. All values
#'   must also appear in \code{fish_data$pop}.} \item{\code{year}}{Integer
#'   variable giving the year the fish spawned (i.e., the brood year). For each
#'   population in \code{fish_data_fwd$pop}, the first year appearing in
#'   \code{fish_data_fwd$year} must be one greater than the last year appearing
#'   in \code{fish_data$year}, i.e.,
#'   \code{min(fish_data_fwd$year[fish_data_fwd$pop==j]) ==
#'   max(fish_data$year[fish_data$pop==j]) + 1}.} \item{\code{A}}{Spawning
#'   habitat size (either stream length or area). Will usually be time-invariant
#'   within a population, but need not be.} \item{\code{F_rate}}{Total harvest
#'   rate (proportion) of natural-origin fish.} \item{\code{B_rate}}{Total
#'   broodstock removal rate (proportion) of natural-origin fish.}
#'   \item{\code{p_HOS}}{Proportion of hatchery-origin spawners.} }
#' @param env_data Optional data frame or named list of data frames whose
#'   variables are time-varying environmental covariates, sequentially ordered
#'   with each row corresponding to a unique year in fish_data. If a named list,
#'   element names correspond to stage- or transition-specific covariate
#'   matrices defined in the Stan model being used. (This is required for
#'   multi-stage models.)
#' @param catch_data Only if model == "IPM_F", a data frame with numeric columns
#'   \describe{ \item{\code{year}}{Year for fishery data. Must be identical to
#'   \code{unique(fish_data$year)}.} \item{\code{R_F_obs}}{Total recruits to the
#'   fishery.} \item{\code{C_obs}}{Total catch.} }
#' @param ages For multi-stage models, a named list giving the fixed ages in
#'   years of all subadult life stages.
#' @param stan_model Character string giving the name of the Stan model being
#'   fit (".stan" filetype extension is not included).
#' @param SR_fun One of \code{"exp"}, \code{"BH"} (the default), or
#'   \code{"Ricker"}, indicating which spawner-recruit function to fit.
#'
#' @return A named list that can be passed to \code{stan} as the \code{data}
#'   argument.
#'
#' @export

stan_data <- function(fish_data, fish_data_fwd = NULL, env_data = NULL, catch_data = NULL, ages = NULL,
                      stan_model, SR_fun = "BH")
{
  fish_data <- as.data.frame(fish_data)
  life_cycle <- strsplit(stan_model, "_")[[1]][2]
  
  if(!is.null(fish_data_fwd))
  {
    if(stan_model != "IPM_SS_pp")
      warning("Argument fish_data_fwd is ignored unless stan_model == 'IPM_SS_pp'.\n")
    N_fwd <- nrow(fish_data_fwd)
    fish_data_fwd <- as.data.frame(fish_data_fwd)
    fish_data_fwd$pop <- factor(fish_data_fwd$pop, levels = levels(factor(fish_data$pop)))
    if(any(!fish_data_fwd$pop %in% fish_data$pop))
      stop("All populations in fish_data_fwd must appear in fish_data.\n")
    year_check <- tapply(fish_data$year, fish_data$pop, max)
    year_check <- year_check[names(year_check) %in% fish_data_fwd$pop]
    year_fwd_check <- tapply(fish_data_fwd$year, fish_data_fwd$pop, min)
    year_fwd_check <- year_fwd_check[names(year_fwd_check) %in% fish_data_fwd$pop]
    if(any(year_fwd_check != year_check + 1))
      stop("First year in fish_data_fwd must equal 1 + last year in fish_data for each population.\n")
    fish_data_fwd$pop <- as.numeric(fish_data_fwd$pop)
    fish_data_fwd$year <- as.numeric(factor(fish_data_fwd$year, 
                                            levels = levels(factor(c(fish_data$year, fish_data_fwd$year)))))
  }
  
  if(is.null(fish_data_fwd))
  {
    N_fwd <- 0
    fish_data_fwd <- data.frame(pop = 1, year = 1, A = 0, F_rate = 0, B_rate = 0, p_HOS = 0)
    fish_data_fwd <- fish_data_fwd[c(),]
  }
  
  fish_data$pop <- as.numeric(factor(fish_data$pop))
  fish_data$year <- as.numeric(factor(fish_data$year))
  fish_data$fit_p_HOS <- as.logical(fish_data$fit_p_HOS)
  
  for(i in c("pop","year","A","fit_p_HOS","B_take_obs"))
    if(any(is.na(fish_data[,i])))
      stop(paste0("Missing values not allowed in fish_data$", i, "\n"))
  
  if(any(is.na(fish_data_fwd)))
    stop("Missing values not allowed in fish_data_fwd.\n")
  
  if(is.null(env_data))
    env_data <- switch(life_cycle,
                       SS = list(matrix(0, max(fish_data$year, fish_data_fwd$year), 0)),
                       SMS = list(M = matrix(0, max(fish_data$year), 0),
                                  MS = matrix(0, max(fish_data$year), 0)),
                       SMaS = list(M = matrix(0, max(fish_data$year), 0),
                                   MS = matrix(0, max(fish_data$year), 0)))
  if(is.matrix(env_data) | is.data.frame(env_data)) env_data <- list(env_data)

  if(any(sapply(env_data, nrow) != max(fish_data$year, fish_data_fwd$year)))
    stop("Length of environmental time series does not equal number of brood years.\n")

  if(any(sapply(env_data, is.na)))
    stop("Missing values are not allowed in environmental covariates.\n")
  
  if(life_cycle != "SS" & any(is.na(ages) | is.null(ages)))
    stop("Multi-stage models must specify age in years for all stages.\n")
  
  if(stan_model == 'IPM_SS_F_pp') 
  {
    if(is.null(catch_data) | any(is.na(catch_data)))
      stop("Missing values are not allowed in total run size and catch with stan_model == 'IPM_SS_F_pp'.\n")
  }
  
  if(life_cycle == "SMaS")
  {
    max_Mage <- max(as.numeric(substring(names(fish_data)[grep("n_Mage", names(fish_data))], 7, 7)))
    max_MSage <- max(as.numeric(substring(names(fish_data)[grep("n_MSage", names(fish_data))], 8, 8))) 
    max_age <- max_Mage + max_MSage
  } else {
    max_age <- max(as.numeric(substring(names(fish_data)[grep("n_age", names(fish_data))], 6, 6)))
  }
  F_rate_check <- tapply(fish_data$F_rate, fish_data$pop, function(x) any(is.na(x[-c(1:max_age)])))
  if(any(F_rate_check))
    stop(paste0("Missing values not allowed in fish_data$F_rate except in years 1:max_age", i, "\n"))
  
  if(any(is.na(fish_data$n_W_obs) != is.na(fish_data$n_H_obs)))
    stop(paste("Conflicting NAs in n_W_obs and n_H_obs in rows", 
               which(is.na(fish_data$n_W_obs) != is.na(fish_data$n_H_obs))), "\n")
  
  age_NA_check <- is.na(fish_data[,grep("n_age", names(fish_data))])
  if(any(!rowSums(age_NA_check) %in% c(0, nrow(age_NA_check))))
    stop(paste("Conflicting NAs in age frequency data in rows", 
               which(!rowSums(age_NA_check) %in% c(0, nrow(age_NA_check))), "\n"))
  
  if(!stan_model %in% c("IPM_SS_np", "IPM_SS_pp","IPM_SS_F_pp","IPM_SMS_np","IPM_SMaS_np", 
                        "RR_SS_np","RR_SS_pp"))
    stop(paste("Stan model", stan_model, "does not exist.\n"))
  
  if(stan_model %in% c("IPM_SS_np","IPM_SS_pp"))
  {
    with(fish_data, {  
      dat <- list(SR_fun = switch(SR_fun, exp = 1, BH = 2, Ricker = 3),
                  N = nrow(fish_data),
                  pop = pop, 
                  year = year,
                  N_X = ncol(env_data[[1]]), 
                  X = as.matrix(env_data[[1]]),
                  N_S_obs = sum(!is.na(S_obs)),
                  which_S_obs = array(which(!is.na(S_obs)), dim = sum(!is.na(S_obs))),
                  S_obs = replace(S_obs, is.na(S_obs) | S_obs==0, 1),
                  N_age = sum(grepl("n_age", names(fish_data))), 
                  max_age = max_age,
                  n_age_obs = as.matrix(fish_data[,grep("n_age", names(fish_data))]),
                  N_H = sum(fit_p_HOS),
                  which_H = array(which(fit_p_HOS), dim = sum(fit_p_HOS)),
                  n_W_obs = array(n_W_obs[fit_p_HOS], dim = sum(fit_p_HOS)),
                  n_H_obs = array(n_H_obs[fit_p_HOS], dim = sum(fit_p_HOS)),
                  A = A,
                  F_rate = replace(F_rate, is.na(F_rate), 0),
                  N_B = sum(B_take_obs > 0),
                  which_B = array(which(B_take_obs > 0), dim = sum(B_take_obs > 0)),
                  B_take_obs = B_take_obs[B_take_obs > 0],
                  N_fwd = N_fwd,
                  pop_fwd = array(fish_data_fwd$pop, dim = nrow(fish_data_fwd)),
                  year_fwd = array(fish_data_fwd$year, dim = nrow(fish_data_fwd)),
                  A_fwd = array(fish_data_fwd$A, dim = nrow(fish_data_fwd)),
                  B_rate_fwd = array(fish_data_fwd$B_rate, dim = nrow(fish_data_fwd)),
                  F_rate_fwd = array(fish_data_fwd$F_rate, dim = nrow(fish_data_fwd)),
                  p_HOS_fwd = array(fish_data_fwd$p_HOS, dim = nrow(fish_data_fwd)))
      
      dat$n_W_obs[is.na(dat$n_W_obs)] <- 0
      dat$n_H_obs[is.na(dat$n_H_obs)] <- 0
      dat$n_age_obs[is.na(dat$n_age_obs)] <- 0
      
      return(dat)
    })
  } else if(stan_model %in% c("IPM_SMS_np")) {
    with(fish_data, {  
      dat <- list(SR_fun = switch(SR_fun, exp = 1, BH = 2, Ricker = 3),
                  N = nrow(fish_data),
                  pop = pop, 
                  year = year,
                  N_X_M = ncol(env_data$M), 
                  X_M = as.matrix(env_data$M),
                  N_X_MS = ncol(env_data$MS), 
                  X_MS = as.matrix(env_data$MS),
                  N_S_obs = sum(!is.na(S_obs)),
                  which_S_obs = array(which(!is.na(S_obs)), dim = sum(!is.na(S_obs))),
                  S_obs = replace(S_obs, is.na(S_obs) | S_obs==0, 1),
                  N_M_obs = sum(!is.na(M_obs)),
                  which_M_obs = array(which(!is.na(M_obs)), dim = sum(!is.na(M_obs))),
                  M_obs = replace(M_obs, is.na(M_obs) | M_obs==0, 1),
                  N_age = sum(grepl("n_age", names(fish_data))), 
                  smolt_age = ages$M,
                  max_age = max_age,
                  n_age_obs = as.matrix(fish_data[,grep("n_age", names(fish_data))]),
                  N_H = sum(fit_p_HOS),
                  which_H = array(which(fit_p_HOS), dim = sum(fit_p_HOS)),
                  n_W_obs = array(n_W_obs[fit_p_HOS], dim = sum(fit_p_HOS)),
                  n_H_obs = array(n_H_obs[fit_p_HOS], dim = sum(fit_p_HOS)),
                  A = A,
                  F_rate = replace(F_rate, is.na(F_rate), 0),
                  N_B = sum(B_take_obs > 0),
                  which_B = array(which(B_take_obs > 0), dim = sum(B_take_obs > 0)),
                  B_take_obs = B_take_obs[B_take_obs > 0])
      
      dat$n_W_obs[is.na(dat$n_W_obs)] <- 0
      dat$n_H_obs[is.na(dat$n_H_obs)] <- 0
      dat$n_age_obs[is.na(dat$n_age_obs)] <- 0
      
      return(dat)
    })
  } else if(stan_model %in% c("IPM_SMaS_np")) {
    with(fish_data, {  
      dat <- list(SR_fun = switch(SR_fun, exp = 1, BH = 2, Ricker = 3),
                  N = nrow(fish_data),
                  pop = pop, 
                  year = year,
                  N_X_M = ncol(env_data$M), 
                  X_M = as.matrix(env_data$M),
                  N_M_obs = sum(!is.na(M_obs)),
                  which_M_obs = array(which(!is.na(M_obs)), dim = sum(!is.na(M_obs))),
                  M_obs = replace(M_obs, is.na(M_obs) | M_obs==0, 1),
                  N_Mage = sum(grepl("n_Mage", names(fish_data))),
                  max_Mage = max_Mage,
                  n_Mage_obs = as.matrix(fish_data[,grep("n_Mage", names(fish_data))]),
                  N_X_MS = ncol(env_data$MS), 
                  X_MS = as.matrix(env_data$MS),
                  N_MSage = sum(grepl("n_MSage", names(fish_data))), 
                  max_MSage = max_MSage,
                  n_MSage_obs = as.matrix(fish_data[,grep("n_MSage", names(fish_data))]),
                  N_S_obs = sum(!is.na(S_obs)),
                  which_S_obs = array(which(!is.na(S_obs)), dim = sum(!is.na(S_obs))),
                  S_obs = replace(S_obs, is.na(S_obs) | S_obs==0, 1),
                  n_GRage_obs = as.matrix(fish_data[,grep("n_GRage", names(fish_data))]),
                  N_H = sum(fit_p_HOS),
                  which_H = array(which(fit_p_HOS), dim = sum(fit_p_HOS)),
                  n_W_obs = array(n_W_obs[fit_p_HOS], dim = sum(fit_p_HOS)),
                  n_H_obs = array(n_H_obs[fit_p_HOS], dim = sum(fit_p_HOS)),
                  A = A,
                  F_rate = replace(F_rate, is.na(F_rate), 0),
                  N_B = sum(B_take_obs > 0),
                  which_B = array(which(B_take_obs > 0), dim = sum(B_take_obs > 0)),
                  B_take_obs = B_take_obs[B_take_obs > 0])
      
      dat$n_W_obs[is.na(dat$n_W_obs)] <- 0
      dat$n_H_obs[is.na(dat$n_H_obs)] <- 0
      dat$n_age_obs[is.na(dat$n_age_obs)] <- 0
      
      return(dat)
    })
  } else if(stan_model %in% c("RR_SS_np","RR_SS_pp")) {
    recon_dat <- run_recon(fish_data)
    which_fit <- which(!is.na(recon_dat$R) & !is.na(recon_dat$S))
    N_fit <- length(which_fit)
    p <- aggregate(recon_dat[,grep("p_age", names(recon_dat))], list(pop = recon_dat$pop), mean, na.rm = TRUE)[,-1]
    
    with(recon_dat, {
      dat <- list(SR_fun = switch(SR_fun, exp = 1, BH = 2, Ricker = 3),
                  N = length(S),
                  pop = pop, 
                  year = year,
                  N_fit = N_fit,
                  which_fit = array(which_fit, dim = N_fit),
                  S = replace(S, S == 0 | is.na(S), 1),
                  R = replace(R, R == 0 | is.na(R), 1),
                  A = A,
                  S_NA = array(as.integer(is.na(S)), dim = length(S)),
                  R_NA = array(as.integer(is.na(R)), dim = length(R)),
                  N_age = sum(grepl("p_age", names(recon_dat))),
                  max_age = max(as.numeric(substring(names(recon_dat)[grep("p_age", names(recon_dat))], 6, 6))),
                  p = as.matrix(p))
      
      return(dat)
    })
  }
}
