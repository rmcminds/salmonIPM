#' Assemble input data for integrated or run-reconstruction spawner-recruit
#' models.
#'
#' @param stan_model Character string giving the name of the **salmonIPM** model being
#'   fit (".stan" filetype extension is not included).
#' @param SR_fun One of `"exp"`, `"BH"` (the default), or `"Ricker"`, 
#' indicating which spawner-recruit function to fit.
#' @param ages For multi-stage models, a named list giving the fixed ages in
#'   years of all subadult life stages.
#' @param par_models  Optional list of two-sided formulas of the form 
#' `theta ~ t1 + ... + tK`, where `theta` is a parameter or state in the model
#' specified by `stan_model` that accepts covariates (see Details for available
#' model-response combinations) and `t1 ... tK` are terms involving variables in
#' `fish_data`. Standard formula syntax such as `:` and `*` may be used;
#' see [stats::formula()].
#' @param scale  Logical indicating whether the model matrices constructed from
#' `fish_data` using the formulas in `par_models` should be scaled to have 
#' column SDs of 1 in addition to being centered (`TRUE`) or centered only (`FALSE`). 
#' @param fish_data Data frame that includes the following `colnames`, in
#'   no particular order except where noted: 
#'   * `pop`  Numeric or character population ID.  
#'   * `year`  Numeric variable giving the year the fish spawned (i.e., the brood year).
#'   * `A`  Spawning habitat size (either stream length or area). 
#'   Will usually be time-invariant within a population, but need not be.   
#'   * `M_obs`  Total number of wild-origin smolts (only needed for models including smolt stage).  
#'   * `tau_M_obs`  If `stan_model == "IPM_LCRchum_pp"`,  known observation error SDs 
#'   for smolt abundance.  
#'   * `downstream_trap`  If `stan_model == "IPM_LCRchum_pp"`, row indices
#'   corresponding to a downstream smolt trap in a different population whose
#'   catch additionally includes the smolts produced in one or more upstream populations,
#'   assuming zero mortality. Each upstream population can have at most one 
#'   downstream trap (in addition to its own, if any) but a trap can have multiple 
#'   upstream populations. If `downstream_trap[i] == j`, `M_downstream[j] <- M[j] + M[i]`; 
#'   if `is.na(downstream_trap[i])` then `M[i]` is not double-counted.  
#'   * `n_Mage[min_Mage]_obs...n_Mage[max_Mage]_obs`  If `life_cycle=="SMaS"`, 
#'   multiple columns of observed smolt age frequencies (i.e., counts), 
#'   where `[min_Mage]` and `[max_Mage]` are the numeral age in years of the youngest 
#'   and oldest smolts, respectively. Note that age is measured in calendar years from 
#'   the brood year (i.e., the Gilbert-Rich system).   
#'   * `S_obs`  Total number (not density) of wild and hatchery-origin spawners.  
#'   * `tau_S_obs`  If `stan_model == "IPM_LCRchum_pp"`, known observation error SDs 
#'   for spawner abundance.  
#'   * `n_age[min_age]_obs...n_age[max_age]_obs`  Multiple columns of
#'   observed spawner age frequencies (i.e., counts), where `[min_age]` and
#'   `[max_age]` are the numeral age in years (total, not ocean age) of the
#'   youngest and oldest spawners, respectively.  
#'   * `n_MSage[min_MSage]_obs...n_MSage[max_MSage]_obs`  If `life_cycle == "SMaS"`, 
#'   multiple columns of observed ocean age frequencies (i.e., counts), 
#'   where `[min_MSage]` and `[max_MSage]` are the youngest and oldest ocean age 
#'   in years, respectively.  
#'   * `n_GRage[min_age]_[min_Mage]_obs...n_GRage[max_age]_[max_Mage]_obs`  If
#'    `life_cycle == "SMaS"`, multiple columns of observed Gilbert-Rich age
#'   frequencies, sorted first by smolt age (`min_Mage:max_Mage`) and then
#'   by total age `min_age:max_age`. For example, a life history with
#'   subyearling or yearling smolts and ocean ages 2:3 would have column names
#'   c("n_GRage_3_1_obs", "n_GRage_4_1_obs", "n_GRage_4_2_obs", "n_GRage_5_2_obs") 
#'   * `n_M_obs`  If `stan_model == "IPM_LCRchum_pp"`, observed frequency of male spawners.
#'   * `n_F_obs`  If `stan_model == "IPM_LCRchum_pp"`, observed frequency of female spawners.
#'   * `p_G_obs`  If `stan_model == "IPM_LCRchum_pp"`, observed proportion (assumed known
#'   without error) of female spawners that are "green", i.e. fully fecund.
#'   * `n_W_obs`  Observed frequency of natural-origin spawners.   
#'   * `n_H_obs`  Observed frequency of hatchery-origin spawners.   
#'   * `fit_p_HOS`  Logical or 0/1 indicating for each row in fish_data whether the 
#'   model should estimate `p_HOS > 0`. This is only required if `model == "IPM"`.  
#'   * `F_rate`  Total harvest rate (proportion) of natural-origin fish.   
#'   * `B_take_obs`  Number of adults taken for hatchery broodstock.   
#'   * `S_add_obs`  If `stan_model == "IPM_LCRchum_pp`, number of adults translocated into 
#'   population. 
#'   * `...`  Additional variables to be used as covariates.  
#' @param fish_data_fwd Only if `stan_model == "IPM_SS_pp"`, optional data frame
#'   with the following columns, representing "forward" or "future"
#'   simulations: 
#'   * `pop`  Numeric or character population ID. All values must also appear in `fish_data$pop`.   
#'   * `year`  Integer variable giving the year the fish spawned (i.e., the brood year). 
#'   For each population in `fish_data_fwd$pop`, the first year appearing in 
#'   `fish_data_fwd$year` must be one greater than the last year appearing in 
#'   `fish_data$year`, i.e., 
#'   `min(fish_data_fwd$year[fish_data_fwd$pop == j]) == max(fish_data$year[fish_data$pop == j]) + 1`.
#'   * `A`  Spawning habitat size (either stream length or area). 
#'   Will usually be time-invariant within a population, but need not be.   
#'   * `F_rate`  Total harvest rate (proportion) of natural-origin fish.   
#'   * `B_rate`  Total broodstock removal rate (proportion) of natural-origin fish.  
#'   * `p_HOS`  Proportion of hatchery-origin spawners.
#'   
#'   Unlike `fish_data`, a given combination of population and
#'   year may occur multiple times, perhaps to facilitate comparisons across
#'   scenarios or "branches" with different inputs (e.g., harvest rate). In this
#'   case, all branches are subjected to the same sequence of process errors in
#'   recruitment and age structure. 
#' @param fecundity_data Only if `stan_model == "IPM_LCRchum_pp"`, data frame with
#' the following columns, representing observations of fecundity with each
#' row corresponding to a female:
#' * `age_E`  Female age in years.   
#' * `E_obs`  Observed fecundity.  
#' @param prior_data Only if `stan_model == "IPM_ICchinook_pp"`, named list
#' with the following elements: 
#' * `s`  Data frame with columns `year`, `mu_prior_D`, `sigma_prior_D`, 
#' `mu_prior_SAR`, `sigma_prior_SAR`, `mu_prior_U`, `sigma_prior_U`, 
#' giving the annual prior means and SDs of logit survival downstream, at sea, 
#' and upstream, respectively.
#' @param age_S_obs Only if `stan_model == "IPM_SSpa_pp"`, a logical or numeric
#'   vector indicating, for each adult age, whether observed total spawner data
#'   includes that age. The default is to treat `S_obs` as including spawners of
#'   all ages.
#' @param age_S_eff Only if `stan_model == "IPM_SSpa_pp"`, a logical or numeric
#'   vector indicating, for each adult age, whether spawners of that age
#'   contribute toward reproduction. This could be used, e.g., to exclude jacks
#'   from the effective breeding population. The default is to include spawners
#'   of all ages.
#' @param conditionGRonMS Only if `stan_model == "IPM_SMaS_np`, logical indicating
#' whether the Gilbert-Rich age frequencies `n_GRage_obs` in `fish_data` are conditioned
#' on ocean age. If `FALSE` (the default) the counts are assumed to be sampled randomly
#' from the population. If `TRUE`, it is assumed that the number of spawners of each ocean
#' age is arbitrary, but smolt (FW) age is randomly sampled within each ocean age; i.e.,
#' in a `smolt age x ocean age` contingency table, the cell frequencies are conditioned
#' on the column totals. 
#'
#' @return A named list that can be passed to `stan` as the `data`
#'   argument.
#'
#' @export

stan_data <- function(stan_model, SR_fun = "BH", ages = NULL, 
                      par_models = NULL, scale = TRUE, fish_data, fish_data_fwd = NULL, 
                      fecundity_data = NULL, prior_data = NULL, 
                      age_S_obs = NULL, age_S_eff = NULL, conditionGRonMS = FALSE)
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
  
  if(stan_model == "IPM_ICchinook_pp") {
    if(is.null(prior_data)) {
      stop(paste("Priors for survival must be specified for model IPM_ICchinook_pp \n"))
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
  
  N <- nrow(fish_data)
  fish_data$pop <- as.numeric(factor(fish_data$pop))
  fish_data$year <- as.numeric(factor(fish_data$year))
  fish_data$fit_p_HOS <- as.logical(fish_data$fit_p_HOS)
  
  if(any(unlist(tapply(fish_data$year, fish_data$pop, diff)) != 1))
    stop(paste0("Non-consecutive years not allowed in fish_data \n"))

  for(i in c("pop","year","A","fit_p_HOS","B_take_obs"))
    if(any(is.na(fish_data[,i])))
      stop(paste0("Missing values not allowed in fish_data$", i, "\n"))
  
  if(any(is.na(fish_data_fwd)))
    stop("Missing values not allowed in fish_data_fwd.\n")
  
  X <- lapply(par_models, function(f) {
    stopifnot(attr(terms(f), "response") == 1) # formulas must be 2-sided
    m <- model.matrix(update(f, NULL ~ . + 1), data = fish_data) # force placeholder intercept 
    if(any(is.na(m))) stop("Missing values are not allowed in environmental covariates.\n")
    m <- subset(m, select = -`(Intercept)`, drop = FALSE)
    m <- scale(m, scale = scale)
    return(m)
  })
  names(X) <- lapply(par_models, function(f) all.vars(f)[1]) # names are responses

  if(stan_model == "IPM_LCRchum_pp") {
    if(is.null(fecundity_data)) {
      stop("Fecundity data must be provided for Lower Columbia chum IPM. \n") 
    } else if(any(is.na(fecundity_data))) {
      stop("Missing values are not allowed in fecundity data. \n")
    }
  }
  
  if(stan_model == "IPM_SSpa_pp" & is.null(age_S_obs))
    age_S_obs <- rep(1, sum(grepl("n_age", names(fish_data))))
  age_S_obs <- as.numeric(age_S_obs)
  
  if(stan_model == "IPM_SSpa_pp" & is.null(age_S_eff))
    age_S_eff <- rep(1, sum(grepl("n_age", names(fish_data))))
  age_S_eff <- as.numeric(age_S_eff)
  
  if(life_cycle != "SS" & any(is.na(ages) | is.null(ages)))
    stop("Multi-stage models must specify age in years for all stages. \n")
  
  if(life_cycle == "SMaS")
  {
    if(!is.logical(conditionGRonMS))
      stop("conditionGRonMS must be TRUE or FALSE for model 'IPM_SMaS_np'. \n")
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
  
  if(stan_model == "IPM_LCRchum_pp" & any(is.na(fish_data$n_M_obs) != is.na(fish_data$n_F_obs)))
    stop(paste("Conflicting NAs in n_M_obs and n_F_obs in rows", 
               which(is.na(fish_data$n_M_obs) != is.na(fish_data$n_F_obs))), "\n")
  
  age_NA_check <- is.na(fish_data[,grep("n_age", names(fish_data))])
  if(any(!rowSums(age_NA_check) %in% c(0, nrow(age_NA_check))))
    stop(paste("Conflicting NAs in age frequency data in rows", 
               which(!rowSums(age_NA_check) %in% c(0, nrow(age_NA_check))), "\n"))
  
  if(!stan_model %in% c("IPM_SS_np","IPM_SS_pp","IPM_SSpa_pp","IPM_SMS_np","IPM_SMS_pp",
                        "IPM_SMaS_np","IPM_LCRchum_pp","IPM_ICchinook_pp",
                        "RR_SS_np","RR_SS_pp"))
    stop(paste("Stan model", stan_model, "does not exist.\n"))
  
  if(stan_model %in% c("IPM_SS_np","IPM_SS_pp","IPM_SSpa_pp"))
  {
    with(fish_data, {  
      dat <- list(
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
        N_X_alpha = ifelse(is.null(X$alpha), 0, ncol(X$alpha)), 
        X_alpha = if(is.null(X$alpha)) matrix(0,N,0) else X$alpha,
        N_X_Rmax = ifelse(is.null(X$Rmax), 0, ncol(X$Rmax)), 
        X_Rmax = if(is.null(X$Rmax)) matrix(0,N,0) else X$Rmax,
        N_X_R = ifelse(is.null(X$R), 0, ncol(X$R)), 
        X_R = if(is.null(X$R)) matrix(0,N,0) else X$R,
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
        n_age_obs = as.matrix(fish_data[,grep("n_age", names(fish_data))]),
        age_S_obs = as.vector(age_S_obs),
        age_S_eff = as.vector(age_S_eff),
        # H/W composition
        N_H = sum(fit_p_HOS),
        which_H = as.vector(which(fit_p_HOS)),
        n_W_obs = as.vector(n_W_obs[fit_p_HOS]),
        n_H_obs = as.vector(n_H_obs[fit_p_HOS])
      )
      
      dat$n_W_obs[is.na(dat$n_W_obs)] <- 0
      dat$n_H_obs[is.na(dat$n_H_obs)] <- 0
      dat$n_age_obs[is.na(dat$n_age_obs)] <- 0
      
      return(dat)
    })
  } else if(stan_model %in% c("IPM_SMS_np","IPM_SMS_pp")) {
    with(fish_data, {  
      dat <- list(
        # info for observed data
        N = N,
        pop = pop, 
        year = year,
        # smolt production
        SR_fun = switch(SR_fun, exp = 1, BH = 2, Ricker = 3),
        A = A,
        N_X_alpha = ifelse(is.null(X$alpha), 0, ncol(X$alpha)), 
        X_alpha = if(is.null(X$alpha)) matrix(0,N,0) else X$alpha,
        N_X_Mmax = ifelse(is.null(X$Mmax), 0, ncol(X$Mmax)), 
        X_Mmax = if(is.null(X$Mmax)) matrix(0,N,0) else X$Mmax,
        N_X_M = ifelse(is.null(X$M), 0, ncol(X$M)), 
        X_M = if(is.null(X$M)) matrix(0,N,0) else X$M,
        # smolt abundance
        N_M_obs = sum(!is.na(M_obs)),
        which_M_obs = as.vector(which(!is.na(M_obs))),
        M_obs = replace(M_obs, is.na(M_obs) | M_obs==0, 1),
        N_age = sum(grepl("n_age", names(fish_data))), 
        smolt_age = ages$M,
        # SAR
        N_X_MS = ifelse(is.null(X$s_MS), 0, ncol(X$s_MS)), 
        X_MS = if(is.null(X$s_MS)) matrix(0,N,0) else X$s_MS,
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
        n_age_obs = as.matrix(fish_data[,grep("n_age", names(fish_data))]),
        # H/W composition
        N_H = sum(fit_p_HOS),
        which_H = as.vector(which(fit_p_HOS)),
        n_W_obs = as.vector(n_W_obs[fit_p_HOS]),
        n_H_obs = as.vector(n_H_obs[fit_p_HOS])
      )
      
      dat$n_W_obs[is.na(dat$n_W_obs)] <- 0
      dat$n_H_obs[is.na(dat$n_H_obs)] <- 0
      dat$n_age_obs[is.na(dat$n_age_obs)] <- 0
      
      return(dat)
    })
  } else if(stan_model %in% c("IPM_SMaS_np")) {
    with(fish_data, {  
      dat <- list(
        # info for observed data
        N = N,
        pop = pop, 
        year = year,
        # smolt production
        SR_fun = switch(SR_fun, exp = 1, BH = 2, Ricker = 3),
        A = A,
        N_X_M = ifelse(is.null(X$M), 0, ncol(X$M)), 
        X_M = if(is.null(X$M)) matrix(0,N,0) else X$M,
        # smolt abundance
        N_M_obs = sum(!is.na(M_obs)),
        which_M_obs = as.vector(which(!is.na(M_obs))),
        M_obs = replace(M_obs, is.na(M_obs) | M_obs==0, 1),
        # smolt age structure
        N_Mage = sum(grepl("n_Mage", names(fish_data))),
        max_Mage = max_Mage,
        n_Mage_obs = as.matrix(fish_data[,grep("n_Mage", names(fish_data))]),
        # SAR
        N_X_MS = ifelse(is.null(X$s_MS), 0, ncol(X$s_MS)), 
        X_MS = if(is.null(X$s_MS)) matrix(0,N,0) else X$s_MS,
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
        n_W_obs = as.vector(n_W_obs[fit_p_HOS]),
        n_H_obs = as.vector(n_H_obs[fit_p_HOS])
      )
      
      dat$n_W_obs[is.na(dat$n_W_obs)] <- 0
      dat$n_H_obs[is.na(dat$n_H_obs)] <- 0
      dat$n_age_obs[is.na(dat$n_age_obs)] <- 0
      
      return(dat)
    })
  } else if(stan_model == "IPM_LCRchum_pp") {
    with(fish_data, {  
      dat <- list(
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
        N_X_psi = ifelse(is.null(X$psi), 0, ncol(X$psi)), 
        X_psi = if(is.null(X$psi)) matrix(0,N,0) else X$psi,
        N_X_Mmax = ifelse(is.null(X$Mmax), 0, ncol(X$Mmax)), 
        X_Mmax = if(is.null(X$Mmax)) matrix(0,N,0) else X$Mmax,
        N_X_M = ifelse(is.null(X$M), 0, ncol(X$M)), 
        X_M = if(is.null(X$M)) matrix(0,N,0) else X$M,
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
        N_X_MS = ifelse(is.null(X$s_MS), 0, ncol(X$s_MS)), 
        X_MS = if(is.null(X$s_MS)) matrix(0,N,0) else X$s_MS,
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
        n_age_obs = as.matrix(fish_data[,grep("n_age", names(fish_data))]),
        n_M_obs = as.vector(n_M_obs),
        n_F_obs = as.vector(n_F_obs),
        p_G_obs = as.vector(p_G_obs),
        # H/W composition
        N_H = sum(fit_p_HOS),
        which_H = as.vector(which(fit_p_HOS)),
        n_W_obs = as.vector(n_W_obs[fit_p_HOS]),
        n_H_obs = as.vector(n_H_obs[fit_p_HOS])
      )
      
      dat$n_W_obs[is.na(dat$n_W_obs)] <- 0
      dat$n_H_obs[is.na(dat$n_H_obs)] <- 0
      dat$n_age_obs[is.na(dat$n_age_obs)] <- 0
      
      return(dat)
    })
  } else if(stan_model == "IPM_ICchinook_pp") {
    with(fish_data, {  
      dat <- list(
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
        N_X_alpha = ifelse(is.null(X$alpha), 0, ncol(X$alpha)), 
        X_alpha = if(is.null(X$alpha)) matrix(0,N,0) else X$alpha,
        N_X_Mmax = ifelse(is.null(X$Mmax), 0, ncol(X$Mmax)), 
        X_Mmax = if(is.null(X$Mmax)) matrix(0,N,0) else X$Mmax,
        N_X_M = ifelse(is.null(X$M), 0, ncol(X$M)), 
        X_M = if(is.null(X$M)) matrix(0,N,0) else X$M,
        # downstream, SAR, upstream survival
        # N_X_D = ncol(env_data$s_D),
        # X_D = as.matrix(env_data$s_D),
        N_prior_D = N_prior_D,
        which_prior_D = as.vector(which_prior_D),
        mu_prior_D = as.vector(mu_prior_D),
        sigma_prior_D = as.vector(sigma_prior_D),
        # N_X_SAR = ncol(env_data$s_SAR),
        # X_SAR = as.matrix(env_data$s_SAR),
        N_prior_SAR = N_prior_SAR,
        which_prior_SAR = as.vector(which_prior_SAR),
        mu_prior_SAR = as.vector(mu_prior_SAR),
        sigma_prior_SAR = as.vector(sigma_prior_SAR),
        # N_X_U = ncol(env_data$s_U),
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
        n_age_obs = as.matrix(fish_data[,grep("n_age", names(fish_data))]),
        # H/W composition
        N_H = sum(fit_p_HOS),
        which_H = as.vector(which(fit_p_HOS)),
        n_W_obs = as.vector(n_W_obs[fit_p_HOS]),
        n_H_obs = as.vector(n_H_obs[fit_p_HOS])
      )
      
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
                  which_fit = as.vector(which_fit),
                  S = replace(S, S == 0 | is.na(S), 1),
                  R = replace(R, R == 0 | is.na(R), 1),
                  A = A,
                  S_NA = as.vector(as.integer(is.na(S))),
                  R_NA = as.vector(as.integer(is.na(R))),
                  N_age = sum(grepl("p_age", names(recon_dat))),
                  max_age = max(as.numeric(substring(names(recon_dat)[grep("p_age", names(recon_dat))], 6, 6))),
                  p = as.matrix(p))
      
      return(dat)
    })
  }
}
