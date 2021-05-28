#' Fits an integrated or run-reconstruction spawner-recruit model.
#'
#' @param model Either `"IPM"` or `"RR"`, indicating whether the data
#'   are intended for an integrated or run-reconstruction model.
#' @param life_cycle Character string indicating which life-cycle model to fit.
#' One of the following options (must be `"SS"` if `model == "RR"`):
#'   * `"SS"`  Spawner-to-spawner (the default)
#'   * `"SMS"`  Spawner-smolt-spawner
#'   * `"SMaS"`  Spawner-smolt-spawner with multiple smolt age classes (currently only
#'   available for `pool_pops == FALSE`)
#'   * `"LCRchum"`  Customized spawner-egg-smolt-spawner model for Lower Columbia chum
#'   (`pool_pops == TRUE`)
#'   * `"ICchinook"`  Customized spawner-smolt-spawner model with downstream, SAR,
#'   and upstream survival (`pool_pops == TRUE`)
#' @param pool_pops Logical defaulting to `TRUE`, indicating whether or not
#'   to treat the different populations as hierarchical rather than
#'   fixed/independent.
#' @param stan_model Character string specifying the **salmonIPM** model to be
#'   fit. A more concise alternative to specifying `model`, `life_cycle`, and `pool_pops` 
#'   (and will override those arguments).
#' @param SR_fun One of `"exp"` (density-independent discrete exponential), 
#'   `"BH"` (Beverton-Holt, the default), or `"Ricker"`, indicating which spawner-recruit
#'   function to fit. Synonyms `"B-H"`, `"bh"`, `"b-h"` and `"ricker"` are also accepted.
#' @param par_models  Optional list of two-sided formulas of the form 
#' `theta ~ t1 + ... + tK`, where `theta` is a parameter or state in the model
#' specified by `stan_model` that accepts covariates (see Details for available
#' model-response combinations) and `t1 ... tK` are terms involving variables in
#' `fish_data`. Standard formula syntax such as `:` and `*` may be used;
#' see [stats::formula()].
#' @param scale  Logical indicating whether the model matrices constructed from
#' `fish_data` using the formulas in `par_models` should be scaled to have 
#' column SDs of 1 in addition to being centered (`TRUE`) or centered only (`FALSE`). 
#' @param ages For multi-stage models, a named list giving the fixed ages in
#'   years of all subadult life stages. (This is not needed for `IPM_SMaS_np` because
#'   in that case smolt age structure is provided in `fish_data`.)
#' @param age_S_obs If `stan_model == "IPM_SS_pp"`, an optional logical or binary numeric
#'   vector indicating, for each adult age, whether observed total spawner data
#'   includes that age. The default is to treat `S_obs` as including spawners of
#'   all ages.
#' @param age_S_eff If `stan_model == "IPM_SS_pp"`, an optional logical or binary numeric
#'   vector indicating, for each adult age, whether spawners of that age
#'   contribute toward reproduction. This could be used, e.g., to exclude jacks
#'   from the effective breeding population. The default is to include spawners
#'   of all ages.
#' @param conditionGRonMS If `stan_model == "IPM_SMaS_np`, logical indicating
#' whether the Gilbert-Rich age frequencies `n_GRage_obs` in `fish_data` are conditioned
#' on ocean age. If `FALSE` (the default) the counts are assumed to be sampled randomly
#' from the population. If `TRUE`, it is assumed that the number of spawners of each ocean
#' age is arbitrary, but smolt (FW) age is randomly sampled within each ocean age; i.e.,
#' in a smolt age [x] ocean age contingency table, the cell frequencies are conditioned
#' on the column totals. 
#' @param fish_data Data frame that includes the following `colnames`, in
#'   no particular order except where noted: 
#'   * `pop`  Numeric or character population ID.  
#'   * `year`  Numeric variable giving the year the fish spawned (i.e., the brood year).
#'   * `A`  Spawning habitat size (either stream length or area). 
#'   Will often be time-invariant within a population, but need not be.   
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
#'   * `n_Mage[min_Mage]_obs...n_Mage[max_Mage]_obs`  If `life_cycle == "SMaS"`, 
#'   multiple columns of observed smolt age frequencies (i.e., counts), 
#'   where `[min_Mage]` and `[max_Mage]` are the numeral age in years of the youngest 
#'   and oldest smolts, respectively. Note that age is measured in calendar years from 
#'   the brood year (i.e., the Gilbert-Rich system).   
#'   * `S_obs`  Total number (not density) of all wild and hatchery-origin spawners.  
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
#'   `c("n_GRage_3_1_obs", "n_GRage_4_1_obs", "n_GRage_4_2_obs", "n_GRage_5_2_obs")` 
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
#'   * `S_add_obs`  If `stan_model == "IPM_LCRchum_pp"`, number of adults translocated into 
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
#' @param fecundity_data If `stan_model == "IPM_LCRchum_pp"`, data frame with
#' the following columns, representing observations of fecundity with each
#' row corresponding to a female:
#' * `age_E`  Female age in years.   
#' * `E_obs`  Observed fecundity.  
#' @param prior_data If `stan_model == "IPM_ICchinook_pp"`, named list
#' with the following elements: 
#' * `s`  Data frame with columns `year`, `mu_prior_D`, `sigma_prior_D`, 
#' `mu_prior_SAR`, `sigma_prior_SAR`, `mu_prior_U`, `sigma_prior_U`, 
#' giving the annual prior means and SDs of logit survival downstream, at sea, 
#' and upstream, respectively.
#' @param init A list of named lists of initial values to be passed to
#'   [rstan::stan()]. If `NULL`, initial values will be automatically
#'   generated from the supplied data using [stan_init()].
#' @param pars Vector of character strings specifying parameters to monitor.
#'   If NULL, default values are used. If a non-default value is supplied, the
#'   user should make sure the parameters requested appear in the model configuration 
#'   specified.
#' @param include Logical scalar defaulting to `TRUE` indicating whether to include
#'  or exclude the parameters given by `pars`. If `FALSE`, only entire 
#'  multidimensional parameters can be excluded, rather than particular elements of them. 
#' @param log_lik Logical scalar indicating whether the pointwise log-likelihood
#'   should be returned for later analysis with [loo::loo()].
#' @param chains Positive integer specifying the number of MCMC chains; 
#'   see [rstan::stan()].
#' @param iter Positive integer specifying the number of iterations for each
#'   chain (including warmup); see [rstan::stan()].
#' @param warmup Positive integer specifying the number of warmup (aka burnin)
#'   iterations per chain. If step-size adaptation is on (which it is by
#'   default), this also controls the number of iterations for which adaptation
#'   is run (and hence these warmup samples should not be used for inference).
#'   The number of warmup iterations should not be larger than `iter`.
#'   See [rstan::stan()].
#' @param thin Positive integer specifying the period for saving samples. The
#'   default is 1, which is usually the recommended value. See [rstan::stan()].
#' @param cores Number of cores to use when executing the chains in parallel.
#'   Defaults to one less than the number of cores available. See [rstan::stan()].
#' @param ... Additional arguments to pass to [rstan::sampling()].
#'   
#' @details 
#' The following parameters and states can be modeled as functions of covariates using
#' the argument `par_models`. The response distribution families are automatically 
#' implemented so there is no need to `log()`- or `qlogis()`-transform the left-hand side
#' of the formula, although such syntax will also work. (This is because the LHS is not
#' found in `fish_data` and is only used to determine the parameter name.) The design
#' matrices passed to the Stan model cannot include an intercept, but it is not necessary 
#' to manually remove it in the RHS; if present by default, [par_model_matrix()] will 
#' automatically remove it.
#' 
#' As with any regression model, the user must ensure the effects specified are 
#' estimable given the design matrix. For example, a spatially varying but 
#' time-invariant predictor would not be identifiable in a `_np` model because
#' populations are modeled independently.
#' 
#' |                    |                             |                                 |                            | **Response (family)**      |                        |                         |                                 |
#' |:-------------------|:---------------------------:|:-------------------------------:|:--------------------------:|:--------------------------:|:----------------------:|:-----------------------:|:-------------------------------:|    
#' | **Model**          | **`alpha` \cr (lognormal)** | **`psi` \cr (logistic normal)** | **`Rmax` \cr (lognormal)** | **`Mmax` \cr (lognormal)** | **`R` \cr (lognormal)**| **`M` \cr (lognormal)** | **`s_MS` \cr (logistic normal)**|
#' | `IPM_SS_np`        | &#x2611;                    | &#x2610;                        | &#x2611;                   | &#x2610;                   | &#x2611;               | &#x2610;                | &#x2610;                        | 
#' | `IPM_SS_pp`        | &#x2611;                    | &#x2610;                        | &#x2611;                   | &#x2610;                   | &#x2611;               | &#x2610;                | &#x2610;                        | 
#' | `IPM_SMS_np`       | &#x2611;                    | &#x2610;                        | &#x2610;                   | &#x2611;                   | &#x2610;               | &#x2611;                | &#x2611;                        | 
#' | `IPM_SMS_pp`       | &#x2611;                    | &#x2610;                        | &#x2610;                   | &#x2611;                   | &#x2610;               | &#x2611;                | &#x2611;                        | 
#' | `IPM_SMaS_np`      | &#x2611;                    | &#x2610;                        | &#x2610;                   | &#x2611;                   | &#x2610;               | &#x2611;                | &#x2611;                        | 
#' | `IPM_ICchinook_pp` | &#x2611;                    | &#x2610;                        | &#x2610;                   | &#x2611;                   | &#x2610;               | &#x2611;                | &#x2610;                        | 
#' | `IPM_LCRchum_pp`   | &#x2610;                    | &#x2611;                        | &#x2610;                   | &#x2611;                   | &#x2610;               | &#x2611;                | &#x2611;                        | 
#' 
#' @return An object of class `stanfit` representing the fitted model. See
#'   [rstan::stan()] for details.
#'
#' @importFrom rstan sampling
#' @export
#' @encoding UTF-8

salmonIPM <- function(model = "IPM", life_cycle = "SS", pool_pops = TRUE, stan_model = NULL, 
                      SR_fun = "BH", par_models = NULL, scale = TRUE, 
                      ages = NULL, age_S_obs = NULL, age_S_eff = NULL, conditionGRonMS = FALSE,
                      fish_data, fish_data_fwd = NULL, fecundity_data = NULL, prior_data = NULL,
                      init = NULL, pars = NULL, include = TRUE, log_lik = FALSE, 
                      chains, iter, warmup, thin = 1, cores = parallel::detectCores() - 1, 
                      ...)
{
  if(is.null(stan_model)) 
    stan_model <- paste(model, life_cycle, ifelse(pool_pops, "pp", "np"), sep = "_")
  if(SR_fun %in% c("B-H","bh","b-h")) SR_fun <- "BH"
  if(SR_fun == "ricker") SR_fun <- "Ricker"
  
  dat <- stan_data(stan_model = stan_model, SR_fun = SR_fun, 
                   par_models = par_models, scale = scale, 
                   ages = ages, age_S_obs = age_S_obs, age_S_eff = age_S_eff, 
                   fish_data = fish_data, fish_data_fwd = fish_data_fwd, 
                   fecundity_data = fecundity_data, prior_data = prior_data, 
                   conditionGRonMS = conditionGRonMS)
  
  if(is.null(pars)) 
  {
    pars <- stan_pars(stan_model)
  } else if(!include) {
    pars <- setdiff(stan_pars(stan_model), pars)
  }
  if(log_lik) pars <- c(pars, "LL")
  
  fit <- rstan::sampling(stanmodels[[stan_model]],
                         data = dat, 
                         init = stan_init(stan_model = stan_model, dat = dat, chains = chains), 
                         pars = pars,
                         chains = chains, cores = cores, 
                         iter = iter, warmup = warmup, thin = thin, 
                         ...)
}