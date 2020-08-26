#' Fits an integrated or run-reconstruction spawner-recruit model.
#'
#' @param fish_data See [stan_data()].
#' @param fish_data_fwd See [stan_data()].
#' @param env_data See [stan_data()].
#' @param fecundity_data See [stan_data()].
#' @param prior_data See [stan_data()].
#' @param ages See [stan_data()].
#' @param age_S_obs See [stan_data()].
#' @param age_S_eff See [stan_data()].
#' @param conditionGRonMS See [stan_data()].
#' @param model Either `"IPM"` or `"RR"`, indicating whether the data
#'   are intended for an integrated or run-reconstruction model.
#' @param life_cycle Character string indicating which life-cycle model to fit.
#'   Currently available options are spawner-to-spawner (`"SS"`, the
#'   default), spawner-to-spawner "harvest model" (`"SS_F"`), or
#'   spawner-smolt-spawner (`"SMS"`).
#' @param pool_pops Logical scalar defaulting to `TRUE`, indicating whether or not
#'   to treat the different populations as hierarchical rather than
#'   fixed/independent.
#' @param stan_model Character string giving the name of the **salmonIPM** model being
#'   fit (".stan" filetype extension is not included).
#' @param SR_fun One of `"exp"`, `"BH"` (the default), or
#'   `"Ricker"`, indicating which spawner-recruit function to fit.
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
#' @param chains Positive integer specifying the number of Markov chains.
#' @param iter Positive integer specifying the number of iterations for each
#'   chain (including warmup).
#' @param warmup Positive integer specifying the number of warmup (aka burnin)
#'   iterations per chain. If step-size adaptation is on (which it is by
#'   default), this also controls the number of iterations for which adaptation
#'   is run (and hence these warmup samples should not be used for inference).
#'   The number of warmup iterations should not be larger than `iter`
#' @param thin Positive integer specifying the period for saving samples. The
#'   default is 1, which is usually the recommended value.
#' @param cores Number of cores to use when executing the chains in parallel.
#'   Defaults to one less than the number of cores available.
#' @param ... Additional arguments to pass to [rstan::sampling()].
#' @return An object of class `stanfit` representing the fitted model. See
#'   [rstan::stan()] for details.
#'
#' @importFrom rstan stan
#'
#' @export

salmonIPM <- function(fish_data, fish_data_fwd = NULL, env_data = NULL, 
                      fecundity_data = NULL, prior_data = NULL,
                      ages = NULL, age_S_obs = NULL, age_S_eff = NULL, conditionGRonMS = FALSE,
                      model, life_cycle = "SS", pool_pops = TRUE, stan_model = NULL, SR_fun = "BH", 
                      init = NULL, pars = NULL, include = TRUE, log_lik = FALSE, 
                      chains, iter, warmup, thin = 1, cores = parallel::detectCores() - 1, ...)
{
  if(is.null(stan_model)) 
  {
    stan_model <- paste(model, life_cycle, ifelse(pool_pops, "pp", "np"), sep = "_")
  } else {
    mlp <- strsplit(stan_model, "_")[[1]]
    model <- mlp[1]
    life_cycle <- mlp[2]
    pool_pops <- mlp[3]
  }
  dat <- stan_data(fish_data = fish_data, fish_data_fwd = fish_data_fwd, env_data = env_data, 
                   fecundity_data = fecundity_data, prior_data = prior_data, 
                   ages = ages, age_S_obs = age_S_obs, age_S_eff = age_S_eff, 
                   conditionGRonMS = conditionGRonMS, stan_model = stan_model, SR_fun = SR_fun)
  
  if(is.null(pars)) 
  {
    pars <- stan_pars(stan_model)
  } else if(!include) {
    pars <- setdiff(stan_pars(stan_model), pars)
  }
  if(log_lik) pars <- c(pars, "LL")
  
  fit <- rstan::sampling(stanmodels[[stan_model]],
              data = dat, 
              init = stan_init(dat, stan_model, chains), 
              pars = pars,
              chains = chains, iter = iter, warmup = warmup, thin = thin, 
              cores = cores, ...)
}