#' Fits an integrated or run-reconstruction spawner-recruit model.
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
#'   matrices defined in the Stan model being used. (This is required if
#'   \code{life_cycle != "SS"}.)
#' @param ages If \code{life_cycle != "SS"}, a named list giving the fixed ages
#'   in years of all subadult life stages.
#' @param age_S_obs Only if `stan_model %in% c("IPM_SS_np","IPM_SS_pp")`, a logical or numeric
#'   vector indicating, for each adult age, whether observed total spawner data
#'   includes that age. The default is to treat `S_obs` as including spawners of all ages.
#' @param model Either \code{"IPM"} or \code{"RR"}, indicating whether the data
#'   are intended for an integrated or run-reconstruction model.
#' @param life_cycle Character string indicating which life-cycle model to fit.
#'   Currently available options are spawner-to-spawner (\code{"SS"}, the
#'   default), spawner-to-spawner "harvest model" (\code{"SS_F"}), or
#'   spawner-smolt-spawner (\code{"SMS"}).
#' @param pool_pops Logical, with default \code{TRUE}, indicating whether or not
#'   to treat the different populations as hierarchical rather than
#'   fixed/independent. Must be TRUE if model == "IPM_F".
#' @param stan_model Character string giving the name of the Stan model being
#'   fit (".stan" filetype extension is not included). If provided,
#'   \code{"stan_model"} overrides \code{"model"}, \code{"life_cycle"}, and
#'   \code{"pool_pops"}.
#' @param SR_fun One of \code{"exp"}, \code{"BH"} (the default), or
#'   \code{"Ricker"}, indicating which spawner-recruit function to fit.
#' @param init A list of named lists of initial values to be passed to
#'   \code{rstan::stan}. If \code{NULL}, initial values will be automatically
#'   generated from the supplied data using \code{salmonIPM::stan_init}.
#' @param pars A vector of character strings specifying parameters to monitor.
#'   If NULL, default values are used. If a non-default value is supplied, it is
#'   the user's responsibility to make sure the parameters requested appear in
#'   the model configuration specified.
#' @param log_lik A logical indicator as to whether the pointwise log-likelihood
#'   should be returned for later analysis with \code{loo:loo}.
#' @param chains A positive integer specifying the number of Markov chains.
#' @param iter A positive integer specifying the number of iterations for each
#'   chain (including warmup).
#' @param warmup A positive integer specifying the number of warmup (aka burnin)
#'   iterations per chain. If step-size adaptation is on (which it is by
#'   default), this also controls the number of iterations for which adaptation
#'   is run (and hence these warmup samples should not be used for inference).
#'   The number of warmup iterations should not be larger than \code{iter}
#' @param thin A positive integer specifying the period for saving samples. The
#'   default is 1, which is usually the recommended value.
#' @param cores Number of cores to use when executing the chains in parallel.
#'   Defaults to one less than the number of cores available.
#' @param ... Additional arguments to pass to \code{stan}.
#' @return An object of class \code{stanfit} representing the fitted model. See
#'   \code{rstan::stan} for details.
#'
#' @importFrom rstan stan
#'
#' @export

salmonIPM <- function(fish_data, fish_data_fwd = NULL, env_data = NULL, ages = NULL, age_S_obs = NULL,  
                      model, life_cycle = "SS", pool_pops = TRUE, stan_model = NULL, SR_fun = "BH", 
                      init = NULL, pars = NULL, log_lik = FALSE, 
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
  dat <- stan_data(fish_data, fish_data_fwd, env_data, catch_data, ages, age_S_obs, stan_model, SR_fun)
  
  if(is.null(pars)) pars <- stan_pars(stan_model)
  if(log_lik) pars <- c(pars, "LL")
  
  fit <- stan(file = file.path(path.package("salmonIPM"), "stan", paste0(stan_model, ".stan")),
              data = dat, 
              init = stan_init(dat, stan_model, chains), 
              pars = pars,
              chains = chains, iter = iter, warmup = warmup, thin = thin, 
              cores = cores, ...)
}