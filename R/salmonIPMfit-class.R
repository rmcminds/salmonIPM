#' Fitted model objects: the `salmonIPMfit` class
#'
#' Models fitted with [salmonIPM-package] are objects of S3 class
#' `salmonIPMfit`. Elements generally mirror the named arguments passed to the
#' constructor function by [salmonIPM()] and are mostly self-explanatory, except
#' `prior.info` is a named list whose elements are lists returned by [priors]
#' functions for all hyperparameters in the model, whether modifiable
#' (user-specified or default) or hard-coded.
#'
#' @param stanfit An object of class [stanfit-class].
#' @param call An object of class [call-class] containing the call to
#'   [salmonIPM()].
#' @param prior.info A named list of priors on `pars`. Each element is itself a
#'   named list representing the prior in the format returned by the functions
#'   in [priors] with an additional attribute `"type"`.
#' @param stan_data Optional named list of input data passed to
#'   [rstan::sampling()] for fitting an IPM, as returned by [stan_data()].
#'   Stored if [salmonIPM()] is called with `save_data = TRUE`, otherwise
#'   `NULL`.
#' @param dims A named list of data dimensions, including the number of cases
#'   (i.e. rows in `fish_data`) `N`, the number of populations `N_pop`, the
#'   number of years `N_year`, and any other model-specific elements.
#' @param pops Character vector giving the unique elements of `fish_data$pop`.
#' @param elapsed_time Wall time (s) to fit the model, as returned by
#'   [get_elapsed_time(stanfit)].
#' @param salmonIPM_version The version of **salmonIPM** used to fit the model.
#' @inheritParams salmonIPM
#'
#' @return A `salmonIPMfit` object with elements described above.
#'
#' @name salmonIPMfit-class
#' @seealso [salmonIPM()], [priors]
salmonIPMfit <- function(stanfit, call, stan_model, model, life_cycle, ages, pool_pops, 
                         SR_fun, RRS, par_models, center, scale, age_S_obs, age_S_eff,
                         conditionGRonMS, prior.info, stan_data, dims, pops, elapsed_time)
{
  out <- list(stanfit = stanfit, call = call, 
              stan_model = stan_model, model = model, life_cycle = life_cycle,
              ages = ages, pool_pops = pool_pops, SR_fun = SR_fun, RRS = RRS,
              par_models = par_models, center = center, scale = scale,
              age_S_obs = age_S_obs, age_S_eff = age_S_eff,
              conditionGRonMS = conditionGRonMS, prior.info = prior.info, 
              stan_data = stan_data, dims = dims, pops = pops, elapsed_time = elapsed_time)
  
  structure(out, class = c("salmonIPMfit", "list"))
}
