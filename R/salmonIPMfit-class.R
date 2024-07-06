#' Fitted model objects: the `salmonIPMfit` class
#' 
#' Models fitted with [salmonIPM-package] are objects of S3 class `salmonIPMfit`. Elements 
#' generally mirror the named arguments passed to the constructor function by [salmonIPM()]
#' and are mostly self-explanatory, except `prior.info` is a named list whose elements are 
#' lists returned by [priors] functions for all hyperparameters in the model, whether
#' modifiable (user-specified or default) or hard-coded.
#' 
#' @param stanfit An object of class `stanfit`.
#' @param call An object of class [call-class] containing the call to [salmonIPM()].
#' @inheritParams salmonIPM
#' 
#' @return A salmonIPMfit object.
#' 
#' @name salmonIPMfit-class

salmonIPMfit <- function(stanfit, call, model, life_cycle, pool_pops, SR_fun, RRS, 
                         par_models, center, scale, prior, age_S_obs, age_S_eff,
                         conditionGRonMS)
{
  out <- list(stanfit = stanfit, call = call, model = model, life_cycle = life_cycle,
              pool_pops = pool_pops, SR_fun = SR_fun, RRS = RRS,
              par_models = par_models, center = center, scale = scale,
              prior = prior, age_S_obs = age_S_obs, age_S_eff = age_S_eff,
              conditionGRonMS = conditionGRonMS)
  
  structure(out, class = c("salmonIPMfit", "list"))
}
