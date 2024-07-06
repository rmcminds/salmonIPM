#' Create a salmonIPMfit object
#' 
#' @param stanfit An object of class [rstan::`stanfit-class`].
#' @param call An object of class `[methods::call-class]` representing the call to [salmonIPM()].
#' @inheritParams salmonIPM
#' 
#' @return A salmonIPMfit object.

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



