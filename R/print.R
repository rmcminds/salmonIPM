#' Print method for `salmonIPMfit` objects
#'
#' Print information on the fitted model and a posterior summary for parameters
#' of interest estimated by the samples included in a [salmonIPMfit] object.
#'
#' @name print.salmonIPMfit
#' @method print salmonIPMfit
#'
#' @param x An object of class [salmonIPMfit].
#' @param pars A character vector specifying (hyper)parameters, states, and/or
#'   quantities of interest ("parameters") to summarize. Parameters can be
#'   explicitly named or one or more shortcuts can be used to specify
#'   hierarchical levels of parameters; see [stan_pars()] for details. The
#'   default is `"hyper"`, i.e. the top-level hyperparameters that are given
#'   priors.
#' @param include Logical scalar defaulting to `TRUE` indicating whether to
#'   include or exclude the parameters given by `pars`. If `FALSE`, only entire
#'   multidimensional parameters can be excluded, rather than particular
#'   elements of them.
#' @param probs A numeric vector of posterior quantiles to print. The default is
#'   `c(0.05, 0.5, 0.95)`, i.e. the median and 90% credible interval.
#' @param digits Number of decimal places to print, defaulting to 2. Applies to
#'   quantities other than the effective sample size, which is always rounded to
#'   the nearest integer.
#' @param ... Currently ignored.
#'
#' @details If `pars` includes correlation matrices, only the lower triangular
#'   elements are returned. This avoids redundant summary output as well as
#'   false positive diagnostic results such as `bulk-ESS` and `Rhat` being `NaN`
#'   for the diagonal elements.
#'
#' @examples
#' # <under construction>
#'
#' @seealso [salmonIPMfit], [diagnostics]
#'
#' @importFrom rstan get_num_divergent get_num_max_treedepth get_bfmi
#' @importFrom rstan check_divergences check_treedepth check_energy
#' @export

print.salmonIPMfit <- function(x, pars = "hyper", include = TRUE, 
                               probs = c(0.05, 0.5, 0.95), digits = 2, ...) 
{
  # Model info
  stan_model <- x$stan_model
  SR_fun <- x$SR_fun
  RRS <- x$RRS
  RRS_out <- paste(RRS, collapse = ", ")
  N_pop <- x$dims$N_pop
  ages <- paste(paste(names(x$ages), x$ages, sep = " = "), collapse = ", ")
  par_models <- x$par_models
  par_models_out <- paste(par_models, collapse = ", ")
  
  cat(
    paste0(
      "salmonIPMfit\n",
      " model type: ", x$model, "\n",
      " life cycle: ", x$life_cycle, 
      if(nchar(ages)) paste0(" (ages: ", ages, ")"),
      "\n SR_fun: ", x$SR_fun, " (RRS: ", RRS_out, ")\n",
      if(nchar(par_models_out)) paste0(" parameter models: ", par_models_out, "\n"),
      " N = ", x$dims$N, " cases in fish_data\n",
      " N_pop = ", N_pop, 
      if(N_pop > 1) paste0(" (", ifelse(x$pool_pops, "partial", "no"), " pooling)"),
      "\n N_year = ", x$dims$N_year, "\n\n"
    )
  )
  
  # Posterior summary
  sdf <- summary(x, pars = pars, include = include, probs = probs)
  rownames(sdf) <- sdf$variable
  sdf <- round(sdf[,-1], digits = digits)
  sdf$ess_bulk <- round(sdf$ess_bulk)
  names(sdf) <- gsub("ess_bulk", "bulk-ESS", names(sdf)) 
  names(sdf) <- gsub("rhat", "Rhat", names(sdf))
  
  print(sdf)
  
  # Sampling info
  chains <- x$stanfit@sim$chains
  iter <- x$stanfit@sim$iter
  warmup <- x$stanfit@sim$warmup
  thin <- x$stanfit@sim$thin
  
  cat(
    paste0(
      "\n", nsamples(x), " posterior draws (chains = ", chains, ", iter = ", iter, 
      ", warmup = ", warmup, ", thin = ", thin, ")\n",
      "Bulk-ESS is a crude measure of effective sample size for location summaries.\n",
      "Rhat is the potential scale reduction factor based on split chains; \n",
      " at convergence Rhat = 1. See ?posterior::diagnostics for details.\n"
    )
  )
  
  if(get_num_divergent(x$stanfit)) {
    check_divergences(x$stanfit)
    message("See ", "https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup")
  }
  if(get_num_max_treedepth(x$stanfit)) {
    check_treedepth(x$stanfit)
    message("See ", "https://mc-stan.org/misc/warnings.html#maximum-treedepth")
  }
  if(any(get_bfmi(x$stanfit) < 0.2)) {
    check_energy(x$stanfit)
    message("See ", "https://mc-stan.org/misc/warnings.html#bfmi-low")
  }
}
