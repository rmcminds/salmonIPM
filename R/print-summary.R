#' Print method for `salmonIPMfit` objects
#'
#' Print basic information on the fitted model and a posterior summary for
#' parameters of interest estimated by the samples included in a [salmonIPMfit]
#' object. This is currently a wrapper for [print.stanfit()] that prepends a
#' summary of the specified **salmonIPM** model to the output of
#' `print.stanfit(object$stanfit)`.
#'
#' @name print.salmonIPMfit
#' @method print salmonIPMfit
#'
#' @param object An object of class [salmonIPMfit].
#' @param pars A character vector specifying (hyper)parameters, states, and/or
#'   quantities of interest ("parameters") to summarize. Parameters can be
#'   explicitly named or one or more shortcuts can be used to specify
#'   hierarchical levels of parameters; see [stan_pars()] for details. The
#'   default is `"hyper"`, i.e. the top-level hyperparameters that are given priors.
#' @param include Logical scalar defaulting to `TRUE` indicating whether to
#'   include or exclude the parameters given by `pars`. If `FALSE`, only entire
#'   multidimensional parameters can be excluded, rather than particular
#'   elements of them.
#' @param probs A numeric vector of posterior quantiles to print. Unlike
#'   [print.stanfit()], the default is `c(0.05,0.5,0.95)`, i.e. the median and
#'   90% credible interval.
#' @param digits_summary Number of significant digits to print, defaulting to 2.
#'   Applies to quantities other than the effective sample size, which is always
#'   rounded to the nearest integer.
#' @param ... Additional arguments passed to [summary.stanfit()].
#'
#' @seealso [print.stanfit()], [salmonIPMfit]
#' @export

print.salmonIPMfit <- function(object, pars = "hyper", include = TRUE, 
                               probs = c(0.05, 0.5, 0.95), digits_summary = 2, ...) 
{
  stan_model <- object$stan_model
  SR_fun <- object$SR_fun
  RRS <- object$RRS
  RRS_out <- paste(RRS, collapse = ", ")
  N_pop <- object$dims$N_pop
  ages <- paste(paste(names(object$ages), object$ages, sep = " = "), collapse = ", ")
  par_models <- object$par_models
  par_models_out <- paste(par_models, collapse = ", ")

  if(all(pars %in% c("all","hyper","group","states","ppd")))
    pars <- stan_pars(stan_model, pars = pars, SR_fun = SR_fun, 
                      RRS = RRS, par_models = par_models)
  if(!include) 
    pars <- setdiff(stan_pars(stan_model, pars = "all", SR_fun = SR_fun, 
                              RRS = RRS, par_models = par_models), 
                    pars)
  
  cat(
    paste0(
      "salmonIPMfit\n",
      " model type: ", object$model, "\n",
      " life cycle: ", object$life_cycle, 
      if(nchar(ages)) paste0(" (ages: ", ages, ")"),
      "\n SR_fun: ", object$SR_fun, " (RRS: ", RRS_out, ")\n",
      if(nchar(par_models_out)) paste0(" parameter models: ", par_models_out, "\n"),
      " N = ", object$dims$N, " cases in fish_data\n",
      " N_pop = ", N_pop, 
      if(N_pop > 1) paste0(" (", ifelse(object$pool_pops, "partial", "no"), " pooling)"),
      "\n N_year = ", object$dims$N_year, "\n\n"
    )
  )
  
  print(object$stanfit, pars = pars, probs = probs, digits_summary = digits_summary, ...)
}


#' Summary method for `salmonIPMfit` objects
#'
#' Summarize the distributions of estimated parameters and derived quantities
#' using the posterior draws in a [salmonIPMfit] object. This is a wrapper for
#' `rstan::summary(object$stanfit)`.
#'
#' @name summary.salmonIPMfit
#' @method summary salmonIPMfit
#'
#' @param object An object of class [salmonIPMfit].
#' @param pars A character vector of parameter names. Defaults to all monitored
#'   parameters as well as the log-posterior (`lp__`).
#' @param probs A numeric vector of posterior quantiles. The default is
#'   `c(0.025,0.25,0.5,0.75,0.975)`.
#' @param use_cache Logical, defaulting to `TRUE`. When `use_cache = TRUE` the summary
#'   quantities for all parameters are computed and cached for future use.
#'   Setting `use_cache = FALSE` can be used to avoid performing the summary
#'   computations for all parameters if `pars` is given as some specific
#'   parameters.
#'   
#' @seealso [rstan::summary,stanfit-method], [salmonIPMfit]
#' @importFrom rstan summary
#' @export

summary.salmonIPMfit <- function(object, pars, probs = c(0.025, 0.25, 0.50, 0.75, 0.975), 
                                 use_cache = TRUE) 
{
  summary(object$stanfit, pars = pars, probs = probs, use_cache = use_cache)
}

