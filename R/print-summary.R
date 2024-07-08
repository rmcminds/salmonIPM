#' Print method for `salmonIPMfit` objects
#'
#' Print basic information on the fitted model and a posterior summary for
#' parameters of interest estimated by the samples included in a [salmonIPMfit]
#' object. This is currently a wrapper for [print.stanfit()] that prepends a
#' summary of the specified **salmonIPM** model to the output of
#' `print.stanfit(object$stanfit)`.
#'
#' @name print.salmonIPMfit
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

#' @rdname print.salmonIPMfit
#' @method print salmonIPMfit
#' @export
print.salmonIPMfit <- function(object, pars = "hyper", include = TRUE, 
                               probs = c(0.05, 0.5, 0.95), digits_summary = 2, ...) 
{
  stan_model <- object$stan_model
  SR_fun <- object$SR_fun
  RRS <- object$RRS
  N_pop <- object$dims$N_pop
  
  if(all(pars %in% c("all","hyper","group","states","ppd")))
    pars <- stan_pars(stan_model, pars = pars, SR_fun = SR_fun, RRS = RRS)
  if(!include) 
    pars <- setdiff(stan_pars(stan_model, pars = "all", SR_fun = SR_fun, RRS = RRS), pars)
  
  cat(
    paste0(
      "salmonIPMfit\n",
      " model type: ", object$model, "\n",
      " life cycle: ", object$life_cycle, "\n",
      " SR_fun: ", object$SR_fun, " (RRS: ", RRS, ")\n",
      " N = ", object$dims$N, " cases in fish_data\n",
      " N_pop = ", N_pop, 
      if(N_pop > 1) paste0(" (", ifelse(object$pool_pops, "partial", "no"), " pooling)"),
      "\n N_year = ", object$dims$N_year, "\n\n"
    )
  )
  
  print(object$stanfit, pars = pars, probs = probs, digits_summary = digits_summary, ...)
}




