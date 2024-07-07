#' Extract posterior means from `salmonIPMfit` object
#'
#' Convenience function to simplify extracting posterior means or single
#' parameters from `stanfit` objects.
#' 
#' @param object A fitted [salmonIPMfit] object.
#' @param pars Character vector with the names of parameters to summarize.
#' 
#' @return A scalar or vector of means of the posterior distribution of `pars`.
#' 
#' @importFrom rstan get_posterior_mean
#'
#' @export
stan_mean <- function(object, pars)
{
  mm <- get_posterior_mean(object$stanfit, pars)
  return(mm[,ncol(mm)])
}


#' Extract posterior samples for a single parameter
#'
#' Wrapper for [extract()] applied to the `stanfit` component of a fitted [salmonIPMfit] object.
#' Unlike [extract()], only returns the first element, not a named list. This is useful for
#' obtaining the samples for a single parameter (of any dimension), though usually equivalent
#' to `as.matrix(object, par)`.
#' 
#' @param object An object of class [salmonIPMfit].
#' @param par Character string giving quantity to return.
#' 
#' @importFrom rstan extract
#'
#' @export
extract1 <- function(object, par)
{
  rstan::extract(object, par)[[1]]
}
