#' Information criteria and cross-validation
#'
#' Compute approximate leave-one-out cross-validation (LOO, LOOIC) scores using
#' the [loo-package]. See [loo::loo()] for more information.
#'
#' @name loo.salmonIPMfit
#' @aliases loo
#'
#' @param x A fitted model object of class [salmonIPMfit].
#' @param ... Additional arguments passed to the S3 generic.
#'
#' @return The structure of the `loo` objects returned is described in
#'   [loo::loo()].
#'   
#' @importFrom loo loo relative_eff
#' @exportS3Method loo::loo
#' @export loo
loo.salmonIPMfit <- function(x, ...) {
  LL <- as.array(x, pars = "LL")
  L <- exp(LL)
  r_eff <- relative_eff(L)
  loo(LL, r_eff = r_eff, ...)
}
