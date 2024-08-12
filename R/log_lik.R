#' Pointwise log-likelihood matrix
#'
#' The `log_lik` method for `salmonIPMfit` objects returns the pointwise
#' log-likelihood matrix.
#'
#' @name log_lik.salmonIPMfit
#' @aliases log_lik
#'
#' @param object An object of class [salmonIPMfit] that includes samples of the
#'   pointwise log-likelihood `LL`, e.g. `log_lik` was set to `TRUE` when
#'   calling [salmonIPM()].
#' @param ... Currently ignored.
#'
#' @return A `N_samples` by `N` matrix containing the pointwise log-likelihood,
#'   where `N_samples` is the number of post-warmup posterior draws and `N` is
#'   the number of cases (i.e. rows) in `fish_data`. Note that this includes all
#'   likelihood components for observations in `fish_data`, but does not include
#'   the likelihood of auxiliary data types. These, as well as the likelihood
#'   components for specific data types in `fish_data`, can be extracted with
#'   [as.matrix()] if they have been monitored using the `pars` argument to
#'   [salmonIPM()].
#'
#' @seealso [salmonIPM()]
#' @importFrom rstantools log_lik
#' @exportS3Method rstantools::log_lik
#' @export log_lik
log_lik.salmonIPMfit <- function(object, ...) {
  as.matrix(object, pars = "LL")
}

#' @rdname log_lik.salmonIPMfit
#' @importFrom rstantools nsamples
#' @exportS3Method rstantools::nsamples
#' @export nsamples
nsamples.salmonIPMfit <- function(object, ...) {
  sum(object$stanfit@sim$n_save - object$stanfit@sim$warmup2)
}
