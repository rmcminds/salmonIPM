#' Create a `draws` object from a `salmonIPMfit` object
#'
#' Convert the posterior samples in a `salmonIPMfit` object to a [draws] format
#' supported by the [posterior-package]. These S3 methods are wrappers
#' for `as_draws_*(as.array(object))`.
#'
#' @name salmonIPMfit-draws
#'
#' @param object An object of class [salmonIPMfit].
#' @param ... Additional arguments passed to [as.array.stanfit()]. Currently
#'   `pars` is the only argument supported.
#'
#' @details To subset variables, iterations, chains, or draws, use
#'   [subset_draws()] after making the `draws` object. Variables can also be
#'   selected by using `...` to pass the `pars` argument to `as.array.stanfit()`.
#'
#' @return A `draws` object from the [posterior-package]. See package
#'   documentation and vignettes for details on working with these objects.
#'
#' @seealso [salmonIPM()], [salmonIPMfit], [draws]

#' @rdname salmonIPMfit-draws
#' @method as_draws salmonIPMfit
#' @importFrom posterior as_draws
#' @exportS3Method posterior::as_draws
#' @export
as_draws.salmonIPMfit <- function(object, ...) {
  as_draws(as.array(object, ...))
}

#' @rdname salmonIPMfit-draws
#' @method as_draws_matrix salmonIPMfit
#' @importFrom posterior as_draws_matrix
#' @exportS3Method posterior::as_draws_matrix
#' @export
as_draws_matrix.salmonIPMfit <- function(object, ...) {
  as_draws_matrix(as.array(object, ...))
}

#' @rdname salmonIPMfit-draws
#' @method as_draws_array salmonIPMfit
#' @importFrom posterior as_draws_array
#' @exportS3Method posterior::as_draws_array
#' @export
as_draws_array.salmonIPMfit <- function(object, ...) {
  as_draws_array(as.array(object, ...))
}

#' @rdname salmonIPMfit-draws
#' @method as_draws_df salmonIPMfit
#' @importFrom posterior as_draws_df
#' @exportS3Method posterior::as_draws_df
#' @export
as_draws_df.salmonIPMfit <- function(object, ...) {
  as_draws_df(as.array(object$stanfit, ...))
}

#' @rdname salmonIPMfit-draws
#' @method as_draws_list salmonIPMfit
#' @importFrom posterior as_draws_list
#' @exportS3Method posterior::as_draws_list
#' @export
as_draws_list.salmonIPMfit <- function(object, ...) {
  as_draws_list(as.array(object$stanfit, ...))
}

#' @rdname salmonIPMfit-draws
#' @method as_draws_rvars salmonIPMfit
#' @importFrom posterior as_draws_rvars
#' @exportS3Method posterior::as_draws_rvars
#' @export
as_draws_rvars.salmonIPMfit <- function(object, ...) {
  as_draws_rvars(as.array(object$stanfit, ...))
}
