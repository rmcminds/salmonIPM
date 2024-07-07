#' Extract posterior samples as a matrix, array or data frame
#' 
#' Post-warmup samples from a [salmonIPMfit] object are converted into standard formats.
#' These S3 methods are wrappers for `as.*.stanfit(object$stanfit)`; see [as.array.stanfit()]
#' for details.
#' 
#' @name as.array.salmonIPMfit
#' 
#' @param object An object of class [salmonIPMfit].
#' @param ... Additional arguments that can be passed to [extract()]. 
#' Currently `pars` is the only argument supported.
#' 
#' @return A `matrix`, `array` or `data.frame` containing the post-warmup samples of
#' the specified quantities.
#' 
#' @seealso [salmonIPM()], [salmonIPMfit], [as.array.stanfit()]
#' 

#' @method as.matrix salmonIPMfit
#' @rdname as.array.salmonIPMfit
#' @export
as.matrix.salmonIPMfit <- function(object, ...) {
  rstan:::as.matrix.stanfit(object$stanfit, ...)
}

#' @method as.array salmonIPMfit
#' @rdname as.array.salmonIPMfit
#' @export
as.array.salmonIPMfit <- function(object, ...) {
  rstan:::as.array.stanfit(object$stanfit, ...)
}

#' @method as.data.frame salmonIPMfit
#' @rdname as.array.salmonIPMfit
#' @export
as.data.frame.salmonIPMfit <- function(object, ...) {
  rstan:::as.data.frame.stanfit(object$stanfit, ...)
}

