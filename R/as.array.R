#' Extract posterior samples as a matrix, array or data frame
#' 
#' Post-warmup samples from a [salmonIPMfit] object are converted into standard formats.
#' These S3 methods are wrappers for `as.*.stanfit(object$stanfit)`; see [as.array.stanfit()]
#' for details.
#' 
#' @name as.array.salmonIPMfit
#' 
#' @param x An object of class [salmonIPMfit].
#' @param pars A character vector specifying (hyper)parameters, states, and/or
#'   quantities of interest ("parameters"). Parameters can be explicitly named
#'   or one or more shortcuts can be used to specify hierarchical levels of
#'   parameters; see [stan_pars()] for details. The default is `"all"`
#'   quantities that were monitored.
#' @param include Logical scalar defaulting to `TRUE` indicating whether to
#'   include or exclude the parameters given by `pars`.
#' @param ... Currently ignored.
#' 
#' @return A `matrix`, `array` or `data.frame` containing the post-warmup samples of
#' the specified quantities.
#' 
#' @seealso [salmonIPM()], [salmonIPMfit], [as.array.stanfit()]
#' @method as.matrix salmonIPMfit
#' @export
as.matrix.salmonIPMfit <- function(x, pars = "all", include = TRUE, ...) {
  pars <- include_pars(pars = pars, stan_model = x$stan_model, 
                       SR_fun = x$SR_fun, RRS = x$RRS, 
                       par_models = x$par_models, include = include)
  as.matrix(x$stanfit, pars = pars, ...)
}

#' @rdname as.array.salmonIPMfit
#' @method as.array salmonIPMfit
#' @export
as.array.salmonIPMfit <- function(x, pars = "all", include = TRUE, ...) {
  pars <- include_pars(pars = pars, stan_model = x$stan_model, 
                       SR_fun = x$SR_fun, RRS = x$RRS, 
                       par_models = x$par_models, include = include)
  as.array(x$stanfit, pars = pars, ...)
}

#' @rdname as.array.salmonIPMfit
#' @method as.data.frame salmonIPMfit
#' @export
as.data.frame.salmonIPMfit <- function(x, pars = "all", include = TRUE, ...) {
  pars <- include_pars(pars = pars, stan_model = x$stan_model, 
                       SR_fun = x$SR_fun, RRS = x$RRS, 
                       par_models = x$par_models, include = include)
  as.data.frame(x$stanfit, pars = pars, ...)
}

