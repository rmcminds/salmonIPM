#' Create a `draws` object from a `salmonIPMfit` object
#'
#' Convert the posterior samples in a `salmonIPMfit` object to a [draws] format
#' supported by the [posterior-package]. These S3 methods are wrappers
#' for `as_draws_*(as.array(object))`.
#'
#' @name draws
#' @aliases salmonIPMfit-draws
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
#' @details To subset variables, iterations, chains, or draws, use
#'   [subset_draws()] after making the `draws` object. 
#'
#' @return A `draws` object from the [posterior-package]. See package
#'   documentation and vignettes for details on working with these objects.
#'
#' @seealso [salmonIPM()], [salmonIPMfit], [draws]
#' @importFrom posterior as_draws
#' @exportS3Method posterior::as_draws
#' @export as_draws
as_draws.salmonIPMfit <- function(x, pars = "all", include = TRUE, ...) {
  pars <- include_pars(pars = pars, stan_model = x$stan_model, 
                       SR_fun = x$SR_fun, RRS = x$RRS, 
                       par_models = x$par_models, include = include)
  as_draws(as.array(x, pars))
}

#' @rdname draws
#' @importFrom posterior as_draws_matrix
#' @exportS3Method posterior::as_draws_matrix
#' @export as_draws_matrix
as_draws_matrix.salmonIPMfit <- function(x, pars = "all", include = TRUE, ...) {
  pars <- include_pars(pars = pars, stan_model = x$stan_model, 
                       SR_fun = x$SR_fun, RRS = x$RRS, 
                       par_models = x$par_models, include = include)
  as_draws_matrix(as.array(x, pars))
}

#' @rdname draws
#' @importFrom posterior as_draws_array
#' @exportS3Method posterior::as_draws_array
#' @export as_draws_array
as_draws_array.salmonIPMfit <- function(x, pars = "all", include = TRUE, ...) {
  pars <- include_pars(pars = pars, stan_model = x$stan_model, 
                       SR_fun = x$SR_fun, RRS = x$RRS, 
                       par_models = x$par_models, include = include)
  as_draws_array(as.array(x, pars))
}

#' @rdname draws
#' @importFrom posterior as_draws_df
#' @exportS3Method posterior::as_draws_df
#' @export as_draws_df
as_draws_df.salmonIPMfit <- function(x, pars = "all", include = TRUE, ...) {
  pars <- include_pars(pars = pars, stan_model = x$stan_model, 
                       SR_fun = x$SR_fun, RRS = x$RRS, 
                       par_models = x$par_models, include = include)
  as_draws_df(as.array(x, pars))
}

#' @rdname draws
#' @importFrom posterior as_draws_list
#' @exportS3Method posterior::as_draws_list
#' @export as_draws_list
as_draws_list.salmonIPMfit <- function(x, pars = "all", include = TRUE, ...) {
  pars <- include_pars(pars = pars, stan_model = x$stan_model, 
                       SR_fun = x$SR_fun, RRS = x$RRS, 
                       par_models = x$par_models, include = include)
  as_draws_list(as.array(x, pars))
}

#' @rdname draws
#' @importFrom posterior as_draws_rvars
#' @exportS3Method posterior::as_draws_rvars
#' @export as_draws_rvars
as_draws_rvars.salmonIPMfit <- function(x, pars = "all", include = TRUE, ...) {
  pars <- include_pars(pars = pars, stan_model = x$stan_model, 
                       SR_fun = x$SR_fun, RRS = x$RRS, 
                       par_models = x$par_models, include = include)
  as_draws_rvars(as.array(x, pars))
}

