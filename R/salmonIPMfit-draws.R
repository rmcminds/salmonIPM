#' Create a `draws` object from a `salmonIPMfit` object
#'
#' Convert the posterior samples in a `salmonIPMfit` object to a [draws] format
#' supported by the [posterior-package]. These S3 methods are wrappers
#' for `as_draws_*(as.array(object))`.
#'
#' @name salmonIPMfit-draws
#' @aliases as_draws
#'
#' @param object An object of class [salmonIPMfit].
#' @param pars A character vector specifying (hyper)parameters, states, and/or
#'   quantities of interest ("parameters"). Parameters can be explicitly named
#'   or one or more shortcuts can be used to specify hierarchical levels of
#'   parameters; see [stan_pars()] for details. The default is `"all"`
#'   quantities that were monitored.
#' @param include Logical scalar defaulting to `TRUE` indicating whether to
#'   include or exclude the parameters given by `pars`.
#'
#' @details To subset variables, iterations, chains, or draws, use
#'   [subset_draws()] after making the `draws` object. 
#'
#' @return A `draws` object from the [posterior-package]. See package
#'   documentation and vignettes for details on working with these objects.
#'
#' @seealso [salmonIPM()], [salmonIPMfit], [draws]

#' @rdname salmonIPMfit-draws
#' @importFrom posterior as_draws
#' @exportS3Method posterior::as_draws
#' @export as_draws
as_draws.salmonIPMfit <- function(object, pars = "all", include = TRUE) {
  pars <- include_pars(pars = pars, stan_model = object$stan_model, 
                       SR_fun = object$SR_fun, RRS = object$RRS, 
                       par_models = object$par_models, include = include)
  as_draws(as.array(object, pars))
}

#' @rdname salmonIPMfit-draws
#' @importFrom posterior as_draws_matrix
#' @exportS3Method posterior::as_draws_matrix
#' @export as_draws_matrix
as_draws_matrix.salmonIPMfit <- function(object, pars = "all", include = TRUE) {
  pars <- include_pars(pars = pars, stan_model = object$stan_model, 
                       SR_fun = object$SR_fun, RRS = object$RRS, 
                       par_models = object$par_models, include = include)
  as_draws_matrix(as.array(object, pars))
}

#' @rdname salmonIPMfit-draws
#' @importFrom posterior as_draws_array
#' @exportS3Method posterior::as_draws_array
#' @export as_draws_array
as_draws_array.salmonIPMfit <- function(object, pars = "all", include = TRUE) {
  pars <- include_pars(pars = pars, stan_model = object$stan_model, 
                       SR_fun = object$SR_fun, RRS = object$RRS, 
                       par_models = object$par_models, include = include)
  as_draws_array(as.array(object, pars))
}

#' @rdname salmonIPMfit-draws
#' @importFrom posterior as_draws_df
#' @exportS3Method posterior::as_draws_df
#' @export as_draws_df
as_draws_df.salmonIPMfit <- function(object, pars = "all", include = TRUE) {
  pars <- include_pars(pars = pars, stan_model = object$stan_model, 
                       SR_fun = object$SR_fun, RRS = object$RRS, 
                       par_models = object$par_models, include = include)
  as_draws_df(as.array(object, pars))
}

#' @rdname salmonIPMfit-draws
#' @importFrom posterior as_draws_list
#' @exportS3Method posterior::as_draws_list
#' @export as_draws_list
as_draws_list.salmonIPMfit <- function(object, pars = "all", include = TRUE) {
  pars <- include_pars(pars = pars, stan_model = object$stan_model, 
                       SR_fun = object$SR_fun, RRS = object$RRS, 
                       par_models = object$par_models, include = include)
  as_draws_list(as.array(object, pars))
}

#' @rdname salmonIPMfit-draws
#' @importFrom posterior as_draws_rvars
#' @exportS3Method posterior::as_draws_rvars
#' @export as_draws_rvars
as_draws_rvars.salmonIPMfit <- function(object, pars = "all", include = TRUE) {
  pars <- include_pars(pars = pars, stan_model = object$stan_model, 
                       SR_fun = object$SR_fun, RRS = object$RRS, 
                       par_models = object$par_models, include = include)
  as_draws_rvars(as.array(object, pars))
}

#' Select parameters to include
#'
#' Helper function to amend a vector of parameter names or hierarchical level
#' shortcuts.
#'
#' @param pars A character vector specifying (hyper)parameters, states, and/or
#'   quantities of interest ("parameters"). Parameters can be explicitly named
#'   or one or more shortcuts can be used to specify hierarchical levels of
#'   parameters; see [stan_pars()] for details. 
#' @param include Logical scalar defaulting to `TRUE` indicating whether to
#'   include or exclude the parameters given by `pars`. 
#' @inheritParams salmonIPM
#' @return A character vector of amended pars.
include_pars <- function(pars, stan_model, SR_fun, RRS, par_models, 
                         include, log_lik = FALSE) {
  if(all(pars %in% c("all","hyper","group","states","ppd")))
    pars <- stan_pars(stan_model, pars = pars, SR_fun = SR_fun, 
                      RRS = RRS, par_models = par_models)
  if(!include) 
    pars <- setdiff(stan_pars(stan_model, pars = "all", SR_fun = SR_fun, 
                              RRS = RRS, par_models = par_models), 
                    pars)
  if(log_lik) pars <- c(pars, "LL")
  return(pars)
}
