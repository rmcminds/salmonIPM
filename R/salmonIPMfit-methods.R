#' Extract posterior samples as a matrix, array or data frame
#' 
#' Post-warmup samples from a [salmonIPMfit] object are converted into standard formats.
#' These S3 methods are wrappers for `as.*.stanfit(object$stanfit)`; see [as.array.stanfit()]
#' for details.
#' 
#' @param object An object of class [salmonIPMfit].
#' @param ... Additional arguments that can be passed to [extract()]. 
#' Currently `pars` is the only argument supported.
#' 
#' @return A `matrix`, `array` or `data.frame` containing the post-warmup samples of
#' the specified quantities.
#' 
#' @seealso [salmonIPM()], [salmonIPMfit], [as.array.stanfit()]
#' @name as.array.salmonIPMfit
#' 
#' @details
#' methods(class = 'stanfit')
#' 
#' extract, get_posterior_mean (S4 generics: rstan -> use extract1, stan_mean instead)    
#' loo, loo_moment_match (S4 generics: rstan -> use S3 generics from rstanarm instead)
#' print
#' summary          
#'
#' methods(class = 'stanreg')
#' 
#' launch_shinystan (S3 generic: shinystan)
#' as_draws, as_draws_array, as_draws_df, as_draws_list, as_draws_matrix, as_draws_rvars (S3 generics: posterior)        
#' log_lik
#' loo, loo_compare (S3 generics: loo)
#' posterior_predict (S3 generic: rstantools)   
#' prior_summary (S3 generic: rstantools)
#' print                  
#' summary

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

