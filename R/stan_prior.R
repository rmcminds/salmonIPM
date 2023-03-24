#' Validate prior specification for an integrated spawner-recruit model.
#' 
#' Takes an optional user-specified prior function and otherwise returns default
#' prior parameters for the given `stan_model`. This function is for internal use. 
#' 
#' @param value A user-specified prior function as described in [`priors`].
#' @param default The default prior function as described in [`priors`].
#'
#' @return A named vector of user-specifiable prior parameters whose elements are part of the
#' `data` argument passed to [rstan::sampling()] when fitting **salmonIPM** models.

stan_prior <- function(value, default)
{
  pd <- eval(default)
  if(is.null(value)) {
    return(as.array(unlist(pd[-1])))
  } else {
    p <- eval(value)
    stopifnot(length(p) == length(pd))
    stopifnot(value$dist == pd$dist)
    return(as.array(unlist(p[-1])))
  }
}
