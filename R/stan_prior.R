#' Validate prior specification for an integrated spawner-recruit model.
#' 
#' Takes an optional user-specified prior function and otherwise returns default
#' prior parameters for the given `stan_model`. This function is for internal use. 
#' 
#' @param value A user-specified prior function as described in [`priors`].
#' @param default The default prior function as described in [`priors`].
#'
#' @return A named array of user-specifiable prior parameters whose elements are part of the
#' `data` argument passed to [rstan::sampling()] when fitting **salmonIPM** models.
#' The `dist` attribute is a character string giving the name of the distribution.

stan_prior <- function(value, default)
{
  pd <- eval(default)
  if(is.null(value)) {
    p <- pd
  } else {
    p <- eval(value)
    stopifnot(all(names(p) == names(pd)))
    stopifnot(value$dist == pd$dist)
  }
  
  return(structure(as.array(unlist(p[-1])), dist = p$dist))
}
