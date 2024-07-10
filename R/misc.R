#' Extract posterior means from `salmonIPMfit` object
#'
#' Convenience function to simplify extracting posterior means or single
#' parameters from `stanfit` objects.
#' 
#' @param object A fitted [salmonIPMfit] object.
#' @param pars Character vector with the names of parameters to summarize.
#' @return A scalar or vector of means of the posterior distribution of `pars`.
#' @importFrom rstan get_posterior_mean
#' @export
stan_mean <- function(object, pars) {
  mm <- get_posterior_mean(object$stanfit, pars)
  return(mm[,ncol(mm)])
}


#' Extract posterior samples for a single parameter
#'
#' Wrapper for [extract()] applied to the `stanfit` component of a fitted [salmonIPMfit] object.
#' Unlike [extract()], only returns the first element, not a named list. This is useful for
#' obtaining the samples for a single parameter (of any dimension), though usually equivalent
#' to `as.matrix(object, par)`.
#' 
#' @param object An object of class [salmonIPMfit].
#' @param par Character string giving quantity to return.
#' @importFrom rstan extract
#' @export
extract1 <- function(object, par) {
  extract(object$stanfit, par)[[1]]
}

#' Validate RRS specification
#'
#' Check whether an RRS specification is consistent with a stan_model.
#'
#' @inheritParams salmonIPM
#' @return Nothing is returned; the function throws an error if the requested
#'   `RRS` is not found in the hyperparameters of `stan_model`.
validate_RRS <- function(stan_model, SR_fun = "BH", RRS) {
  pool_pops <- switch(strsplit(stan_model, "_")[[1]][3], np = FALSE, pp = TRUE)
  stanmodel <- gsub("iter", "", stan_model) # same Stan code for iteroparity
  smtext <- strsplit(stanmodels[[stanmodel]]@model_code, "\\n")[[1]]
  
  pars <- stan_pars(stan_model = stan_model, 
                    pars = ifelse(pool_pops, "group", "hyper"), 
                    SR_fun = SR_fun)
  
  RRS_opts <- sapply(pars, function(.par) {
    regex <- paste0(.par, "_[HW];")  # finds H/W parameter declaration 
    return(ifelse(any(grepl(regex, smtext)), .par, NA))
  })
  
  RRS_opts <- c("none", unique((RRS_opts)))
  RRS_check <- RRS %in% RRS_opts
  
  if(!all(RRS_check))
    stop(paste(RRS[!RRS_check], collapse = ", "), 
         " is not a valid RRS specification in ", stan_model, ".\n",
         "  Available options are: ", paste(na.omit(RRS_opts), collapse = ", "), ".")
}


