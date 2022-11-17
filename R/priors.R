#' Prior distributions
#' 
#' @name priors
#' 
#' @description The functions described on this page are used to specify priors
#' on selected (hyper)parameters in **salmonIPM** models. 
#'   
#'   The default priors used in the various models
#'   are intended to be weakly informative in that they provide moderate
#'   regularization and help stabilize computation. Priors on scaling parameters, e.g.
#'   `Rmax` or `mu_Rmax`, are automatically adjusted to be weakly informative but
#'   consistent with overall population density. For many applications these
#'   defaults will perform well, but if external information not included in 
#'   `fish_data` is available, it can be incorporated via user-specified priors on
#'   key parameters. See **Details** for a table of available prior options.
#' 
#' @param mean Prior mean for normal or generalized normal distribution.
#' @param sd Prior standard deviation for normal distribution.
#' @param scale Prior scale for generalized normal distribution. Equivalent to
#' `alpha` in [gnorm::gnorm], but renamed to avoid confusion with the spawner-recruit
#' intrinsic productivity parameter. 
#' @param shape Prior shape for generalized normal distribution. Equivalent to `beta`
#' in [gnorm::gnorm].
#' @param meanlog,sdlog Prior log-scale mean and standard deviation, respectively, for 
#' lognormal distribution; see [stats::Lognormal].
#' @param a,b Prior shape parameters for the Beta distribution. Equivalent to
#' `shape1` and `shape2`, respectively, in [stats::Beta].
#' @param concentration Vector of shape parameters for the Dirichlet distribution.
#' Equivalent to `alpha` in [gtools::dirichlet], but renamed to avoid confusion with 
#' the spawner-recruit intrinsic productivity parameter.
#'
#' @details The table below shows the parameters in each model that can be given
#' user-specified priors and the corresponding distributions. Note that
#' users can modify the prior parameters but not the distribution families; attempting
#' to do the latter will result in an error.
#' 
#' The generalized normal density with `shape >> 1` is useful as a "soft-uniform" prior
#' to regularize the posterior away from regions of parameter space that may cause
#' computational or sampling problems. In the case of spawner and smolt observation error
#' log-SDs, the default prior bounds them &#8819; 0.1.
#' 
#' For parameters that are modeled as functions of covariates (see [salmonIPM::salmonIPM]),
#' the specified prior applies when all predictors are at their sample means.
#' 
#' |                    |             |                |              |            |              |            | **Parameter** |            |           |             |          |            |            |
#' |:-------------------|:-----------:|:--------------:|:------------:|:----------:|:------------:|:----------:|:-------------:|:----------:|:---------:|:-----------:|:--------:|:----------:|:----------:|
#' | **Model**          | **`alpha`** | **`mu_alpha`** | **`mu_psi`** | **`Rmax`** | **`mu_Rmax`**| **`Mmax`** | **`mu_Mmax`** | **`mu_MS`**| **mu_p**  | **`mu_SS`** | **`tau`**| **`tau_S`**| **`tau_M`**|
#' | `IPM_SS_np`        | lognormal   |                |              | lognormal  |              |            |               |            | dirichlet |             | gnormal  |            |            |
#' | `IPM_SSiter_np`    | lognormal   |                |              | lognormal  |              |            |               |            | dirichlet | beta        | gnormal  |            |            |
#' | `IPM_SS_pp`        |             | normal         |              |            | normal       |            |               |            | dirichlet |             | gnormal  |            |            |
#' | `IPM_SSiter_pp`    |             | normal         |              |            | normal       |            |               |            | dirichlet | beta        | gnormal  |            |            |
#' | `IPM_SMS_np`       | lognormal   |                |              |            |              | lognormal  |               | beta       | dirichlet |             |          | gnormal    | gnormal    |
#' | `IPM_SMS_pp`       |             | normal         |              |            |              |            | normal        | beta       | dirichlet |             |          | gnormal    | gnormal    |
#' | `IPM_SMaS_np`      | lognormal   |                |              |            |              | lognormal  |               |            |           |             |          | gnormal    | gnormal    |
#' | `IPM_ICchinook_pp` |             |                |              |            |              |            | normal        |            | dirichlet |             |          | gnormal    | gnormal    |
#' | `IPM_LCRchum_pp`   |             |                | beta         |            |              |            | normal        | beta       | dirichlet |             |          | gnormal    | gnormal    |
#' 
#' @return A named list to be used internally by the **salmonIPM** model-fitting functions.
#'   
NULL

#' @rdname priors
#' @export
normal <- function(mean = 0, sd = 1) 
{
  stopifnot(sd > 0)
  list(dist = "normal", mean = mean, sd = sd)
}

#' @rdname priors
#' @export
gnormal <- function(mean = 0, scale = 1, shape = 1)
{
  stopifnot(scale > 0 && shape > 0)
  list(dist = "gnormal", mean = mean, scale = scale, shape = shape)
}

#' @rdname priors
#' @export
lognormal <- function(meanlog = 0, sdlog = 1) 
{
  stopifnot(sdlog > 0)
  list(dist = "lognormal", meanlog = meanlog, sdlog = sdlog)
}

#' @rdname priors
#' @export
beta <- function(a = 1, b = 1)
{
  stopifnot(a > 0 && b > 0)
  list(dist = "beta", a = a, b = b)
}

#' @rdname priors
#' @export
dirichlet <- function(concentration = 1) 
{
  stopifnot(all(concentration > 0))
  list(dist = "dirichlet", concentration = concentration)
}




