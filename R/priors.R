#' Prior distributions
#' 
#' @name priors
#' 
#' @description These functions are used to specify priors
#' on selected (hyper)parameters in **salmonIPM** models. 
#'   
#'   The default priors used in the various models
#'   are intended to be weakly informative, in that they provide moderate
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
#' `alpha` in [`gnorm::gnorm`], but renamed to avoid confusion with the spawner-recruit
#' intrinsic productivity parameter. 
#' @param shape Prior shape for generalized normal distribution. Equivalent to `beta`
#' in [`gnorm::gnorm`].
#' @param meanlog,sdlog Prior log-scale mean and standard deviation, respectively, for 
#' lognormal distribution; see [`stats::Lognormal`].
#' @param a,b Prior shape parameters for the Beta distribution. Equivalent to
#' `shape1` and `shape2`, respectively, in [`stats::Beta`].
#' @param concentration Vector of shape parameters for the Dirichlet distribution.
#' Equivalent to `alpha` in [`gtools::dirichlet`], but renamed to avoid confusion with 
#' the spawner-recruit intrinsic productivity parameter.
#'
#' @details The table below shows the parameters in each model that can be given
#' user-specified priors and the corresponding distributions. Note that
#' users can modify the prior parameters but not the distribution families; attempting
#' to do the latter will result in an error.
#' 
#' Priors for parameters that are bounded on the positive real line 
#' (e.g. `tau`, `tau_S` and `tau_M`) are automatically left-truncated at zero.
#' 
#' For parameters that are modeled as functions of covariates (see [salmonIPM::salmonIPM()]),
#' the specified prior applies when all predictors are at their sample means.
#' 
#' The generalized normal density with `shape >> 1` is useful as a platykurtic "soft-uniform" 
#' prior to regularize the posterior away from regions of parameter space that may cause
#' computational or sampling problems. In the case of spawner and smolt observation error
#' log-SDs, the default prior bounds them &#8819; 0.1.
#' 
#' |                    |                         |                         |                     |                        |                        |                        | **Parameter (PDF)**     |                    |                        |                    |                     |                                   |
#' |:-------------------|:-----------------------:|:-----------------------:|:-----------------------:|:------------------:|:----------------------:|:----------------------:|:-----------------------:|:------------------:|:----------------------:|:------------------:|:-------------------:|:---------------------------------:|
#' | **Model**          | `alpha` \cr `lognormal` | `mu_alpha` \cr `normal` | `mu_psi` \cr `beta` | `Rmax` \cr `lognormal` | `mu_Rmax` \cr `normal` | `Mmax` \cr `lognormal` | `mu_Mmax` \cr `normal`  | `mu_MS` \cr `beta` | `mu_p` \cr `dirichlet` | `mu_SS` \cr `beta` | `tau` \cr `gnormal` | `tau_S` \cr `tau_M` \cr `gnormal` |
#' | `IPM_SS_np`        | &#x2611;                | &#x2610;                | &#x2610;            | &#x2611;               | &#x2610;               | &#x2610;               | &#x2610;                | &#x2610;           | &#x2611;               | &#x2610;           | &#x2611;            | &#x2610;                          |
#' | `IPM_SSiter_np`    | &#x2611;                | &#x2610;                | &#x2610;            | &#x2611;               | &#x2610;               | &#x2610;               | &#x2610;                | &#x2610;           | &#x2611;               | &#x2611;           | &#x2611;            | &#x2610;                          |
#' | `IPM_SS_pp`        | &#x2610;                | &#x2611;                | &#x2610;            | &#x2610;               | &#x2611;               | &#x2610;               | &#x2610;                | &#x2610;           | &#x2611;               | &#x2610;           | &#x2611;            | &#x2610;                          |
#' | `IPM_SSiter_pp`    | &#x2610;                | &#x2611;                | &#x2610;            | &#x2610;               | &#x2611;               | &#x2610;               | &#x2610;                | &#x2610;           | &#x2611;               | &#x2611;           | &#x2611;            | &#x2610;                          |
#' | `IPM_SMS_np`       | &#x2611;                | &#x2610;                | &#x2610;            | &#x2610;               | &#x2610;               | &#x2611;               | &#x2610;                | &#x2611;           | &#x2611;               | &#x2610;           | &#x2610;            | &#x2611;                          |
#' | `IPM_SMS_pp`       | &#x2610;                | &#x2611;                | &#x2610;            | &#x2610;               | &#x2610;               | &#x2610;               | &#x2611;                | &#x2611;           | &#x2611;               | &#x2610;           | &#x2610;            | &#x2611;                          |
#' | `IPM_SMaS_np`      | &#x2611;                | &#x2610;                | &#x2610;            | &#x2610;               | &#x2610;               | &#x2611;               | &#x2610;                | &#x2610;           | &#x2610;               | &#x2610;           | &#x2610;            | &#x2611;                          |
#' | `IPM_ICchinook_pp` | &#x2610;                | &#x2611;                | &#x2610;            | &#x2610;               | &#x2610;               | &#x2610;               | &#x2611;                | &#x2610;           | &#x2611;               | &#x2610;           | &#x2610;            | &#x2611;                          |
#' | `IPM_LCRchum_pp`   | &#x2610;                | &#x2610;                | &#x2611;            | &#x2610;               | &#x2610;               | &#x2610;               | &#x2611;                | &#x2611;           | &#x2611;               | &#x2610;           | &#x2610;            | &#x2610;                          |
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




