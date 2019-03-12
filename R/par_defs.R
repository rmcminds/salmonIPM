#' Generates a table of model parameters, their dimensions, and their 
#' definitions.
#'
#' @param stan_model Character string giving the name of the Stan model being
#'   fit (".stan" filetype extension is not included).   
#'
par_defs <- function(stan_model) {
  hdrs <- c("Param", "Dims", "Defn")
  full_tbl <- list(
    c("alpha", "N_pops x 1", "Intrinsic productivity"),
    c("mu_alpha", "scalar", "Hyper-mean of log intrinsic productivity"),
    c("sigma_alpha", "scalar", "Hyper-SD of log intrinsic productivity"),
    c("Rmax", "N_pops x 1", "Asymptotic recruitment"),
    c("mu_Rmax", "scalar", "Hyper-mean of log asymptotic recruitment"),
    c("sigma_Rmax", "scalar", "Hyper-SD of log asymptotic recruitment"),
    c("beta", "N_pop x N_X", "Regression coefs for log productivity anomalies"),
    c("rho", "N_pop x 1", "AR(1) coefs for log productivity anomalies"),
    c("sigma", "N_pop x 1", "SD of process errors"),
    c("mu_p", "N_pop x 1", "Mean age distribution of the popn"),
    c("sigma_p", "N_pop x (N_age-1)", "SDs of log-ratio cohort age distribution"),
    c("R_p", "N_pop x N_pop", "Correlation matrix of within-popn cohort log-ratio age distns"),
    c("p", "(N_pop x N_year) x N_age", "Year-specific cohort age distributions"),
    c("p_HOS", "N_H x 1", "True prop of HOS in years in with H fish present"),
    c("B_rate_all", "(N_pop x N_year) x 1", "True broodstock take rate in all years"),
    c("tau", "(N_pop x N_year) x 1", "SD of observation errors of total spawners"),
    c("S", "(N_pop x N_year) x 1", "True total spawner abundance"),
    c("R", "(N_pop x N_year) x 1", "True recruit abundance (not density) by brood year"),
    c("q", "(N_pop x N_year) x N_age", "True spawner age distributions"),
    c("", "", ""),
    c("", "", ""),
  )
  idx <- stan_pars(stan_model)
}