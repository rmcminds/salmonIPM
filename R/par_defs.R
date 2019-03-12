#' Returns a list of the parameters returned by default for the specified
#' model.
#'
#' @param stan_model One of two options:
#' \enumerate{
#'   \item A \code{character} string giving the name of the Stan model being fit
#'   (without the ".stan" filetype extension)
#'   \item An object of class \code{stanfit} as returned by a call to
#'   \code{\link{salmonIPM}}:
#' }
#' 
#' The default argument \code{stan_model = NULL} returns a list of all of the 
#' parameters returned by default, across all model forms.
#'
#' @return Data frame with columns for parameter name, its dimensions, and its
#' definition
#' 
#' @export
par_defs <- function(stan_model = NULL) {
  par_list <- list(
    c("alpha", "N_pop x 1", "Intrinsic productivity"),
    c("mu_alpha", "scalar", "Hyper-mean of log intrinsic productivity"),
    c("sigma_alpha", "scalar", "Hyper-SD of log intrinsic productivity"),
    c("Rmax", "N_pop x 1", "Asymptotic recruitment"),
    c("mu_Rmax", "scalar", "Hyper-mean of log asymptotic recruitment"),
    c("sigma_Rmax", "scalar", "Hyper-SD of log asymptotic recruitment"),
    c("rho_alphaRmax", "scalar", "Correlation between log(alpha) and log(Rmax)"),
    c("beta", "N_pop x N_X", "Regression coefs for log productivity anomalies"),
    c("rho", "N_pop x 1", "AR(1) coefs for log productivity anomalies"),
    c("sigma", "N_pop x 1", "SD of process errors"),
    c("mu_p", "N_pop x 1", "Mean age distribution of the popn"),
    c("sigma_p", "N_pop x (N_age-1)", "SDs of log-ratio cohort age distribution"),
    c("R_p", "N_pop x N_pop", "Correlation matrix of within-popn cohort log-ratio age distns"),
    c("p", "N x N_age", "Year-specific cohort age distributions"),
    c("p_HOS", "N_H x 1", "True prop of hatchery-origin spawners (HOS) in years with H fish present"),
    c("B_rate_all", "N x 1", "True broodstock take rate in all years"),
    c("tau", "N x 1", "SD of observation errors of total spawners"),
    c("S", "N x 1", "True total spawner abundance"),
    c("R", "N x 1", "True recruit abundance (not density) by brood year"),
    c("q", "N x N_age", "True spawner age distributions"),
    c("beta_phi", "N_X x 1", "Regression coefs for log productivity anomalies"),
    c("sigma_phi", "scalar", "Hyper-SD of brood year log productivity anomalies"),
    c("rho_phi", "scalar", "AR(1) coef for log productivity anomalies"),
    c("sigma_gamma", "(N_age-1) x 1", "Among-pop SD of mean log-ratio age distributions"),
    c("R_gamma", "(N_age-1) x (N_age-1)", "Among-pop correlation matrix of mean log-ratio age distns"),
    c("beta_M", "N_pop x N_X_M", "Regression coefs for spawner-smolt productivity"),
    c("rho_M", "N_pop x 1", "AR(1) coefs for spawner-smolt productivity"),
    c("sigma_M", "N_pop x 1", "SD of spawner-smolt process error"),
    c("mu_MS", "N_pop x 1", "Mean smolt-to-adult return (SAR)"),
    c("beta_MS", "N_pop x N_X_MS", "Regression coefs for SAR"),
    c("rho_MS", "N_pop x 1", "AR(1) coefs for SAR"),
    c("sigma_MS", "N_pop x 1", "SD of SAR process error"),
    c("s_MS", "N x 1", "True SAR by outmigration year"),
    c("tau_M", "N_pop x 1", "SD of smolt observation error"),
    c("tau_S", "N_pop x 1", "SD of spawner observation error"),
    c("M", "N x 1", "True smolt abundance (not density) by outmigration year"),
    c("mu_p_M", "N_Mage x N_pop", "Popn mean smolt age distributions"),
    c("sigma_p_M", "N_pop x (N_Mage-1)", "SDs of log-ratio cohort smolt age distribution"),
    c("R_p_M", "[(N_Mage-1) x (N_Mage-1)] x N_pop", "Correlation matrices of log-ratio smolt age distns"),
    c("p_M", "N x N_Mage", "True smolt age distributions by brood year"),
    c("q_M", "N x N_Mage", "True smolt age distributions by calendar year"),
    c("R_MS", "(N_Mage x N_Mage) x N_pop", "Correlation matrices of logit SAR by smolt age"),
    c("mu_p_MS", "N_MSage x N_pop x N_Mage", "Popn mean ocean age distributions for each smolt age"),
    c("sigma_p_MS", "(N_MSage-1) x N_pop x N_Mage", "SDs of log-ratio ocean age for each smolt age"),
    c("R_p_MS", "[(N_Mage*(N_MSage-1)) x (N_Mage*(N_MSage-1))] x N_pop", "Correlation matrices of log-ratio ocean age distns"),
    c("p_MS", "N_MSage x N x N_Mage", "True ocean age distns by outmigration year"),
    c("q_MS", "N x N_MSage", "True ocean age distns of spawners"),
    c("q_GR", "N x N_MSage", "true Gilbert-Rich age distns of spawners"),
    c("R_hat", "N x 1", "Expected recruit abundance (not density) by brood year"),
    c("S_sim", "N x 1", "Simulated number of spawners"),
    c("R_sim", "N x 1", "Simulated number of recruits"),
    c("gamma", "N_pop x (N_age-1)", "Popn mean log-ratio age distributions"),
    c("phi", "N_year_all x 1", "log of productivity anomalies by brood year")
  )
  full_tbl <- do.call(rbind, par_list)
  colnames(full_tbl) <- c("Parameter/state", "Dimensions", "Definition")
  if(is.null(stan_model)) {
    full_tbl <- full_tbl
  } else {
    cls <- class(stan_model)
    if(cls == "character") {
      idx <- full_tbl[,"Parameter/state"] %in% stan_pars(stan_model)
    } else {
      if(cls == "stanfit") {
        idx <- full_tbl[,"Parameter/state"] %in% stan_pars(stan_model@model_name)
      } else {
        stop("stan_model must be a character string or stanfit object, but was ",
             capture.output(dput(cls)),
             call. = TRUE)
      }
    }
    full_tbl <- full_tbl[idx,]
  }
  return(print(full_tbl, quote = FALSE))
}