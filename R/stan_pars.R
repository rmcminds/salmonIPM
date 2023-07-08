#' Return the (hyper)parameters and states in a specified **salmonIPM** model
#'
#' @param stan_model Character string giving the name of the model.
#'  See [salmonIPM()] for details.
#' @param pars An optional character vector specifying one or more hierarchical levels 
#' of parameters. Options are `"all"` (the default), `"hyper"` (top-level parameters that
#' are given priors), `"group"` (`pop`- or `year`-level parameters shared by multiple states),
#' `"states"` (the lowest level, corresponding to unique rows in `fish_data`), and
#' `"ppd"` (only if `model == "RR"`, observation-level predictions drawn from the posterior
#' predictive distribution).
#' 
#' Using pars="beta" will restrict the displayed parameters to only the regression coefficients (without the intercept). "alpha" can also be used as a shortcut for "(Intercept)". If the model has varying intercepts and/or slopes they can be selected using pars = "varying".
#'
#' @return Character vector with names of selected parameters and states
#' 
#' @export

stan_pars <- function(stan_model = c("IPM_SS_np","IPM_SSiter_np","IPM_SS_pp","IPM_SSiter_pp",
                                     "IPM_SMS_np","IPM_SMS_pp","IPM_SMaS_np",
                                     "IPM_LCRchum_pp","IPM_ICchinook_pp",
                                     "RR_SS_np","RR_SS_pp"), 
                      pars = c("all","hyper","group","states","ppd")) 
{
  stan_model <- match.arg(stan_model)
  pars <- match.arg(pars, several.ok = TRUE)
  
  par_list <- list( 
    IPM_SS_np = list(
      hyper = c("alpha","beta_alpha","Rmax","beta_Rmax","beta_R","rho_R","sigma_R",
                "mu_p","sigma_p","R_p","tau"),
      states = c("R","p","S","q","p_HOS")
    ),
    
    IPM_SSiter_np = list(
      hyper = c("alpha","beta_alpha","Rmax","beta_Rmax","beta_R","rho_R","sigma_R",
                "mu_p","sigma_p","R_p","mu_SS","beta_SS","rho_SS","sigma_SS","tau"),
      states = c("R","p","s_SS","S","q","p_HOS")
    ),
    
    IPM_SS_pp = list(
      hyper = c("mu_alpha","beta_alpha","sigma_alpha","mu_Rmax","beta_Rmax","sigma_Rmax",
                "rho_alphaRmax","beta_R","sigma_year_R","rho_R","sigma_R",
                "mu_p","sigma_pop_p","R_pop_p","sigma_p","R_p","tau"),
      group = c("alpha","Rmax","eta_year_R","mu_pop_alr_p"),
      states = c("R","p","S","q","p_HOS")
    ),
    
    IPM_SSiter_pp = list(
      hyper = c("mu_alpha","beta_alpha","sigma_alpha","mu_Rmax","beta_Rmax","sigma_Rmax",
                "rho_alphaRmax","beta_R","sigma_year_R","rho_R","sigma_R",
                "mu_p","sigma_pop_p","R_pop_p","sigma_p","R_p",
                "mu_SS","beta_SS","rho_SS","sigma_year_SS","sigma_SS","tau"),
      group = c("alpha","Rmax","eta_year_R","mu_pop_alr_p","eta_year_SS"),
      states = c("R","p","s_SS","S","p_HOS","q")
    ),
    
    IPM_SMS_np = list(
      hyper = c("alpha","beta_alpha","Mmax","beta_Mmax","beta_M","rho_M","sigma_M",
                "mu_MS","beta_MS","rho_MS","sigma_MS", "mu_p","sigma_p","R_p","tau_M","tau_S"),
      states = c("M","s_MS","p","S","q","p_HOS")
    ),
    
    IPM_SMS_pp = list(
      hyper = c("mu_alpha","beta_alpha","sigma_alpha","mu_Mmax","beta_Mmax","sigma_Mmax",
                "rho_alphaMmax","beta_M","rho_M","sigma_year_M","sigma_M","tau_M",
                "mu_MS","beta_MS","rho_MS","sigma_year_MS","sigma_MS",
                "mu_p","sigma_pop_p","R_pop_p","sigma_p","R_p","tau_S"),
      group = c("alpha","Mmax","eta_year_M","eta_year_MS","mu_pop_alr_p"),
      states = c("M","s_MS","p","S","q","p_HOS")
    ),
    
    IPM_SMaS_np = list(
      hyper = c("alpha","beta_alpha","Mmax","beta_Mmax","beta_M","rho_M","sigma_M",
                "mu_p_M","sigma_p_M","R_p_M","tau_M",
                "mu_MS","beta_MS","rho_MS","sigma_MS","R_MS",
                "mu_p_MS","sigma_p_MS","R_p_MS","tau_S"),
      states = c("p_M","M","q_M","s_MS","p_MS","S","q_MS","q_GR","p_HOS")
    ),
    
    IPM_ICchinook_pp = list(
      hyper = c("mu_alpha","beta_alpha","sigma_alpha","mu_Mmax","beta_Mmax","sigma_Mmax",
                "rho_alphaMmax","beta_M","rho_M","sigma_M",
                "mu_D","beta_D","rho_D","sigma_D",
                "mu_SAR","beta_SAR","rho_SAR","sigma_SAR",
                "mu_U","beta_U","rho_U","sigma_U",
                "mu_p","sigma_pop_p","R_pop_p","sigma_p","R_p","tau_S"),
      group = c("alpha","Mmax","s_D","SAR","s_U","mu_pop_alr_p"),
      states = c("M","p","p_HOS","S","q")
    ),
    
    IPM_LCRchum_pp = list(
      hyper = c("mu_E","sigma_E","delta_NG","mu_psi","beta_psi","sigma_psi",
                "mu_Mmax","beta_Mmax","sigma_Mmax",
                "beta_M","rho_M","sigma_year_M","sigma_M",
                "mu_MS","beta_MS","rho_MS","sigma_year_MS","sigma_MS",
                "mu_p","sigma_pop_p","R_pop_p","sigma_p","R_p",
                "mu_F","sigma_pop_F","sigma_F","P_D",
                "mu_tau_M","sigma_tau_M","mu_tau_S","sigma_tau_S"),
      group = c("psi","Mmax","eta_year_M","eta_year_MS","mu_pop_alr_p"),
      states = c("M","tau_M","s_MS","p","p_F","S","tau_S","q","q_F","q_O","p_HOS")
    ),
    
    RR_SS_np = list(
      hyper = c("alpha","Rmax","rho_R","sigma_R","R_hat"),
      ppd = c("S_sim","R_sim")
    ),
    
    RR_SS_pp = list(
      hyper = c("mu_alpha","sigma_alpha","mu_Rmax","sigma_Rmax","rho_alphaRmax",
                "rho_R","sigma_year_R","sigma_R"),
      group = c("alpha","Rmax","eta_year_R"),
      ppd = c("R_hat","S_sim","R_sim")
    )
  )
  
  if(identical(pars, "all")) {
    return(unlist(par_list[[stan_model]]))
  } else {
    return(unlist(par_list[[stan_model]][pars]))
  }
}
