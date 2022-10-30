#' Helper function that returns the parameters and states that a given
#' model should save.
#'
#' @param stan_model Character string giving the name of the Stan model being
#'   fit (`".stan"` filetype extension is not included).   
#'
#' @return Character vector with names of parameters and states that Stan will
#' save.
#' 
#' @export
stan_pars <- function(stan_model) {
  pars <- list( 
    IPM_SS_np = c("alpha","beta_alpha","Rmax","beta_Rmax","beta_R","rho_R","sigma_R",
                  "mu_p","sigma_p","R_p","p",
                  "p_HOS","tau","S","R","q"),
    IPM_SSiter_np = c("alpha","beta_alpha","Rmax","beta_Rmax","beta_R","rho_R","sigma_R",
                      "mu_p","sigma_p","R_p","p",
                      "mu_SS","beta_SS","rho_SS","sigma_SS","s_SS",
                      "p_HOS","tau","S","R","q"),
    IPM_SS_pp = c("mu_alpha","beta_alpha","sigma_alpha","alpha",
                  "mu_Rmax","beta_Rmax","sigma_Rmax","Rmax","rho_alphaRmax",
                  "beta_R","sigma_year_R","rho_R","eta_year_R","sigma_R",
                  "mu_p","sigma_pop_p","R_pop_p","mu_pop_alr_p","sigma_p","R_p","p",
                  "p_HOS","tau","S","R","q"),
    IPM_SMS_np = c("alpha","beta_alpha","Mmax","beta_Mmax","beta_M","rho_M","sigma_M",
                   "mu_MS","beta_MS","rho_MS","sigma_MS","s_MS",
                   "mu_p","sigma_p","R_p","p",
                   "p_HOS","tau_M","tau_S","S","M","q"),
    IPM_SMS_pp = c("mu_alpha","beta_alpha","sigma_alpha","alpha",
                   "mu_Mmax","beta_Mmax","sigma_Mmax","Mmax","rho_alphaMmax",
                   "beta_M","rho_M","sigma_year_M","eta_year_M","sigma_M",
                   "mu_MS","beta_MS","rho_MS","sigma_year_MS","eta_year_MS","sigma_MS","s_MS",
                   "mu_p","sigma_pop_p","R_pop_p","mu_pop_alr_p","sigma_p","R_p","p",
                   "p_HOS","tau_M","tau_S","S","M","q"),
    IPM_SMaS_np = c("alpha","beta_alpha","Mmax","beta_Mmax",
                    "beta_M","rho_M","sigma_M","tau_M","M",
                    "mu_p_M","sigma_p_M","R_p_M","p_M","q_M",
                    "mu_MS","beta_MS","rho_MS","sigma_MS","R_MS","s_MS",
                    "mu_p_MS","sigma_p_MS","R_p_MS","p_MS","q_MS",
                    "q_GR","p_HOS","tau_S","S"),
    IPM_LCRchum_pp = c("mu_E","sigma_E","delta_NG","mu_psi","beta_psi","sigma_psi","psi",
                       "mu_Mmax","beta_Mmax","sigma_Mmax","Mmax",
                       "beta_M","rho_M","sigma_year_M","eta_year_M","sigma_M",
                       "mu_MS","beta_MS","rho_MS","sigma_year_MS","eta_year_MS","sigma_MS","s_MS",
                       "mu_p","sigma_pop_p","R_pop_p","mu_pop_alr_p","sigma_p","R_p","p",
                       "mu_F","sigma_pop_F","sigma_F","p_F","p_origin",
                       "mu_tau_M","sigma_tau_M","tau_M","mu_tau_S","sigma_tau_S","tau_S",
                       "M","S","q","q_F"),
    IPM_ICchinook_pp = c("mu_alpha","beta_alpha","sigma_alpha","alpha",
                         "mu_Mmax","beta_Mmax","sigma_Mmax","Mmax","rho_alphaMmax",
                         "beta_M","rho_M","sigma_M","M",
                         "mu_D","beta_D","rho_D","sigma_D","s_D",
                         "mu_SAR","beta_SAR","rho_SAR","sigma_SAR","SAR",
                         "mu_U","beta_U","rho_U","sigma_U","s_U",
                         "mu_p","sigma_pop_p","R_pop_p","mu_pop_alr_p","sigma_p","R_p","p",
                         "p_HOS","tau_S","S","q"),
    RR_SS_np = c("alpha","Rmax","rho_R","sigma_R","R_hat","S_sim","R_sim"),
    RR_SS_pp = c("mu_alpha","sigma_alpha","alpha",
                 "mu_Rmax","sigma_Rmax","Rmax","rho_alphaRmax",
                 "rho_R","sigma_year_R","eta_year_R","sigma_R",
                 "R_hat","S_sim","R_sim")
  )
  return(pars[[stan_model]])
}
