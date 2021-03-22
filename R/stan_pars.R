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
    IPM_SS_np = c("alpha","Rmax","beta","rho","sigma",
                  "mu_p","sigma_p","R_p","p",
                  "p_HOS","B_rate_all","tau","S","R","q"),
    IPM_SS_pp = c("mu_alpha","sigma_alpha","alpha",
                  "mu_Rmax","sigma_Rmax","Rmax","rho_alphaRmax",
                  "beta_phi","sigma_phi","rho_phi","phi",
                  "mu_p","sigma_gamma","R_gamma","gamma",
                  "sigma_p","R_p","p",
                  "p_HOS","B_rate_all",
                  "sigma","tau","S","R","q"),
    IPM_SSpa_pp = c("mu_alpha","sigma_alpha","alpha",
                  "mu_Rmax","sigma_Rmax","Rmax","rho_alphaRmax",
                  "beta_phi","sigma_phi","rho_phi","phi",
                  "mu_p","sigma_gamma","R_gamma","gamma",
                  "sigma_p","R_p","p",
                  "p_HOS","B_rate_all",
                  "sigma","tau","S","R","q"),
    IPM_SMS_np = c("alpha","Rmax","beta_M","rho_M","sigma_M",
                   "mu_MS","beta_MS","rho_MS","sigma_MS","s_MS",
                   "mu_p","sigma_p","R_p","p","p_HOS","B_rate_all",
                   "tau_M","tau_S","S","M","q"),
    IPM_SMS_pp = c("mu_alpha","sigma_alpha","alpha",
                   "mu_Rmax","sigma_Rmax","Rmax","rho_alphaRmax",
                   "beta_phi_M","rho_phi_M","sigma_phi_M","phi_M","sigma_M",
                   "mu_MS","beta_phi_MS","rho_phi_MS","sigma_phi_MS","phi_MS",
                   "sigma_MS","s_MS",
                   "mu_p","sigma_gamma","R_gamma","gamma",
                   "sigma_p","R_p","p",
                   "p_HOS","B_rate_all",
                   "tau_M","tau_S","S","M","q"),
    IPM_SMaS_np = c("alpha","Rmax","beta_M","rho_M","sigma_M","tau_M","M",
                    "mu_p_M","sigma_p_M","R_p_M","p_M","q_M",
                    "mu_MS","beta_MS","rho_MS","sigma_MS","R_MS","s_MS",
                    "mu_p_MS","sigma_p_MS","R_p_MS","p_MS","q_MS",
                    "q_GR","p_HOS","B_rate_all","tau_S","S"),
    IPM_LCRchum_pp = c("mu_E","sigma_E","mu_psi","sigma_psi","psi",
                       "mu_Mmax","sigma_Mmax","Mmax","rho_psiMmax",
                       "beta_M","rho_M","sigma_year_M","eta_year_M","sigma_M",
                       "mu_MS","beta_MS","rho_MS","sigma_year_MS","eta_year_MS",
                       "sigma_MS","s_MS",
                       "mu_p","sigma_pop_p","R_pop_p","eta_pop_p","sigma_p","R_p","p",
                       "p_HOS","B_rate_all",
                       "mu_tau_M","sigma_tau_M","tau_M","mu_tau_S","sigma_tau_S","tau_S",
                       "M","S","q"),
    IPM_ICchinook_pp =  c("mu_alpha","sigma_alpha","alpha",
                          "mu_Rmax","sigma_Rmax","Rmax","rho_alphaRmax",
                          "beta_M","rho_M","sigma_M","M",
                          "mu_D","beta_D","rho_D","sigma_D","s_D",
                          "mu_SAR","beta_SAR","rho_SAR","sigma_SAR","SAR",
                          "mu_U","beta_U","rho_U","sigma_U","s_U",
                          "mu_p","sigma_gamma","R_gamma","gamma",
                          "sigma_p","R_p","p",
                          "p_HOS","B_rate_all",
                          "tau_S","S","q"),
    RR_SS_np = c("alpha","Rmax","rho","sigma","R_hat","S_sim","R_sim"),
    RR_SS_pp = c("mu_alpha","sigma_alpha","alpha",
                 "mu_Rmax","sigma_Rmax","Rmax","rho_alphaRmax",
                 "rho_phi","sigma_phi","phi","sigma",
                 "R_hat","S_sim","R_sim")
  )
  return(pars[[stan_model]])
}
