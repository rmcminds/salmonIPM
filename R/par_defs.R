#' Definitions of parameters and states in a **salmonIPM** model
#'
#' Returns a data frame of parameters with their types, dimensions and
#' definitions in the specified model.
#'
#' @param stan_model Character string giving the name of the model. See
#'   [salmonIPM()] for details.
#' @param pars An optional character vector specifying one or more hierarchical
#'   levels of parameters. Options are `"all"` (the default), `"hyper"`
#'   (top-level hyperparameters that are given priors), `"group"` (`pop`- or
#'   `year`-level parameters shared by multiple states), `"states"` (the lowest
#'   level, corresponding to unique rows in `fish_data`), and `"ppd"` (only if
#'   `model == "RR"`, observation-level predictions drawn from the posterior
#'   predictive distribution).
#' @param object A [salmonIPMfit] object. If this is provided then `SR_fun`,
#'   `RRS` and `par_models` are not needed and will be ignored; their values are
#'   extracted from `object`.
#' @inheritParams salmonIPM
#'
#' @return Data frame with columns listing the parameters and/or states and
#'   their definitions, types and dimensions.
#'
#' @importFrom dplyr tibble
#'
#' @export
par_defs <- function(stan_model = c("IPM_SS_np","IPM_SSiter_np","IPM_SS_pp","IPM_SSiter_pp",
                                    "IPM_SMS_np","IPM_SMS_pp","IPM_SMaS_np",
                                    "IPM_LCRchum_pp","RR_SS_np","RR_SS_pp"), 
                     pars = c("all","hyper","group","states","ppd"), 
                     SR_fun = "BH", RRS = "none", par_models = NULL, object = NULL) 
{
  if(!is.null(object)) {
    stopifnot("salmonIPMfit" %in% class(object))
    stan_model <- object$stan_model
  } else {
    stan_model <- match.arg(stan_model)
  }
  pars <- stan_pars(stan_model = stan_model, pars = pars, SR_fun = SR_fun, RRS = RRS,
                    par_models = par_models, object = object)
  stanmodel <- gsub("iter", "", stan_model)  # same Stan code for iteroparity
  
  # Parse model code and extract parameter declarations
  smtext <- strsplit(stanmodels[[stanmodel]]@model_code, "\\n")[[1]]
  pd <- data.frame(par = pars, def = NA, type = NA)

  for(.par in pars) {
    # .+           leading space(s) and type declaration, then space before .par
    # (?: = .+)?   optional ")?" non-capturing group "(?:" for declaration-assignment
    # ; +.*        end of statement, one or more spaces, zero or more other chars e.g. "//?"
    # // (.+)      capture group for comment 
    dregex <- paste0(".+ ", .par, "(?: = .+)?; +.*// (.+)")
    pd$def[pd$par == .par] <- gsub(dregex, "\\1", grep(dregex, smtext, value = TRUE)[1])
    #  *    leading space(s) 
    # (.+)  capture group for type declaration, then space before .par
    tregex <- paste0(" *(.+) ", .par, "(?: = .+)?;.+")
    pd$type[pd$par == .par] <- gsub(tregex, "\\1", grep(tregex, smtext, value = TRUE)[1])
  }
  
  return(tibble(pd))
}
