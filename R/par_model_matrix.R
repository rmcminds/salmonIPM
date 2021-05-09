#' Helper function to construct covariate model matrices from a list of two-sided formulas.
#' 
#' @param par_models A list of two-sided formulas of the form 
#' `theta ~ t1 + ... + tK`, where `theta` is a parameter or state in a `salmonIPM` model
#' that accepts covariates and `t1 ... tK` are terms involving variables in
#' `fish_data`; see [salmonIPM()] for details.
#' @param scale  Logical indicating whether the main effects in model matrices 
#' constructed from `fish_data` using the formulas in `par_models` should be scaled to have 
#' column SDs of 1 in addition to being centered (`TRUE`) or centered only (`FALSE`).
#' Interactions are computed from centered and (possibly) scaled predictors.
#' @param fish_data See [salmonIPM()]. In particular, the columns `...` may be used on
#' the right-hand side of formulas in `par_models`.
#' 
#' @return A named list of model matrices corresponding to the elements of `par_models`.
#' 
#' @export

par_model_matrix <- function(par_models, scale = TRUE, fish_data)
{
  X <- lapply(par_models, function(f) {
    stopifnot(attr(terms(f), "response") == 1) # formulas must be 2-sided
    ff <- update(f, NULL ~ . + 1) # force placeholder intercept
    mf <- model.frame(ff, data = fish_data, na.action = na.fail) # no missing covariates
    mf <- as.data.frame(scale(mf, scale = scale))
    mm <- model.matrix(ff, data = mf) 
    mm <- subset(mm, select = -`(Intercept)`, drop = FALSE)
    return(mm)
  })
  names(X) <- lapply(par_models, function(f) all.vars(f)[1]) # names are responses
  return(X)
}
