#' Helper function to construct covariate model matrices from a list of two-sided formulas.
#' 
#' @param par_models A list of two-sided formulas of the form 
#' `theta ~ t1 + ... + tK`, where `theta` is a parameter or state in a `salmonIPM` model
#' that accepts covariates and `t1 ... tK` are terms involving variables in
#' `fish_data`; see [salmonIPM()] for details.
#' @param center Logical indicating whether the terms in model matrices 
#' constructed from `fish_data` using the formulas in `par_models` should be centered.
#' It is usually recommended to use the default (`TRUE`) so the baseline parameter
#' estimate applies when predictors are at their sample means, but in some cases such
#' as factor predictors `center = FALSE` may be appropriate.
#' @param scale  Logical indicating whether the terms in model matrices 
#' constructed from `fish_data` using the formulas in `par_models` should be scaled to have 
#' column SDs of 1.
#' @param fish_data See [salmonIPM()]. In particular, the columns `...` may be used on
#' the right-hand side of formulas in `par_models`.
#' 
#' @return A named list of model matrices corresponding to the elements of `par_models`.
#' 
#' @seealso [salmonIPM()] for specifying covariate effects via `par_models`, 
#'   [stan_data()] for assembling input data including covariate model matrices
#' 
#' @export

par_model_matrix <- function(par_models, center = TRUE, scale = TRUE, fish_data)
{
  X <- lapply(par_models, function(f) {
    stopifnot(attr(terms(f), "response") == 1) # formulas must be 2-sided
    ff <- update(f, NULL ~ .)
    mf <- model.frame(ff, data = fish_data, na.action = na.fail) # no missing covariates
    mm <- model.matrix(ff, data = mf)
    mm <- mm[, setdiff(colnames(mm), "(Intercept)"), drop = FALSE] # remove intercept
    mm <- scale(mm, center = center, scale = scale)
    return(mm)
  })
  names(X) <- lapply(par_models, function(f) all.vars(f)[1]) # names are responses
  return(X)
}
