#' Summarize the priors used in a `salmonIPM` model
#'
#' Print a summary of the priors on hyperparameters in a fitted [salmonIPMfit]
#' model object.
#'
#' @name prior_summary.salmonIPMfit
#' @aliases prior_summary
#'
#' @param object An object of class [salmonIPMfit] with prior information stored
#'   in `prior.info`.
#' @param digits Number of decimal places to print.
#' @param ... Currently ignored.
#'
#' @return A character vector containing the prior summary information.
#'
#' @details This function is called mainly for its side effect of printing the
#'   prior summary. The character vector to be printed is returned invisibly,
#'   but the same information can be found in more usable form in
#'   `object$prior.info`. Priors are identified as user-specifiable or
#'   hard-coded in the printed output.
#'
#'   Note that correlation matrices are given
#'   [`lkj_corr`](https://mc-stan.org/docs/functions-reference/correlation_matrix_distributions.html)
#'   priors, but the actual Stan implementation uses the more efficient lower
#'   Cholesky factor version of the LKJ distribution.
#'
#' @seealso [prior_summary()], [priors], [salmonIPM()], [salmonIPMfit]
#'
#' @importFrom rstantools prior_summary
#' @exportS3Method rstantools::prior_summary
#' @export prior_summary

prior_summary.salmonIPMfit <- function(object, digits = 2, ...) {
  prior.info <- object$prior.info
  pars <- names(prior.info)
  types <- sapply(prior.info, attr, "type")
  user_priors <- as.list(object$call)$priors
  user_pars <- pars %in% sapply(user_priors, function(f) all.vars(f)[1])
  user <- ifelse(types == "user", paste0("[", ifelse(user_pars, "x", " "), "]  "), "     ")
  names(user) <- pars
  pad <- sapply(max(nchar(pars)) - nchar(pars), function(n) paste(rep(" ", n), collapse = ""))
  names(pad) <- pars
  
  prsum <- sapply(pars, function(.par) {
    prinfo <- prior.info[[.par]]
    dist <- prinfo$dist
    args <- lapply(prinfo[setdiff(names(prinfo), "dist")], round, digits = digits)
    args <- paste(paste0(names(args), " = ", args), collapse = ", ")
    paste0(user[.par], pad[.par], .par, " ~ ", dist, "(", args, ")\n") 
  })
  prsum <- unlist(prsum)
  
  cat("Priors for model ", object$stan_model, "\n", 
      "[ ] user-specifiable (default) \n",
      "[x] set by user \n\n",
      prsum,
      sep = "")
  
  invisible(prsum)
}

#' Extract prior info from Stan-formatted data and model code
#'
#' Helper function to extract both user-specifiable and hard-coded prior data 
#' and assemble it in the appropriate form for the `prior.info` element of a
#' [salmonIPMfit] object.
#'
#' @param stan_data Named list of input data passed to [rstan::sampling()] for
#'   fitting an IPM, as returned by [stan_data()].
#' @param stanmodel An object of S4 class `stanmodel` as constructed by
#'   [stan_model()] or by **rstantools** during installation of
#'   [salmonIPM-package] and stored internally in `salmonIPM:::stanmodels`.
#' @param pars Character vector of hyperparameters for which priors are
#'   specified.
#'
#' @return A named list of priors on `pars`. Each element is itself a named list
#'   representing the prior in the format returned by the functions in [priors]
#'   with an additional attribute `"type"`.
#'
#' @details Priors include (1) Those that are user-specifiable through the
#'   `salmonIPM(priors)` argument (whether modified or defaults) and returned in
#'   Stan format by [stan_data()], (2) those that are explicitly hard-coded in the
#'   `model` block of `stanmodel`, and (3) those that are implicitly defined by bounds
#'   on `parameter` declarations.
#'
#'   To extract type (1), `get_prior_info()` calls the [priors] functions using the
#'   contents of `stan_data`.
#'
#'   To extract type (2), `get_prior_info()` parses `stanmodel` into lines of text
#'   and then uses regex to pull out the sampling statement for each parameter
#'   as a string, which corresponds directly to a [priors] function call. 
#'   
#' @importFrom utils relist

get_prior_info <- function(stan_data, stanmodel, pars) 
{
  pars <- pars[!grepl("delta_.*alpha|delta_.*max", pars)] # H/W discounts don't have priors
  
  # User-specifiable priors from stan_data
  prior_data <- stan_data[paste0("prior_", pars)]
  prior_data <- prior_data[!sapply(prior_data, is.null)]
  names(prior_data) <- gsub("prior_", "", names(prior_data))
  user_priors <- lapply(prior_data, function(x) {
    structure(do.call(attr(x, "dist"), relist(x)), type = "user")
  })
  
  # Hard-coded priors from stanmodel
  #
  # .*      any character, 0 or more (e.g. indent spaces, any transformation function)
  # .par    hyperparameter name
  # .*      any character, 0 or more (e.g. closing paren of transformation, space) 
  # ~ *     this line contains a sampling statement (may be 0 or more spaces after "~")
  # ([^;]+) pdf capture group: any character except ";", 1 or more
  # ;       end of sampling statement
  # .*      any character, 0 or more (e.g. comment)
  hard_pars <- setdiff(pars, names(user_priors)) 
  code_pars <- gsub("R_", "L_", hard_pars)   # correlation matrices -> Cholesky factors
  smtext <- strsplit(stanmodel@model_code, "\\n")[[1]]
  hard_priors <- lapply(code_pars, function(.par) {
    regex <- paste0(".*", .par, ".*~ *([^;]+);.*")
    prtext <- gsub(regex, "\\1", grep(regex, smtext, value = TRUE))
    prtext <- gsub("_cholesky", "", prtext)  # Cholesky factors -> correlation matrices
    if(length(prtext)) structure(as.list(eval(parse(text = prtext))), type = "hard")
  })
  names(hard_priors) <- hard_pars
  hard_priors <- hard_priors[sapply(hard_priors, length) > 0]

  # Implicit priors declared by bounds in stanmodel
  #
  # .*          any character, 0 or more (e.g. indent spaces, variable type)
  # <lower *= * beginning of bounds (may or may not be spaces around "=")
  # (.+)        lb capture group: 1 or more characters 
  # , *         comma separating bounds followed by 0 or more spaces
  # upper *= *  upper bound declaration
  # (.+)        ub capture group 
  # >.*         end of bounds followed by 0 or more characters (e.g. dims, spaces)
  # .par        hyperparameter name
  # ;.*         end of declaration followed by 0 or more characters (e.g. comment)
  bounded_pars <- setdiff(pars, c(names(user_priors), names(hard_priors)))
  bounded_priors <- lapply(bounded_pars, function(.par) {
    if(grepl("rho_alpha.max", .par)) {  # hack b/c declared parameter is L_alpha[R/M]max
      structure(uniform(-1,1), type = "bounded")
    } else {
      regex <- paste0(".*<lower *= *(.+), *upper *= *(.+)>.*", .par, ";.*") 
      bounds <- gsub(regex, "\\1 \\2", grep(regex, smtext, value = TRUE))
      bounds <- as.numeric(unlist(strsplit(bounds, " ")))
      if(length(bounds)) structure(do.call(uniform, as.list(bounds)), type = "bounded")
    }
  })
  names(bounded_priors) <- bounded_pars
  bounded_priors <- bounded_priors[sapply(bounded_priors, length) > 0]
  
  all_priors <- c(user_priors, hard_priors, bounded_priors)
  return(all_priors[match(pars, names(all_priors))])
}
