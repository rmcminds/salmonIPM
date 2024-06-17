#' The 'salmonIPM' package.
#'
#' @name salmonIPM-package
#' @aliases salmonIPM-package
#' @useDynLib salmonIPM, .registration = TRUE
#' @import methods
#' @import stats
#' @import Rcpp
#' @import gnorm
#' @importFrom utils capture.output head tail
#' @importFrom RcppParallel RcppParallelLibs CxxFlags
#' @importFrom rstan sampling
#' @importFrom mvtnorm rmvnorm
#' @keywords internal
#'
#' @references
#' Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.19.3. https://mc-stan.org
#'
"_PACKAGE"
