#' Spawner-recruit functions
#' 
#' Compute recruitment given spawner abundance, a spawner-recruit function and parameters.
#' 
#' @param RRS `r lifecycle::badge("experimental")` 
#'   A character string or vector of strings naming parameters of the function specified by 
#'   `SR_fun` (i.e., `"alpha"` and/or `"Rmax"`) that differ between wild- and hatchery-origin 
#'   spawners, such that the relative reproductive success of hatchery spawners is not equal
#'   to 1. The default is `"none"`.
#' @param alpha If `RRS == "none"`, numeric vector, matrix, data frame, [posterior::rvar()],
#'   or `draws` of intrinsic productivity (i.e., recruits per spawner at zero spawner 
#'   density; slope of the spawner-recruit function at the origin).
#' @param alpha_W,alpha_H If `"alpha" %in% RRS`, numeric vectors, matrices, data frames, 
#'   [posterior::rvar()] or `draws` objects of intrinsic productivity for wild- and hatchery-
#'   origin spawners, respectively.
#' @param Rmax If `RRS == "none"`, numeric vector, matrix, data frame, [posterior::rvar()], 
#'   or `draws` of maximum recruitment per unit of habitat (length or area). This corresponds to 
#'   the asymptote of the Beverton-Holt or the mode of the Ricker.
#' @param Rmax_W,Rmax_H If `"Rmax" %in% RRS`, numeric vectors, matrices, data frames, 
#'   [posterior::rvar()] or `draws` objects of maximum recruitment per unit habitat 
#'   for wild- and hatchery- origin spawners, respectively.
#' @param S If `RRS == "none"`, numeric vector, matrix, data frame, [posterior::rvar()], 
#'   or `draws` of spawner abundance.
#' @param S_W,S_H If `RRS != "none"`, numeric vectors, matrices, data frames, 
#'   [posterior::rvar()] or `draws` objects of wild- and hatchery-origin spawner abundance,
#'   respectively.
#' @param A A numeric vector, matrix, data frame, [posterior::rvar()], or `draws` of 
#'   spawning habitat size (either stream length or area), used to standardize `Rmax`. 
#'   The default is 1, in which case `Rmax` is in units of abundance (which is also density).
#' @param R_per_S Logical indicating whether to return recruits per spawner rather than
#'   recruits (the default).
#' @inheritParams salmonIPM
#'   
#' @return A vector, matrix, data frame, [posterior::rvar()] or `draws`, depending on the 
#'   argument types, containing either recruits or recruits per spawner.
#'   
#' @details The **salmonIPM** package uses a nonstandard parameterization of the Ricker
#'   model by the maximum recruitment `Rmax`. This is typically better identified by
#'   data than the carrying capacity or per capita density dependence, and it 
#'   facilitates a common interpretation and priors with the Beverton-Holt.
#'   
#'   Note that the functions for the `RRS != "none"` case are written in their most general
#'   form, with both `alpha` and `Rmax` differing between wild and hatchery spawners. If
#'   only one parameter is specified in `RRS`, then the `_W` and `_H` values of the other
#'   parameter are equal and the expression can be simplified.
#'   
#'   `RRS == "none"`
#'   |                       |                                                                       |
#'   |----------------------:|:----------------------------------------------------------------------|
#'   | Discrete exponential: | `alpha*S`                                                             |
#'   | Beverton-Holt:        | `alpha*S / (1 + alpha*S/(A*Rmax))`                                    |
#'   | Ricker:               | `alpha*S*exp(-alpha*S / (exp(1)*A*Rmax))`                             |
#'
#'   `RRS != "none"`
#'   |                       |                                                                                                     |
#'   |----------------------:|:----------------------------------------------------------------------------------------------------|
#'   | Discrete exponential: | `alpha_W*S_W + alpha_H*S_H`                                                                         |
#'   | Beverton-Holt:        | `(alpha_W*S_W + alpha_H*S_H) / (1 + alpha_W*S_W/(A*Rmax_W) + alpha_H*S_H/(A*Rmax_H))`               |
#'   | Ricker:               | `(alpha_W*S_W + alpha_H*S_H) * exp(-alpha_W*S_W/(exp(1)*A*Rmax_W) - alpha_H*S_H/(exp(1)*A*Rmax_H))` |
#'   
#'   Calculations are vectorized and elements of shorter arguments are recycled 
#'   as necessary.
#'   
#' @examples 
#' alpha <- 3
#' Rmax <- 1000
#' S <- 500
#' 
#' # default is Beverton-Holt
#' SR(alpha = alpha, Rmax = Rmax, S = S) 
#'
#' # approximately Rmax
#' SR(alpha = alpha, Rmax = Rmax, S = 1e6) 
#'
#' # scale Rmax by habitat area
#' SR(alpha = alpha, Rmax = Rmax, S = S, A = 0.1) 
#'
#' # discrete exponential ignores Rmax
#' SR(SR_fun = "exp", alpha = alpha, Rmax = Rmax, S = S) 
#'
#' # vectorization with recycling
#' SR(alpha = rep(alpha, 10), Rmax = rep(Rmax, 10), S = matrix(S, 10, 4)) 
#'
#' # return recruits per spawner
#' SR(alpha = alpha, Rmax = Rmax, S = S, R_per_S = TRUE) 
#'
#' # plot Ricker and show Rmax
#' curve(SR(SR_fun = "Ricker", alpha = alpha, Rmax = Rmax, S = x), from = 0, to = 2000,
#'       xlab = "Spawners", ylab = "Recruits", main = "Ricker")
#' abline(h = Rmax, lty = 2)
#'
#' # differential hatchery / wild relative reproductive success
#' alpha_W <- 3
#' alpha_H <- 2
#' Rmax_W <- 1000
#' Rmax_H <- 800
#' S_W <- 400
#' S_H <- 100
#' # compare to BH result above
#' SR(RRS = c("alpha","Rmax"), alpha_W = alpha_W, alpha_H = alpha_H, 
#'    Rmax_W = Rmax_W, Rmax_H = Rmax_H, S_W = S_W, S_H = S_H)
#'
#' @seealso [salmonIPM()] for fitting models, [simIPM()] for simulating data
#' 
#' @importFrom posterior is_rvar
#' @export

SR <- function(SR_fun = c("BH","B-H","bh","b-h","Ricker","ricker","exp"), RRS = "none", 
               alpha = NULL, alpha_W = NULL, alpha_H = NULL, 
               Rmax = NULL, Rmax_W = NULL, Rmax_H = NULL,
               S = NULL, S_W = NULL, S_H = NULL, A = 1, R_per_S = FALSE) 
{
  SR_fun <- match.arg(SR_fun)
  if(SR_fun == "DI") SR_fun <- "exp"
  if(SR_fun %in% c("B-H","bh","b-h")) SR_fun <- "BH"
  if(SR_fun == "ricker") SR_fun <- "Ricker"
  stopifnot(all(RRS %in% c("none","alpha","Rmax")))

  # Use most general form of unequal-RRS model to avoid repetition
  if(identical(RRS, "alpha")) {
    Rmax_W <- Rmax
    Rmax_H <- Rmax
  }
  if(identical(RRS, "Rmax")) {
    alpha_W <- alpha
    alpha_H <- alpha
  }
  
  if(is.data.frame(alpha)) alpha <- as.matrix(alpha)
  if(is.data.frame(alpha_W)) alpha_W <- as.matrix(alpha_W)
  if(is.data.frame(alpha_H)) alpha_H <- as.matrix(alpha_H)
  if(is.data.frame(Rmax)) Rmax <- as.matrix(Rmax)
  if(is.data.frame(Rmax_W)) Rmax_W <- as.matrix(Rmax_W)
  if(is.data.frame(Rmax_H)) Rmax_H <- as.matrix(Rmax_H)
  if(is.data.frame(S)) S <- as.matrix(S)
  if(is.data.frame(S_W)) S_W <- as.matrix(S_W)
  if(is.data.frame(S_H)) S_H <- as.matrix(S_H)
  if(is.data.frame(A)) A <- as.matrix(A)
  
  if(identical(RRS, "none")) { # relative reproductive success == 1
    R <- switch(
      SR_fun,
      exp = alpha*S,
      BH = alpha*S / (1 + alpha*S/(A*Rmax)),
      Ricker = alpha*S*exp(-alpha*S/(exp(1)*A*Rmax))
    ) 
  } else {  # relative reproductive success != 1
    S <- S_W + S_H
    R <- switch(
      SR_fun,
      exp = alpha_W*S_W + alpha_H*S_H,
      BH = (alpha_W*S_W + alpha_H*S_H) / (1 + alpha_W*S_W/(A*Rmax_W) + alpha_H*S_H/(A*Rmax_H)),
      Ricker = (alpha_W*S_W + alpha_H*S_H) * 
        exp(-alpha_W*S_W/(exp(1)*A*Rmax_W) - alpha_H*S_H/(exp(1)*A*Rmax_H)),
    )
  }
  
  return(if(R_per_S) R/S else R) 
}

