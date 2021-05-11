#' Spawner-recruit functions
#' 
#' Compute recruitment given spawner abundance, a spawner-recruit function and parameters.
#' 
#' @param SR_fun One of `"exp"` (density-independent discrete exponential), 
#'   `"BH"` (Beverton-Holt, the default), or `"Ricker"`. 
#'   Synonyms `"B-H"`, `"bh"`, `"b-h"` and `"ricker"` are also accepted.
#' @param alpha A numeric vector, matrix or data frame of intrinsic productivity 
#'   (i.e., recruits per spawner at zero spawner density; slope of the spawner-recruit 
#'   function at the origin.).
#' @param Rmax A numeric vector, matrix or data frame of maximum recruitment per unit 
#'   of habitat (length or area). This corresponds to the asymptote of the Beverton-Holt 
#'   or the mode of the Ricker.
#' @param S A numeric vector, matrix or data frame of spawner abundance.
#' @param A A numeric vector, matrix or data frame of spawning habitat size 
#'   (either stream length or area), used to standardize `Rmax`. The default is 1, 
#'   in which case `Rmax` is in units of abundance rather than density.
#' @param R_per_S Logical indicating whether to return recruits per spawner rather than
#'   recruits (the default). 
#'   
#' @details The `salmonIPM` package uses a nonstandard parameterization of the Ricker
#'   model by the maximum recruitment `Rmax`. This is typically better identified by
#'   data than the carrying capacity or per capita density dependence, and it 
#'   facilitates a common interpretation and priors with the Beverton-Holt.
#'   
#'   |                       |                                                         |
#'   |----------------------:|:--------------------------------------------------------|
#'   | Discrete exponential: | `R = alpha * S`                                         |
#'   | Beverton-Holt:        | `R = alpha * S / (1 + alpha * S / (A * Rmax))`          |
#'   | Ricker:               | `R = alpha * S * exp(-alpha * S / (A * exp(1) * Rmax))` |
#'   
#'   Calculations are vectorized and elements of shorter arguments are recycled 
#'   as necessary.
#'   
#' @examples 
#' alpha <- 3
#' Rmax <- 1000
#' S <- 500
#' 
#' SR(alpha = alpha, Rmax = Rmax, S = S) # default is Beverton-Holt
#' SR(alpha = alpha, Rmax = Rmax, S = 1e6) # approximately Rmax
#' SR(alpha = alpha, Rmax = Rmax, S = S, A = 0.1) # scale Rmax by habitat area
#' SR(SR_fun = "exp", alpha = alpha, Rmax = Rmax, S = S) # discrete exponential ignores Rmax
#' SR(alpha = rep(alpha, 10), Rmax = rep(Rmax, 10), S = matrix(S, 10, 4)) # vectorization with recycling
#' SR(alpha = alpha, Rmax = Rmax, S = S, R_per_S = TRUE) # return recruits per spawner
#' 
#' curve(SR(SR_fun = "Ricker", alpha = alpha, Rmax = Rmax, S = x), from = 0, to = 2000,
#'       xlab = "Spawners", ylab = "Recruits", main = "Ricker")
#' abline(h = Rmax, lty = 2)
#' 
#' @seealso [salmonIPM()] for fitting models, [simIPM()] for simulating data
#' 
#' @export

SR <- function(SR_fun = "BH", alpha, Rmax, S, A = 1, R_per_S = FALSE) 
{
  stopifnot(is.character(SR_fun))
  stopifnot(is.numeric(alpha) | is.data.frame(alpha))
  stopifnot(is.numeric(Rmax) | is.data.frame(Rmax))
  stopifnot(is.numeric(S) | is.data.frame(S))
  stopifnot(is.numeric(A) | is.data.frame(A))
  
  if(is.data.frame(alpha)) alpha <- as.matrix(alpha)
  if(is.data.frame(Rmax)) alpha <- as.matrix(Rmax)
  if(is.data.frame(S)) alpha <- as.matrix(S)
  if(is.data.frame(A)) alpha <- as.matrix(A)
  if(SR_fun %in% c("B-H","bh","b-h")) SR_fun <- "BH"
  if(SR_fun == "ricker") SR_fun <- "Ricker"
  
  R <- switch(SR_fun,
              exp = alpha * S,
              BH = alpha * S / (1 + alpha * S / (A * Rmax)),
              Ricker = alpha * S * exp(-alpha * S / (A * exp(1) * Rmax)),
              stop("Error: ", SR_fun, " is not an available spawner-recruit function"))
  
  return(if(R_per_S) R/S else R) 
}
