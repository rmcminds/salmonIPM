#' Plot prior and posterior distributions of model hyperparameters
#'
#' Plot posterior draws, analytical prior distributions, and (optionally) true
#' values of hyperparameters in a fitted **salmonIPM** model.
#'
#' @name plot_prior_posterior
#'
#' @param object An object of class [salmonIPMfit] with prior information stored
#'   in `prior.info`.
#' @param pars Character vector specifying hyperparameters to plot. The default
#'   is all top-level hyperparameters, i.e. those parameters that are given
#'   priors.
#' @param include Logical scalar defaulting to `TRUE` indicating whether to
#'   include or exclude the parameters given by `pars`. If `FALSE`, only entire
#'   multidimensional parameters can be excluded, rather than particular
#'   elements of them. This is most likely to be useful for excluding large
#'   correlation matrices from the plot.
#' @param true Named list containing true hyperparameter values used to generate
#'   the pseudo-data to which `object` was fitted. See the `pars` argument to
#'   [simIPM()] for the structure of this list.
#'
#' @return A [ggplot] object. If the result is not assigned to an object, it
#'   will be automatically plotted.
#'
#' @details For each hyperparameter in `pars`, the posterior draws are plotted
#'   as a density histogram with the prior PDF overlaid. If `true` is provided,
#'   the known parameter value used to simulate the pseudo-data for the model
#'   fit is shown as a vertical line.
#'
#' @seealso [prior_summary()], [priors], [salmonIPM()], [salmonIPMfit]
#'
#' @importFrom dplyr %>% select mutate across matches rename_with reframe group_by
#' @importFrom tidyr pivot_longer replace_na
#' @importFrom posterior as_rvar mutate_variables
#' @importFrom distributional dist_normal dist_lognormal dist_beta dist_uniform
#'   dist_truncated
#' @importFrom ggplot2 ggplot aes after_stat geom_histogram geom_line geom_vline
#'   scale_x_continuous labs alpha margin element_blank element_text facet_wrap
#'   vars theme theme_classic
#' @importFrom ggdist dlkjcorr_marginal rlkjcorr_marginal plkjcorr_marginal
#'   qlkjcorr_marginal
#' @export

plot_prior_posterior <- function(object, pars = NULL, include = TRUE, true = NULL) {
  stopifnot("salmonIPMfit" %in% class(object))
  prior.info <- object$prior.info
  if(is.null(pars)) {
    pars <- names(prior.info)
  } else {
    pars <- if(include) pars else setdiff(names(prior.info), pars)
  }

  # Posterior
  drf <- as_draws_df(object, pars = pars) %>% as.data.frame() %>% 
    rename_with(~ sub("(^alpha.+|^Rmax.+|^Mmax.+)", "log(\\1)", .x)) %>%
    mutate(across(contains("log("), log)) %>%
    draws_df_lower_tri()  # only lower triangle of correlation matrices

  post <- drf %>% 
    pivot_longer(cols = !starts_with("."), names_to = "par", values_to = "value") %>% 
    mutate(par = factor(par, levels = unique(par)))

  # Priors
  priors <- lapply(levels(post$par), function(par) {
    .par <- gsub("log\\((.+)\\)", "\\1", gsub("\\[.*\\]", "", par))
    .indx <- suppressWarnings(as.numeric(gsub(".+\\[.*(\\d+)\\]", "\\1", par)))
    prinfo <- prior.info[[.par]]
    dist <- prinfo$dist
    args <- prinfo[setdiff(names(prinfo), "dist")]
    args <- switch(dist, 
                   uniform = list(min = args$lb, max = args$ub),
                   lognormal = list(mu = args$meanlog, sigma = args$sdlog),
                   beta = list(shape1 = args$a, shape2 = args$b),
                   dirichlet = list(shape1 = args$concentration[.indx],
                               shape2 = sum(args$concentration[-.indx])),
                   lkj_corr = c(args, K = object$stanfit@par_dims[[.par]][2]),
                   args)
    dist <- switch(dist, lognormal = "normal", dirichlet = "beta", 
                   lkj_corr = "lkjcorr_marginal", dist)
    pr <- do.call(paste0("dist_", dist), args)
    bounds <- replace_na(attr(prinfo, "bounds"), Inf)
    if(!is.null(bounds)) pr <- dist_truncated(pr, lower = bounds[1], upper = bounds[2])
    return(pr)
  })
  
  priors <- data.frame(par = unique(post$par), prior = do.call("c", priors)) %>% 
    mutate(x0 = apply(drf[,par], 2, min), x1 = apply(drf[,par], 2, max)) %>% 
    group_by(par) %>% reframe(x = seq(x0, x1, length = 100), px = unlist(density(prior, x)))
  
  # True values
  if(is.null(true)) {
    true <- data.frame(par = unique(post$par), value = as.numeric(NA))
  } else {
    if("sigma_R" %in% pars && !object$pool_pops && !is.null(true$sigma_year_R)) 
      true$sigma_R <- true$sigma_year_R
    true <- lapply(true[pars], as_rvar) %>% as_draws_df() %>% as.data.frame() %>%
      draws_df_lower_tri() %>% setNames(names(drf)) %>%
      mutate(across(contains("log("), log)) %>%
      pivot_longer(!starts_with("."), names_to = "par", values_to = "value") %>%
      mutate(par = factor(par, levels = unique(par)))
  }
  
  # Plot
  gg <- post %>% ggplot() + 
    geom_histogram(aes(value, after_stat(density)), bins = 20, 
                   col = "white", fill = alpha("slategray4", 0.5)) +
    geom_line(data = priors, aes(x, px), col = "black", lwd = 0.8) +
    geom_vline(data = true, aes(xintercept = value), na.rm = TRUE, lwd = 0.8, lty = 2) +
    facet_wrap(vars(par), scales = "free", strip.position = "bottom") + 
    theme_classic() + 
    theme(strip.placement = "outside", strip.background = element_blank(),
          strip.text = element_text(size = 11, margin = margin(b = 3, t = 0)),
          axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text = element_text(size = 10), axis.text.y = element_blank()) +
    labs(x = "", y = "")
  
  gg
}

#' Wrapper for [gnorm::gnorm()] distribution using [distributional::dist_wrap()]
#'
#' @inheritParams priors
#'
#' @importFrom distributional dist_wrap
dist_gnormal <- function(mean = 0, scale = 1, shape = 1) {
  dist_wrap("gnorm", mu = mean, alpha = scale, beta = shape)
}

#' Wrapper for [ggdist::lkjcorr_marginal()] distribution using
#' [distributional::dist_wrap()]
#'
#' @param K Dimension of the correlation matrix. Must be greater than or equal to 2.
#' @inheritParams priors
#'
#' @importFrom distributional dist_wrap
dist_lkjcorr_marginal <- function(K, eta = 1) {
  dist_wrap("lkjcorr_marginal", K = K, eta = eta)
}
