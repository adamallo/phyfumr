# R/posterior_prior_comparison.R

#Code related to the comparison of posterior and prior volumes

#Numerical volume estimation

#' Finds a density such that the cumulative probability of an empirical density
#' function is above a given value
#'
#' @details KDE-based
#'
#' @param x a vector of samples from the distribution of interest
#' @param prob cumulative probability at which to find the density threshold
#' @param ... other params for [ks::kde()]
#'
#' @returns a list with the density threshold and the kde object

get_kde_threshold <- function(x, prob = 0.95, ...) {
  kde_obj <- ks::kde(x = matrix(x, ncol = 1), ...)
  dx <- diff(kde_obj$eval.points[1:2])
  sorted_density <- sort(kde_obj$estimate, decreasing = TRUE)
  cum_mass <- cumsum(sorted_density * dx)
  threshold <- sorted_density[which(cum_mass >= prob)[1]]
  return(list(threshold = threshold, kde = kde_obj))
}

#' Estimates the volume above a density threshold in an empirical density function
#'
#' @param kde_obj [ks::kde()] object
#' @param d density threshold
#'
#' @seealso [an_volume_above()] for its analytical counterpart
#'
#' @returns the probability of the random variable being between the bounds
#'   determined by d
#' @export

est_volume_above <- function(kde_obj, d) {
  dx <- diff(kde_obj$eval.points[1:2])
  sum(kde_obj$estimate[kde_obj$estimate > d]) * dx
}

#Analytical volume calculation

#' Volume of a PDF above a given density (analytical)
#'
#' This function returns the probability of a continuous random variable being
#' between two bounds determined by a given density. It assumes an unimodal
#' distribution
#'
#' @param dfun density function
#' @param qfun quantile function
#' @param d threshold density
#' @param support support
#' @param tol tolerance from asymptotes
#' @param ... parameters for qfun and dfun
#' @seealso [est_volume_above()] for its numerical counterpart
#'
#' @returns the probability of the random variable being between the bounds
#'   determined by d
#' @export

an_volume_above <- function(dfun, qfun, d, support = c(-Inf, Inf), tol = 1e-8, ...){
  bounds <- an_support_above(dfun, qfun, d, support, tol, ...)
  if (all(!is.na(bounds))) return (stats::integrate(dfun, lower = bounds[1], upper = bounds[2], ...)$value)
  else return(0)
}


#' Probability interval of a PDF over a given density (analytical)
#'
#' Assumes the PDF is unimodal and continuous
#'
#' @param dfun density function
#' @param qfun quantile function
#' @param d threshold density
#' @param support support
#' @param tol tolerance from asymptotes
#' @seealso [an_volume_above()] which finds the PDF's volume within this
#'   interval
#'
#' @returns A 2-element array with the lower and upper values of the interval.
#'   NAs if the PDF volume over d is 0.

an_support_above <- function(dfun, qfun, d, support = c(-Inf, Inf), tol = 1e-8) {

  # Auto-support inference (fallback to quantiles if Inf)
  auto_support <- c(
    if (is.finite(support[1])) support[1] else qfun(tol),
    if (is.finite(support[2])) support[2] else qfun(1 - tol)
  )

  # Check density at mode
  x_mode <- stats::optimize(function(x) -dfun(x), interval = auto_support)$minimum
  d_mode <- dfun(x_mode)

  if (d_mode < d) {
    warning("Prior mass = 0, returning empty bounds")
    return(c(NA, NA))
  }

  # Get quantile limits from support
  p_lower <- if (is.finite(support[1])) suppressWarnings(
    stats::uniroot(function(p) qfun(p) - auto_support[1],
            interval = c(tol, 1 - tol))$root
  ) else tol

  p_upper <- if (is.finite(support[2])) suppressWarnings(
    stats::uniroot(function(p) qfun(p) - auto_support[2],
            interval = c(tol, 1 - tol))$root
  ) else 1 - tol


  # Find roots safely
  find_bound <- function(p_from, p_to) {
    tryCatch({
      root <- stats::uniroot(function(p) dfun(qfun(p)) - d,
                      lower = p_from, upper = p_to, extendInt = "no")$root
      qfun(root)
    }, error = function(e) NA)
  }

  x1 <- find_bound(p_lower, 0.5)
  x2 <- find_bound(0.5, p_upper)

  if (is.na(x1) && is.na(x2)) {
    warning("Could not determine support interval for density > d.")
    return(c(NA, NA))
  } else if (is.na(x1)) {
    return(c(auto_support[1], x2))  # e.g. mode at left edge
  } else if (is.na(x2)) {
    return(c(x1, auto_support[2]))  # e.g. mode at right edge
  } else {
    return(c(x1, x2))
  }
}

#Others

#' Calculates the Kullback–Leibler (KL) divergence from an empirical distribution to a probability distribution with known density function
#' @param posteriorKDE [ks::kde()] object with the KDE estimation of the posterior
#'   samples
#' @param dfun density function that takes a vector of values and returns their
#'   density
#'
#' @returns KL divergence in bits
#' @export
kl_divergence_prior_to_posterior <- function(posteriorKDE,dfun) {

  #Density estimates from KDE
  x_vals <- posteriorKDE$eval.points
  p_post <- posteriorKDE$estimate

  #Prior density of the same points from analytical function
  p_prior <- dfun(posteriorKDE$eval.points)

  # Avoid log(0) or division by 0
  valid <- (p_post > 0) & (p_prior > 0)

  #Calculation
  kl <- sum(p_post[valid] * log2(p_post[valid] / p_prior[valid])) * mean(diff(posteriorKDE$eval.points))

  return(kl)
}

#' Compares the densities of two distributions to deem them comparable
#'
#' Used in a Bayesian context to compare Prior and Posterior distributions.
#' Empirically, we found that if the density threshold for a high credible
#' interval (e.g. 95%) does not have any volume in the prior (i.e., the prior is
#' much more disperse than the posterior), the data has informed the posterior
#' enough and the model is quite robust to changes in the prior.
#'
#' @param x a vector of samples from the distribution of interest
#' @param dfun density function
#' @param qfun quantile function
#' @param prob cumulative probability at which to find the density threshold
#' @param support support
#' @param tol tolerance from asymptotes
#' @param ... other params for [ks::kde()]
#' @seealso [compare_e_a()] to actually calculate the amount of change in volume
#'   over a density threshold (continuous take of this binary statistic)
#'
#' @returns logical indicating if the distributions are comparable (True) or not
#'   (False)
#' @export

is_dispersion_comparable_e_a <- function(x,dfun,qfun,prob=0.95,support = c(-Inf, Inf),tol=1e-8,...){

  #Empirical density threshold of the objective distribution
  d <- get_kde_threshold(x,prob,...)$threshold

  # Auto-support inference (fallback to quantiles if Inf)
  auto_support <- c(
    if (is.finite(support[1])) support[1] else qfun(tol),
    if (is.finite(support[2])) support[2] else qfun(1 - tol)
  )

  #Maximal analytical density of the reference distribution
  ref_mode <- stats::optimize(function(x) -dfun(x), interval = auto_support)$minimum
  d_ref_mode <- dfun(ref_mode)

  return(d_ref_mode >= d)
}

#' Compares an empirical and analytical distribution
#'
#' Used in a Bayesian context to compare Prior and Posterior distributions.
#'
#' @param x a vector of samples from the distribution of interest
#' @param dfun density function
#' @param qfun quantile function
#' @param prob cumulative probability at which to find the density threshold
#' @param support support
#' @param tol tolerance from asymptotes
#' @param ... other params for [ks::kde()]
#' @seealso [is_dispersion_comparable_e_a] for a binary take of this metric, and
#'   [kl_divergence_prior_to_posterior] for the other metric used here
#'
#' @returns list with the posterior volume, prior volume, density threshold,
#'   volume ratio, update info (log2 transformation of the former), and
#'   Kullback–Leibler (KL) divergence from the analytical distribution to the
#'   empirical (False)
#' @export

compare_e_a <- function(x, dfun, qfun, prob = 0.95, support = c(-Inf, Inf), tol=1e-8, ...) {
  # Estimate density threshold
  kdeResult <- get_kde_threshold(x, prob,...)
  threshold <- kdeResult$threshold

  # Prior volume above threshold (analytically)
  prior_vol <- an_volume_above(dfun,qfun,threshold,support)

  # KL divergence
  kldiv <- kl_divergence_prior_to_posterior(kdeResult$kde,dfun)

  list(
    posterior_volume = prob,
    prior_volume = prior_vol,
    volume_density_threshold = threshold,
    volume_ratio = prob/prior_vol,
    update_info = log2(prob/prior_vol),
    kl_divergence = kldiv
  )
}

