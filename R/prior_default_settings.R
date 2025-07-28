# R/prior_default_settings.R

#' Quantile function for the half Student's T distribution
#'
#'
#' LaplacesDemon qhalft implementation does not work well for our purposes
#' here (there's some stochasticity in the results and their variance is quite
#' large)
#'
#' @param p as in [LaplacesDemon::qhalft()]
#' @param nu as in [LaplacesDemon::qhalft()]
#' @param scale as in [LaplacesDemon::qhalft()]
#'
#' @returns the quantile of the half Student's T with p = p
#' @keywords internal

qhalft <- function(p, nu, scale) {
  scale * stats::qt((p + 1) / 2, nu)
}

#' Prior Settings
#'
#' Database with default prior configurations.
#' @details It is used internally by [get_dfun_prior()] and [get_qfun_prior()]
#'   and can be obtained and modified using [user_priors()] and checked with
#'   [check_priors()]. It is automatically loaded on package load in zzz.R, and
#'   is generated during development using
#'   data_generation/make_prior_default_settings.R (only when needed, using a
#'   hash). The source object is default_prior_settings below.
#'
#' @format A named list of lists. The first level indicates the configuration
#'   settings (e.g., O for original), the second level indicates the parameter
#'   name, and the third , the dfun, qfun, and a list of parameters needed to
#'   run those functions. If the defaults are not valid, they should contain a
#'   valid element set to FALSE.
#' @keywords internal
"default_prior_settings"

default_prior_settings <- list(
  T = list(
    flipflop.lambda = list(
      dfun = "LaplacesDemon::dhalft",
      qfun = "phyfumr::qhalft",
      params = list(nu = 1.5, scale = 2.7)
    ),
    flipflop.mu = list(
      dfun = "LaplacesDemon::dhalft",
      qfun = "phyfumr::qhalft",
      params = list(nu = 1.5, scale = 0.13)
    ),
    clock.rate = list(
      dfun = "LaplacesDemon::dhalft",
      qfun = "phyfumr::qhalft",
      params = list(nu = 1.5, scale = 0.13)
    ),
    flipflop.gamma = list(
      dfun = "LaplacesDemon::dhalft",
      qfun = "phyfumr::qhalft",
      params = list(nu = 1.5, scale = 0.13)
    ),
    flipflop.MDbias = list(
      dfun = "stats::dlnorm",
      qfun = "stats::qlnorm",
      params = list(meanlog = 0, sdlog = 1)
    ),
    flipflop.rlambda = list(
      dfun = "stats::dlnorm",
      qfun = "stats::qlnorm",
      params = list(meanlog = 3, sdlog = 1.5)
    ),
    errorModel.kappaScale = list(
      dfun = "stats::dlnorm",
      qfun = "stats::qlnorm",
      params = list(meanlog = 4.56, sdlog = 0.3)
    ),
    errorModel.etaOffset = list(
      dfun = "stats::dbeta",
      qfun = "stats::qbeta",
      params = list(shape1 = 95, shape2 = 5)
    ),
    errorModel.deltaOffset = list(
      dfun = "stats::dbeta",
      qfun = "stats::qbeta",
      params = list(shape1 = 5, shape2 = 95)
    ),
    luca.branch = list(
      dfun = "stats::dunif",
      qfun = "stats::qunif",
      params = list(min = 0, max = Inf)
    ),
    constant.popSize = list(
      dfun = "stats::dlnorm",
      qfun = "stats::qlnorm",
      params = list(meanlog = 2.298218, sdlog = 2.148)
    )
  ),
  L = list(
    flipflop.lambda = list(
      dfun = "stats::dlnorm",
      qfun = "stats::qlnorm",
      params = list(meanlog = 2.077, sdlog = 0.672)
    ),
    flipflop.mu = list(
      dfun = "stats::dlnorm",
      qfun = "stats::qlnorm",
      params = list(meanlog = -2.528, sdlog = 0.672)
    ),
    clock.rate = list(
      dfun = "stats::dlnorm",
      qfun = "stats::qlnorm",
      params = list(meanlog = -2.528, sdlog = 0.672)
    ),
    flipflop.gamma = list(
      dfun = "stats::dlnorm",
      qfun = "stats::qlnorm",
      params = list(meanlog = -2.528, sdlog = 0.672)
    ),
    flipflop.MDbias = list(
      dfun = "stats::dlnorm",
      qfun = "stats::qlnorm",
      params = list(meanlog = 0, sdlog = 1)
    ),
    flipflop.rlambda = list(
      dfun = "stats::dlnorm",
      qfun = "stats::qlnorm",
      params = list(meanlog = 3, sdlog = 1.5)
    ),
    errorModel.kappaScale = list(
      dfun = "stats::dlnorm",
      qfun = "stats::qlnorm",
      params = list(meanlog = 4.56, sdlog = 0.3)
    ),
    errorModel.etaOffset = list(
      dfun = "stats::dbeta",
      qfun = "stats::qbeta",
      params = list(shape1 = 95, shape2 = 5)
    ),
    errorModel.deltaOffset = list(
      dfun = "stats::dbeta",
      qfun = "stats::qbeta",
      params = list(shape1 = 5, shape2 = 95)
    ),
    luca.branch = list(
      dfun = "stats::dunif",
      qfun = "stats::qunif",
      params = list(min = 0, max = Inf)
    ),
    constant.popSize = list(
      dfun = "stats::dlnorm",
      qfun = "stats::qlnorm",
      params = list(meanlog = 2.298218, sdlog = 2.148)
    )
  ),
  H = list(
    flipflop.lambda = list(
      dfun = "stats::dlnorm",
      qfun = "stats::qlnorm",
      params = list(meanlog = -0.770, sdlog = 1.241)
    ),
    flipflop.mu = list(
      dfun = "stats::dlnorm",
      qfun = "stats::qlnorm",
      params = list(meanlog = -3.835, sdlog = 1.241)
    ),
    clock.rate = list(
      dfun = "stats::dlnorm",
      qfun = "stats::qlnorm",
      params = list(meanlog = -3.835, sdlog = 1.241)
    ),
    flipflop.gamma = list(
      dfun = "stats::dlnorm",
      qfun = "stats::qlnorm",
      params = list(meanlog = -3.835, sdlog = 1.241)
    ),
    flipflop.MDbias = list(
      dfun = "stats::dlnorm",
      qfun = "stats::qlnorm",
      params = list(meanlog = 0, sdlog = 1)
    ),
    flipflop.rlambda = list(
      dfun = "stats::dlnorm",
      qfun = "stats::qlnorm",
      params = list(meanlog = 3, sdlog = 1.5)
    ),
    errorModel.kappaScale = list(
      dfun = "stats::dlnorm",
      qfun = "stats::qlnorm",
      params = list(meanlog = 4.56, sdlog = 0.3)
    ),
    errorModel.etaOffset = list(
      dfun = "stats::dbeta",
      qfun = "stats::qbeta",
      params = list(shape1 = 95, shape2 = 5)
    ),
    errorModel.deltaOffset = list(
      dfun = "stats::dbeta",
      qfun = "stats::qbeta",
      params = list(shape1 = 5, shape2 = 95)
    ),
    luca.branch = list(
      dfun = "stats::dunif",
      qfun = "stats::qunif",
      params = list(min = 0, max = Inf)
    ),
    constant.popSize = list(
      dfun = "stats::dlnorm",
      qfun = "stats::qlnorm",
      params = list(meanlog = 2.298218, sdlog = 2.148)
    )
  ),
  O = list(
    flipflop.lambda = list(
      dfun = "LaplacesDemon::dhalfnorm",
      qfun = "LaplacesDemon::qhalfnorm",
      params = list(scale = sqrt(pi / 2) / 1)
    ),
    flipflop.mu = list(
      dfun = "LaplacesDemon::dhalfnorm",
      qfun = "LaplacesDemon::qhalfnorm",
      params = list(scale = sqrt(pi / 2) / 0.05)
    ),
    clock.rate = list(
      dfun = "LaplacesDemon::dhalfnorm",
      qfun = "LaplacesDemon::qhalfnorm",
      params = list(scale = sqrt(pi / 2) / 0.05)
    ),
    flipflop.gamma = list(
      dfun = "LaplacesDemon::dhalfnorm",
      qfun = "LaplacesDemon::qhalfnorm",
      params = list(scale = sqrt(pi / 2) / 0.05)
    ),
    flipflop.MDbias = list(
      dfun = "stats::dlnorm",
      qfun = "stats::qlnorm",
      params = list(meanlog = 0, sdlog = 1)
    ),
    flipflop.rlambda = list(
      dfun = "stats::dlnorm",
      qfun = "stats::qlnorm",
      params = list(meanlog = 3, sdlog = 1.5)
    ),
    errorModel.kappaScale = list(
      dfun = "stats::dlnorm",
      qfun = "stats::qlnorm",
      params = list(meanlog = 4.56, sdlog = 0.3)
    ),
    errorModel.etaOffset = list(
      dfun = "stats::dbeta",
      qfun = "stats::qbeta",
      params = list(shape1 = 95, shape2 = 5)
    ),
    errorModel.deltaOffset = list(
      dfun = "stats::dbeta",
      qfun = "stats::qbeta",
      params = list(shape1 = 5, shape2 = 95)
    ),
    luca.branch = list(
      dfun = "stats::dunif",
      qfun = "stats::qunif",
      params = list(min = 0, max = Inf),
      valid = F
    ),
    constant.popSize = list(
      dfun = "stats::dlnorm",
      qfun = "stats::qlnorm",
      params = list(meanlog = 2.298218, sdlog = 2.148)
    )
  )
)
