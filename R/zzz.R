# R/zzz.R

# Package environment to cache data and store defaults
.phyfumr_env <- new.env(parent = emptyenv())

# Making the default prior parameters available in the .pkg_env
.onLoad <- function(libname, pkgname) {
  utils::data("default_prior_settings", package = pkgname, envir = .phyfumr_env)
  utils::data("default_model_params", package = pkgname, envir = .phyfumr_env)

  # Checking if parallelism using pbapply is available
  .phyfumr_env[["use_parallel"]] <- requireNamespace("pbapply", quietly = TRUE)

  # Setting rwty.processors to 1 by default since we may call their functions in
  # parallel
  rwty.processors <<- 1

  # list to store CSV data
  .phyfumr_env[["loadedCSVs"]] <- list()

  # Alternative phyfum parameterizations
  .phyfumr_env[["absolute_rates_params"]] <- c("flipflop.mu","flipflop.lambda","flipflop.gamma")
  .phyfumr_env[["relative_rates_params"]] <- c("clock.rate","flipflop.rlambda","flipflop.MDbias")

  # Phyfum log output information
  .phyfumr_env[["not_params"]] <- c("chain","state","treedistance")

  # System detection
  .phyfumr_env[["sys"]] <- Sys.info()[["sysname"]]
  .phyfumr_env[["has_cairo"]] <- capabilities("cairo")

  # To use system's tail for efficiency, not in Windows
  .phyfumr_env[["system_tail"]] <- Sys.which("tail") != ""

  # Avoiding spurious Rplots.pdf
  if(!interactive()) grDevices::pdf(NULL)
}

# Silencing devtools::check with data.table syntax
utils::globalVariables(c("loadedCSVs",
                         "valid",
                         "patient",
                         "samplename",
                         "locus",
                         "chain",
                         "param",
                         "meanP",
                         "medianP",
                         "state",
                         "treedistance",
                         "x",
                         ".",
                         "burnin",
                         "likelihood",
                         "flipflop.mu",
                         "flipflop.lambda",
                         "flipflop.gamma",
                         "clock.rate",
                         "flipflop.rlambda",
                         "flipflop.MDbias",
                         "problem",
                         "Rhat",
                         "N",
                         "best",
                         "cond",
                         "extreme_best",
                         "extreme_valid",
                         "lML",
                         "max_s",
                         "method_order",
                         "minY",
                         "min_s",
                         "s",
                         "selected",
                         "valid_method"))

