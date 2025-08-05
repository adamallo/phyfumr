# R/zzz.R
# Package environment to cache data and store defaults
.phyfumr_env <- new.env(parent = emptyenv())

# list to store CSV data
.phyfumr_env[["loadedCSVs"]] <- list()

# Making the default prior parameters available in the .pkg_env
.onLoad <- function(libname, pkgname) {
  utils::data("default_prior_settings", package = pkgname, envir = .phyfumr_env)
  utils::data("default_model_params", package = pkgname, envir = .phyfumr_env)
}

# Silencing devtools::check with data.table syntax
utils::globalVariables(c("loadedCSVs",
                         "valid",
                         "patient",
                         "samplename",
                         "locus"))
