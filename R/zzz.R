# R/zzz.R
# Environment with default prior parameters
.pkg_env <- new.env(parent = emptyenv())

#' @importFrom utils data
.onLoad <- function(libname, pkgname) {
  data("default_prior_settings", package = pkgname, envir = .pkg_env)
}
