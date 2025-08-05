# R/model_default_params.R

#' Model Settings
#'
#' Database with default model parameter values based on empirical data
#'
#' @details It is used internally by [qc_csv()] and can be obtained and modified
#'   using [get_model_params()]. It is automatically loaded on package load
#'   in zzz.R.
#'
#' @format A named list with values for
#'   errorModel.deltaOffset,errorModel.etaOffset, errorModel.kappaScale, and
#'   alignment.stemCells
#' @keywords internal
"default_model_params"

#TODO: re-generate this object with the latest version of the analysis
