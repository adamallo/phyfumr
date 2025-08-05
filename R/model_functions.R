# R/model_functions.R

#' Returns the default model settings for user exploration and modification
#'
#' @returns A named list with values for
#'   errorModel.deltaOffset,errorModel.etaOffset, errorModel.kappaScale, and
#'   alignment.stemCells
#' @seealso [qc_csv()]
#' @export

get_model_params <- function(){
  return(.phyfumr_env$default_model_params)
}
