# R/aaa-utils.R

#' Uses pblapply or lapply depending on their availability
#'
#' @inherit base::lapply
#' @inheritDotParams pbapply::pblapply cl
#'
#' @keywords internal

slapply <- function(X, FUN, ...) {
  if(is.null(.phyfumr_env[["use_parallel"]]) || !.phyfumr_env[["use_parallel"]]){
    ellipsis_args <- list(...)
    in_args <- c(ellipsis_args,list(X = X, FUN = FUN))
    in_args[["cl"]] <- NULL
    return(do.call(lapply,in_args))
  } else {
    return(pbapply::pblapply(X = X, FUN = FUN, ...))
  }
}
