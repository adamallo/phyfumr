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

#' Run a bundled script
#'
#' Used as Rscript -e 'phyfumr::run("scriptName.R")' --args arg1 arg2
#'
#' @param script Name of the script (file in `exec/`)
#' @export
#'

run <- function(script) {
  path <- system.file("exec", script, package = "phyfumr")
  if (path == "") stop("Script not found: ", script, call. = FALSE)

  # capture all args after --args
  args <- commandArgs(trailingOnly = TRUE)

  # run the bundled script with those args
  res <- system2(command = R.home("bin/Rscript"),
                 args = c(path, args))
  invisible(res)
}
