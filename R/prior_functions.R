# R/prior_functions.R

#' Returns the default Prior settings for user exploration and modification
#'
#' @details The default settings have improper luca.branch priors so this function is needed to modify it if needed
#'
#' @returns list of lists with prior settings. The first level indicates the named prior configuration, the second the parameter name.
#' @export

user_priors <- function(){
  return(.phyfumr_env$default_prior_settings)
}

#' Finds invalid prior settings
#'
#' Returns the list of configs with a list of invalid priors within
#'
#' @param prior_settings in case you need to check user prior configuration instead of the default (edge use-case)
#'
#' @returns list of configs with a list of invalid priors within
#' @export

check_priors <- function(prior_settings=NULL){
  if(is.null(prior_settings)){
    these_priors <- prior_settings
  } else {
    these_priors <- user_priors()
  }
  Filter(length, lapply(these_priors, function(config) {
    names(Filter(function(param) isFALSE(param$valid), config))
  }))
}

#' Internal prior-setting checking tool for a given parameter and prior config
#'
#' @param param name of the parameter
#' @param config prior settings used, O by default
#' @param prior_settings user-specified prior settings, normally set by
#'   modification of the default prior settings using [user_priors()]
#' @seealso [check_priors()]
#'
#' @returns prior_settings

check_param_prior <- function(param,config="O",prior_settings=NULL){
  these_settings <- NULL
  if(!is.null(prior_settings)) {
    these_settings <- prior_settings[[config]][[param]]
  }
  if(is.null(these_settings)){
    if(!is.null(prior_settings))
      warning("User-specified prior settings were not set properly and will not be used, using the default instead")
    these_settings <- .phyfumr_env$default_prior_settings[[config]][[param]]
  }
  if(is.null(these_settings))
    return(function(x) NA_real_)
  if(!is.null(these_settings[["valid"]]) && these_settings[["valid"]]==F){
    stop(paste0("This prior is not valid.",
                " Use user_priors to obtain the default priors and modify them accordingly.",
                "You can obtain a list of invalid priors using the invalid_priors function"))
  }
  return(these_settings)
}

#' Obtains function from function name including package/environment name
#'
#' @param x the name of the function. E.g., stats::qfun
#'
#' @returns function

get_function <- function(x){
  namespace <- gsub("^([^:]*):::?[^:]*$","\\1",x)
  function_name <- gsub("^[^:]*:::?([^:]*)$","\\1",x)
  return(get(x = function_name,envir = asNamespace(namespace)))
}

#' Returns the density function of the Prior for a given parameter
#'
#' @inheritParams check_param_prior
#' @seealso [check_priors()]
#'
#' @return function that returns the density of the parameter value x, or
#'   function that returns NA_real_ if it is unknown
#' @export

get_dfun_prior <- function(param,config="O",prior_settings=NULL) {
  these_settings <- check_param_prior(param,config,prior_settings)
  dfun <- get_function(these_settings$dfun)
  #Preventing lazy evaluation issues
  force(these_settings)
  force(dfun)
  if(!is.null(dfun)){
    return(function(x) do.call(dfun,c(list(x=x),these_settings$params)))
  } else {
    return(NULL)
  }
}

#' Returns the quantile function of the Prior for a given parameter
#'
#' @inherit get_dfun_prior params seealso
#'
#' @return function that returns the quantile of the parameter probability p, or
#'   function that returns NA_real_ if it is unknown
#' @export

get_qfun_prior <- function(param,config="O",prior_settings=NULL) {
  these_settings <- check_param_prior(param,config,prior_settings)
  qfun <- get_function(these_settings$qfun)
  #Preventing lazy evaluation issues
  force(these_settings)
  force(qfun)
  if(!is.null(qfun)){
    return(function(p) do.call(qfun,c(list(p=p),these_settings$params)))
  } else {
    return(NULL)
  }
}
