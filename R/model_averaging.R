# R/integrate_over_S.R

#' Resample posterior using a sample map
#'
#' @param samples list of posterior samples
#' @param map 2D array with a row per final posterior sample and column per
#'   input posterior sample, with elements indicating the ith position of the
#'   input posterior sample to be used as output, or NA otherwise. Only one
#'   non-na element per row is expected
#'
#' @returns posterior sample (array)
#' @keywords internal

apply_samplemap <- function(samples,map){
  rowSums(sapply(1:length(samples),function(i){
    samples[[i]][map[,i]]
    }),na.rm=T)
}

#' Calculate posterior probabilities given all log joint probabilities
#'
#' @param ljps array of the log joint probabilities (log marginal likelihoods +
#'   log prior) of each of the posterior samples (x)
#'
#' @return probabilities
#' @keywords internal
calculate_model_averaging_posteriors <- function(ljps) {
  log_marginal <- matrixStats::logSumExp(ljps)
  log_probs <- ljps-log_marginal
  exp(log_probs)
}



resample_softmax <- function(trees_files,
                             log_files = NULL,
                             ...){
message("Loading traces...")
sink(nullfile()) #shut up!
chains <- invisible(lapply(1:length(trees_files),FUN=function(i_file){
  fix_nexus(trees_files[i_file],...)
  chain <- rwty::load.trees(trees_files[i_file],format="BEAST",logfile = log_files[i_file]);
  chain}))
sink()
message("Done\n")

}

resample_softmax_chains <- function(chains,
                             ljps = NULL,
                             lpps = NULL,
                             burnin_p = 0){

  chains <- ensure_list_chains(chains)

  ##TODO remove burnin

  if(!is.null(ljps) && length(chains) != length(ljps))
     stop("ERROR: ljps must be NULL (calculated here) or of the same length as chains")

  if(!is.null(ljps)){
    probs <- calculate_model_averaging_posteriors(ljps)
  } else {    ##estimate using HME
    lMLEs <- sapply(chains,function(chain){
      if(is.null(chain[["ptable"]]) || is.null(chain[["ptable"]][,"likelihood"]))
         return(NaN)
      LaplacesDemon::LML(LL=chain$ptable$likelihood,method="HME")$LML})
    if(!is.null(lpps)){
      if(length(lpps) != length(ljps))
        stop("ERROR: lpps must be of the same length as ljps and chains")
      lpjs <- lMLEs + lpps
    } else {
      lpjs <- lMLEs
      warning("lMLEs are calculated here and log prior probabilities for each model were not provided, thus, we are assuming an improper uniform prior")
    }
    probs <- calculate_model_averaging_posteriors(ljps)
  }

  nSamples <- unique(sapply(chains,FUN=function(x){length(x$trees)}))



  return(resampled_chains)
}

#' Resample posteriors according to probabilities
#'
#' @param x list of posterior samples. They do not have to have the same length
#' @param probs array of sampling probabilities (one per posterior sample x)
#' @param burnin_p proportion of samples in each posterior sample to be removed
#'   as burnin
#' @param method one of "jump" or "thin". Jump requires all posterior samples to
#'   be of the same length
#'
#' @returns a list with the posterior sample and a 2D array with a row per final
#'   posterior sample and column per input posterior sample, with elements
#'   indicating the ith position of the input posterior sample to be used as
#'   output, or NA otherwise. Only one non-na element per row is expected
#' @keywords internal

resample_posterior_array <- function(x,
                                     probs,
                                     burnin_p=0.1,
                                     method=c(NULL,"jump","thin")){
  n_x <- length(x)
  method <- match.arg(method)

  if(length(x) != length(probs))
    stop("ERROR: probs must be of the same length as chains")

  valid_x <- lapply(x,function(x,burnin_p){
    last_i <- length(x)
    first_i <- round(last_i*burnin_p)
    x[first_i:last_i]
  },burnin_p=burnin_p)

  lengths <- sapply(valid_x,length)

  if(length(unique(lengths))>1 | method=="thin") {
    if(method=="jump")
      warning("The selected method \"jump\" will not be used since posterior samples are of different lenght")

    max_length <- max(lengths)
    which_x <- t(stats::rmultinom(max_length,1,probs)) #resampling
    cum_samples <- apply(which_x,2,cumsum) #number of samples needed from each posterior sample (col) at each final possible length (row)
    last_sample <- min(sapply(1:ncol(cum_samples),function(i){suppressWarnings(min(which(cum_samples[,i] == lengths[i])))})) #length of the maximum length possible given the sampling
    n_samples <- cum_samples[last_sample,] #final number of samples per posterior sample
    final_x <- which_x[1:last_sample,]#final samples from each posterior sample

    which_i <- sapply(1:n_x,function(i){#subsampling each of the posteriors to get the appropriate number of samples (thinning)
      round(seq(from=1,to=lengths[i],length.out=n_samples[i]))})

    final_map <- sapply(1:n_x,function(i){#applying the sampled positions
      this_x <- ifelse(final_x[,i]==0,NA,final_x[,i])#using NAs because 0 could be a valid sample
      this_x[!is.na(this_x)] <- which_i[[i]]
      this_x})

    return(list(sample=get_samples(valid_x,final_map),
                map=final_map))
  } else {
    #resample jump
    #We set the most probable posterior sample as default, and copy values from other chains where needed based on their probabilities
    def_x <- which(probs==max(probs))
    if(length(def_x)>1){#If several have the same
      def_x <- sample(x=def_x,size=1)
    }
    out_x <- x[[def_x]]

    #Replace when needed
    x_i <- apply(stats::rmultinom(lengths[1],1,probs),2,FUN=function(x){which(x==1)})
    for (i in 1:length(x_i)) {
      if(x_i[i] != def_x) {
        out_x[i] <- x[[out_x[i]]][i]
      }
    }
    which_i <- 1:length(valid_x)
    return(list(sample=out_x,map=t(sapply(1:length(x_i),function(x)ifelse(x_i[x]==which_i,x,NA)))))
  }
}



