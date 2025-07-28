# R/distance_to_origin.R

#Code related to the calculation of distances between LUCA and the samples

#' Calculates empirical frequencies given the tip likelihoods
#'
#' @inheritParams get_entropy
#' @seealso [get_peak_likelihoods()]
#'
#' @returns empirical frequencies
#' @export

get_empirical_freqs <- function(likelihoods,ml=F) {
  flat_likelihoods <- do.call(rbind, apply(likelihoods, MARGIN = 3, identity, simplify = F)) # we calculate frequencies across samples and loci
  flat_probs <- flat_likelihoods / rowSums(flat_likelihoods) # Transforms state (theoretical b-value peaks) likelihoods into peak assignment probabilities
  cleanflat_probs <- flat_probs[stats::complete.cases(flat_probs),,drop=F] #complete.cases checks by row (i.e., if any element of the row is missing the row is missing), returns a array[nrow]
  if(ml==F){
    freqs <- colSums(cleanflat_probs)/nrow(cleanflat_probs)
  } else {
    freqs <- colSums((cleanflat_probs==matrixStats::rowMaxs(cleanflat_probs))+0)/nrow(cleanflat_probs)
  }
  return(freqs)
}

#' Calculates the distance between the origin and the data
#'
#' It uses the Earth's Mover Distance between the assumed starting frequencies
#' (0.5 each for fully methylated and fully demethylated) and the empirical
#' frequencies
#'
#' @details The relative Earth's Mover Distance is normalized by the maximum
#'   theoretical distance. Importantly, this is an overestimation since, given a
#'   specific phyfum model, the maximum expected distance is that of the
#'   equilibrium frequencies of the substitution process. We may add this in the
#'   near future but we have not for now.
#'
#' @param freqs empirical frequencies
#'
#' @returns list with the Earth's Mover Distance (emd) and the relative Earth's
#'   Mover Distance (rEMD)
#'
emd_to_origin <- function(freqs){
  freqBs <- seq(from=0,to=1,length.out=length(freqs))
  origin_freqs <- rep(0,length.out=length(freqs))
  origin_freqs[1] <- 0.5
  origin_freqs[length(origin_freqs)] <- 0.5
  emd <- emdist::emdw(freqBs,origin_freqs,freqBs,freqs)
  maxDivFreq <- rep(0,length.out=length(freqs)) #this is the maximum divergence frequency under any model, not this specific model
  midFreq <- length(freqs)/2
  maxDivFreq[ceiling(midFreq)] <- 1
  return(list(emd=emd,rEMD=emd/emdist::emdw(freqBs,origin_freqs,freqBs,maxDivFreq)))
}
