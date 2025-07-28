# R/entropy.R

#Code related to the estimation of entropy under the phyfum

#' Calculates the b-peak likelihoods given a fCpG dataset and phyfum model
#'
#' @param b_data b-value alignment (i.e., matrix (or 2-dimensional array) with a
#'   column per fCpG site and row per sample)
#' @param delta flipflop.delta parameter of the phyfum model (methylated
#'   background noise)
#' @param eta flipflop.eta parameter of the phyfum model (de-methylated
#'   background negative noise)
#' @param kappa flipflop.kappa parameter of the phyfum model (b concentration)
#' @param s s parameter of the phyfum model (number of stem cells)
#'
#' @returns 3D array with \code{[sample,peak,locus]} dimmentions with the b-peak
#'   likelihoods for each sample and fCpG locus
#' @export

get_peak_likelihoods <- function(b_data,delta,eta,kappa,s){
  #state_count <- (0.5 * (s + 1) * (s + 2))
  peak_count <- 2 * s + 1
  loci_count <- nrow(b_data)
  sample_count <- ncol(b_data)

  ideal_beta <- NULL
  beta <- NULL
  transformed_alpha <- NULL
  transformed_beta <- NULL
  lgamma_alpha <- NULL
  lgamma_beta <- NULL
  lgamma_alphaplusbeta <- NULL

  for (peak in 1:peak_count){
    ideal_beta[peak] = (peak - 1) / (peak_count -1)
    beta[peak] = (eta - delta) * ideal_beta[peak] + delta
    transformed_alpha[peak] = beta[peak] * kappa;
    transformed_beta[peak] =  (1 - beta[peak]) * kappa;
    lgamma_alpha[peak] = lgamma(transformed_alpha[peak]);
    lgamma_beta[peak] = lgamma(transformed_beta[peak]);
    lgamma_alphaplusbeta[peak] = lgamma(transformed_alpha[peak] + transformed_beta[peak]);
  }

  likelihoods <- array(NA,dim = c(loci_count,peak_count,sample_count))
  for (locus in 1:loci_count){
    for (peak in 1:peak_count){
      for (sample in 1:sample_count){
        likelihoods[sample,peak,locus] = exp((transformed_alpha[peak] - 1) * log(b_data[sample,locus]) + (transformed_beta[peak] - 1) * log1p(-b_data[sample,locus]) - lgamma_alpha[peak] - lgamma_beta[peak] + lgamma_alphaplusbeta[peak]);
      }
    }
  }

  return(likelihoods)
}

#' Calculates the entropy of a single locus using the marginalized approach
#'
#' @param mat matrix \code{[sample,peak]}
#'
#' @returns the entropy of a single locus
#' @seealso [get_entropy()]
#' @keywords internal

calculate_locus_entropy <- function(mat) {
  rSums <- rowSums(mat)
  if(any(rSums==0)) #This site is not valid because a b value is not valid with the current model
    return(NA)
  freqs <- colSums(mat / rSums)/nrow(mat)
  return(-sum(freqs*log2(freqs)))
}

#' Calculates the entropy of a single locus using the maximum likelihood peak
#' assignment
#'
#' @inherit calculate_locus_entropy

calculate_locus_entropy_ml <- function(mat) {
  rMaxs <- matrixStats::rowMaxs(mat)
  if(any(rMaxs==0)) #This site is not valid because a b value is not valid with the current model
    return(NA)
  freqs <- colSums((mat==rMaxs)+0)/nrow(mat)
  logFreqs <- log2(freqs)
  logFreqs[logFreqs == -Inf] <- 0 ## lim x -> 0 p log p = 0
  return(-sum(freqs*logFreqs))
}

#' Calculates the b-value entropy under they phyfum model
#'
#' @details The full/marginalized likelihood calculation does not work very well
#'   empirically and thus we recommend using the ML implementation.
#'
#' @param likelihoods 3D array with \code{[sample,peak,locus]} dimensions with
#'   the b-peak likelihoods for each sample and fCpG locus as returned by
#'   [get_peak_likelihoods()]
#' @param ml Logical that specifies whether the entropy is calculated using
#'   maximum likelihood peak assignments or marginalized over all possible peaks
#'
#' @returns expected entropy
#' @export

get_entropy <- function(likelihoods,ml=T) {
  if(ml==T){
    h <- mean(apply(likelihoods, MARGIN = 3, calculate_locus_entropy_ml,simplify = T),na.rm=T)
  } else {
    h <- mean(apply(likelihoods, MARGIN = 3, calculate_locus_entropy,simplify = T),na.rm=T)
  }
  return(h)
}

#' Compute full substitution saturation entropy
#'
#' This version only works well with small n and length(freqs)
#'
#' @details This version uses full enumeration of multinomial compositions and
#'   thus requires a lot of memory when freqs or n are large. The recursive
#'   version [get_fss_entropy_ml()] is usually preferred. Filtering compositions
#'   with very small probability speeds up the computation without big effect in
#'   the estimation.
#'
#' @param freqs Vector of category probabilities (i.e., frequencies, sums to 1)
#' @param n Total count (number of sequences)
#' @param prob_threshold Absolute multinomial probability cutoff (optional). Set
#'   to 0 to explore all compositions (Default).
#'
#' @returns list with the expected full saturation ML entropy estimator,
#'   expected full saturation ML entropy estimator without the
#'   skipped-probability normalization, expected full saturation skipped
#'   probability, probability threshold, and max_error
#' @seealso [[get_fss_entropy_ml()]]
#' @export

get_fss_entropy_ml_full_composition <- function(freqs, n, prob_threshold = 0) {
  if (!requireNamespace("partitions", quietly = TRUE)) {
    stop("Package 'partitions' is required. Please install it.")
  }

  counts <- t(partitions::compositions(n, length(freqs)))
  multinomDs <- apply(counts, 1, stats::dmultinom, prob = freqs)

  # Apply threshold to remove counts with very small probabilities
  countsToKeep <- multinomDs >= prob_threshold
  skipped_prob <- sum(multinomDs[!countsToKeep])

  counts <- counts[countsToKeep,]
  multinomDs <- multinomDs[countsToKeep]

  # Entropy calculation
  countFreqs <- counts / n
  logCountFreqs <- log2(countFreqs)
  logCountFreqs[logCountFreqs == -Inf] <- 0
  countEntropies <- rowSums(countFreqs * logCountFreqs) * -1

  entropy <- sum(multinomDs * countEntropies)
  norm_entropy <- if (skipped_prob >= 1) NA else entropy / (1 - skipped_prob)

  return(list(fssmle = norm_entropy, fssmleU = entropy, skipped_prob = skipped_prob, prob_threshold = prob_threshold))
}

#' Estimate a Multinomial Probability Threshold for Entropy Approximation
#'
#' This function estimates a probability threshold for skipping low-probability
#' compositions in multinomial-based entropy calculations. It ensures that the
#' total probability mass of all skipped compositions is less than or equal to
#' the specified error tolerance (`epsilon`).
#'
#' @inheritParams get_fss_entropy_ml_full_composition
#' @param epsilon Maximum total probability mass allowed to be skipped. Default
#'   is 0.05.
#' @param sample_size Number of compositions to sample for estimating the
#'   threshold. Default is 1M.
#'
#' @returns A numeric scalar giving the probability threshold. Compositions with
#'   multinomial probabilities lower than this value can be safely skipped while
#'   keeping total skipped probability mass ≤ `epsilon`.
#' @seealso [get_fss_entropy_ml()]
#' @keywords internal

estimate_prob_threshold <- function(freqs, n, epsilon = 0.05, sample_size = 1000000) {
  k <- length(freqs)
  samples <- stats::rmultinom(sample_size, size = n, prob = freqs)
  probs <- apply(samples, 2, stats::dmultinom, prob = freqs)
  sorted_probs <- sort(probs)
  cum_probs <- cumsum(sorted_probs)
  idx <- which(cum_probs <= epsilon)
  if (length(idx) == 0) {
    message("All sampled probabilities exceed epsilon; no compositions will be skipped.")
    return(0)
  }
  return(sorted_probs[max(idx)])
}

#' Compute full substitution saturation entropy using recursion
#'
#' This version only works well with small n and length(freqs)
#'
#' @details This version uses recursion with to calculate multinomial
#'   compositions on the fly and thus require much less memory than
#'   [get_fss_entropy_ml_full_composition()]. Filtering compositions with very
#'   small probability speeds up the computation without big effect in the
#'   estimation.
#'
#' @inherit get_fss_entropy_ml_full_composition params return
#' @seealso [[get_fss_entropy_ml_full_composition()]]
#' @export

get_fss_entropy_ml <- function(freqs, n, prob_threshold = NULL) {
  nFreqs <- length(freqs)

  total <- 0
  skipped_prob <- 0
  generate <- function(current, left, i) {
    if (i == nFreqs) {
      current[i] <- left
      prob <- stats::dmultinom(current, prob = freqs)
      if (prob < prob_threshold) {
        skipped_prob <<- skipped_prob + prob # superassignment operator
        return()
      }
      countFreqs <- current / n
      logCountFreqs <- log2(countFreqs)
      logCountFreqs[!is.finite(logCountFreqs)] <- 0
      entropy <- -sum(countFreqs * logCountFreqs)
      total <<- total + prob * entropy # superassignment operator
    } else {
      for (val in 0:left) {
        current[i] <- val
        generate(current, left - val, i + 1)
      }
    }
  }

  generate(integer(nFreqs), n, 1)
  norm_entropy <- if (skipped_prob >= 1) NA else total / (1 - skipped_prob)

  return(list(fssmle = norm_entropy, fssmleU = total, skipped_prob = skipped_prob, prob_threshold = prob_threshold))
}


#LEGACY CODE
############

#' Simulates b-values
#'
#' Needed for the full-probability fss estimation
#'
#' This is old code not currently used and has memory and speed issues unless S
#' is very due small since it pre-computes all the counts and does not discard
#' cases with very small multinomial probability
#'
#' @param counts partitions::compositions(n,length(freqs))
#' @param delta flipflop.delta parameter of the phyfum model (methylated
#'   background noise)
#' @param eta flipflop.eta parameter of the phyfum model (de-methylated
#'   background negative noise)
#' @param kappa flipflop.kappa parameter of the phyfum model (b concentration)
#'
#' @returns b-value alignment (i.e., matrix (or 2-dimensional array) with a
#'   column per fCpG site and row per sample) like the one needed by
#'   [get_peak_likelihoods()]
#' @seealso [get_fss_entropy()]
#' @keywords internal, unsupported

simulate_betas <- function(counts, delta, eta, kappa) {
  peak_count <- ncol(counts)
  nSeqs <- unique(rowSums(counts))
  nSim <- unique(colSums(counts))

  betas <- (eta - delta) * (seq_len(peak_count) - 1) / (peak_count - 1) + delta
  simulatedBetas <- matrix(stats::rbeta(nSim * peak_count, betas * kappa, (1 - betas) * kappa), nrow = peak_count, byrow = F)

  lastPeakBetaValueI <- rep(1,peak_count)

  finalData <- array(NA,c(nSeqs,nrow(counts)))
  tempSeq <- array(NA,nSeqs)
  iTempSeq <- NA
  nBetasToCopyThisPeak <- NA

  for (iCount in 1:nrow(counts)) {
    iTempSeq <- 1
    for (iPeak in 1:peak_count){
      nBetasToCopyThisPeak <- counts[iCount,iPeak]
      if(nBetasToCopyThisPeak == 0)
        next
      tempSeq[iTempSeq:(iTempSeq-1+nBetasToCopyThisPeak)] <- simulatedBetas[iPeak,lastPeakBetaValueI[iPeak]:(lastPeakBetaValueI[iPeak]+nBetasToCopyThisPeak-1)]
      iTempSeq <- iTempSeq + nBetasToCopyThisPeak
      lastPeakBetaValueI[iPeak] <- lastPeakBetaValueI[iPeak] + nBetasToCopyThisPeak
      if(iTempSeq == peak_count)
        break
    }
    finalData[,iCount] <- tempSeq
  }

  return(finalData)
}

#' Compute full substitution saturation entropy
#'
#' Old version that should not be used currently. Use [get_fss_entropy_ml()]
#' instead.
#'
#' @details This is an old version that allows to use ML or full-probability.
#'   The full-probability only uses compositions and simulation and does not
#'   include probability filtering
#'
#' @inherit get_fss_entropy_ml params return
#' @inheritParams simulate_betas
#' @inheritDotParams get_fss_entropy_ml prob_threshold
#' @inheritParams get_entropy
#'
#' @seealso [get_fss_entropy_ml()]
#' @keywords internal, unsupported


get_fss_entropy <- function(freqs,n,delta,eta,kappa,ml=F,...){
  if(ml==T)
    return(get_fss_entropy_ml(freqs,n,...))
  counts <- t(partitions::compositions(n,length(freqs)))
  betas <- simulate_betas(counts,delta,eta,kappa)
  #partials <- get_peak_likelihoods(betas,delta,eta,kappa,length(freqs)) #TODO check this. get_peak_likelihoods expects s not peak_count!
  partials <- get_peak_likelihoods(betas,delta,eta,kappa,(length(freqs)-1)/2)
  hs <- apply(partials, MARGIN = 3, calculate_locus_entropy,simplify = T)
  multinomDs <- apply(counts,1,stats::dmultinom,prob=freqs)
  return(sum(hs*multinomDs,na.rm=T))
}

