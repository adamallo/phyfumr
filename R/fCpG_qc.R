# R/fCpG_qc.R

# Functions to carry out the QC of a fCpG dataset

#' QC of selected fCpG sites
#'
#' @details This function can runs in parallel (per patient) if pbapply is
#'   available
#'
#' @param csvs csv file(s) containing b-values. A row per CpG site and column
#'   per sample, with row and column names.
#' @param patientData data.table with columns for samplename and patient. The
#'   samplename columns should correspond with the column/sample name in the CSV
#'   or XML files
#' @param max_rel_emd_to_origin maximum mean relative earth's mover distance to
#'   the origin
#' @param max_rel_entropy maximum mean relative entropy
#' @param model_params user_defined model parameters used for the computation.
#'   Obtain the defaults using [get_model_params()] and modify them to your
#'   needs.
#' @param min_multinomial_p minimum multinomial probability to be included in
#'   the calculation of the full saturation entropy
#' @param ml_entropies whether the entropy calculation should be carried out
#'   using maximum likelihood state asignment or not. See [get_entropy()]
#' @param cl cl settings for [pbapply::pblapply()] if available
#' @param ... arguments for estimate_prob_threshold
#'
#' @returns data.table with distance to the origin and relative entropy data per
#'   patient, including if they are suspicious of being problematic (outside the
#'   max* thresholds)
#' @export
#'

#TODO optimize thresholds
qc_csv <- function(csvs,
                   patientData,
                  max_rel_emd_to_origin=0.5,
                  max_rel_entropy=0.95,
                  model_params=NULL,
                  min_multinomial_p=0.95,
                  ml_entropies=T,
                  cl=1,
                  ...){
  if(ml_entropies==F)
    stop("ERROR: option not implemented yet. Use ml_entropies=T")

  if(is.null(model_params)){
    model_params <- get_model_params()
  }

  load_Bdata_from_CSVs(csvs)
  patients <- patientData[,unique(patient)]

  infun <- function(patient) {
    theData <- get_Bdata_from_CSVs(csvs,
                                   patient,
                                   patientData)
    if(is.null(theData))
      stop(paste("Missing data for patient ",patient))

    n_samples <- nrow(theData)
    peak_likelihoods <- get_peak_likelihoods(theData,
                                             delta = model_params$errorModel.deltaOffset,
                                             eta = model_params$errorModel.etaOffset,
                                             kappa = model_params$errorModel.kappaScale,
                                             s = model_params$alignment.stemCells)
    rm(theData)
    hml <- get_entropy(likelihoods = peak_likelihoods,
                       ml = ml_entropies)
    e_freqs <- get_empirical_freqs(likelihoods = peak_likelihoods,
                                   ml = ml_entropies)
    rm(peak_likelihoods)
    fss <- get_fss_entropy_ml(freqs = e_freqs,
                              n = n_samples,
                              prob_threshold = estimate_prob_threshold(freqs = e_freqs,
                                                                       n = n_samples,
                                                                       epsilon = min_multinomial_p,
                                                                       ...))
    distance_to_origin <- emd_to_origin(freqs = e_freqs)

    return(data.table::data.table(patient = patient,
                                  Iml = hml/fss$fssmle,
                                  remdO = distance_to_origin$remdO,
                                  emdO = distance_to_origin$emdO,
                                  flag = (hml/fss$fssmle > max_rel_entropy) + (distance_to_origin$remdO > max_rel_emd_to_origin)))

  }

  if(requireNamespace("pbapply", quietly = TRUE)){
    results <- data.table::rbindlist(pbapply::pblapply(patients, infun, cl = cl))
  } else {
    results <- data.table::rbindlist(lapply(patients, infun))
  }
  return(results)
}
