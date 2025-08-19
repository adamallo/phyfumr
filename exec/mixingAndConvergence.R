#!/usr/bin/env Rscript

library(phyfumr)

#Parsing command-line arguments and checking input
#######################################
args <- commandArgs(trailingOnly = TRUE) #Comment for debugging, uncomment to run with Rscript

if (length(args) <=3) {
  stop("Usage: script burnin_proportion n_cores baseOutdir sample1.tree ... sampleN.tree")
} else {
  ##ForPablo: Remove if we don't want this to be verbose
  print(paste0("Running the script in the trace files (",paste(collapse =", ",args[-c(1:3)]),") with a burnin of ",args[1]," proportion of the trees and ",args[2]," processors"))
}

burnin_p <- as.numeric(args[1]) #ForPablo: the default should be 0.1
out_dir <- args[3]
files <- args[-c(1:3)]


if(!all(file.exists(files))){
  stop("ERROR: Not all input files exist. ",paste(files,collapse=", "))
}
#######################################

#Preparing output
#######################################
run_name <- gsub(x = basename(files[1]),pattern = ".trees$",replacement = "")
print(paste0("Saving outputs in ",out_dir)) #ForPablo: remove?
#######################################

capture.output(mcmc_qc_patient(files,
                               out_dir = out_dir,
                               plot_dir = out_dir,
                               n_cores = 1,
                               backup = T
                               ),
               file = nullfile())
