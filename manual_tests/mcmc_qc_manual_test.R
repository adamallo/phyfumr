library(phyfumr)
baseDir <- "~/Downloads/short_worflow_res/test"
inDir <- "~/Documents/Ciencia/Postdoc/projects/flipFlop/beFirstRun/"
trees_files <- c(paste(sep="/",inDir,"666/169flipflopS8.trees"),paste(sep="/",inDir,"999/169flipflopS8.trees"))
out_dir <- "~/test/out_dir"
plot_dir <- "~/test/plot_dir"
convergence <- mcmc_qc_patient(trees_files,out_dir = out_dir,plot_dir = plot_dir)
