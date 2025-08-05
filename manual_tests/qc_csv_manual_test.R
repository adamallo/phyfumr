library(data.table)
sample_data <- fread("~/projects/flipFlop/newG22new/Gabbutt2022-sampleinfo-corrected.csv")
sample_data <- sample_data[!patient%in%c(NA,"US"),]
dirs <- "~/projects/flipFlop/latestG22/tissueUnspecificSelection/nf/data/"
files <- list.files(dirs,pattern = "*.csv",full.names = T)
results <- qc_csv(files,sample_data,cl = 8)
delete_cached_CSV_data()
