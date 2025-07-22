# data_generation/make_prior_default_settings.R
devtools::load_all()
hash <- digest::digest(file = "R/prior_default_settings.R")
hash_file <- "data_generation/prior_default_settings.hash"

if (!file.exists(hash_file) || readLines(hash_file) != hash) {
  source("R/prior_default_settings.R")
  usethis::use_data(default_prior_settings, overwrite = TRUE)
  writeLines(hash, hash_file)
}
