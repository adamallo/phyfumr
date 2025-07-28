# R/data_parsing.R

#' Reads b-values from a Phyfum XML file
#'
#' @param xml_file xml file
#' @param matrix whether the return should be a data table (column by sample) or an alignment-like matrix (column by fCpG). All phyfumr functions expect matrix=T (Default)
#'
#' @returns b-value dataset in data.table or matrix form
#' @export

get_Bdata_from_XML <- function(xml_file,matrix=T){
  if(!file.exists(xml_file)) {
    warning(paste0(xml_file," does not exits"))
    return(NULL)
  }

  xmlContent <- XML::xmlToList(XML::xmlParse(xml_file))

  bData <- data.table::rbindlist(lapply(xmlContent$afalignment,function(thisAF){
    if(is.list(thisAF) && !is.null(thisAF$text)){
      data.table::data.table(b=as.numeric(unlist(strsplit(split = ",",x = gsub(pattern = "\n *",replacement = "",thisAF$text)))),sample=thisAF$taxon)}
  }))

  if(matrix==T) {
    bData[, locus := (seq_len(.N)),by=sample]
    bData <- t(as.matrix(data.table::dcast(bData,locus~sample,value.var = "b")[,`:=`(locus=NULL)]))
  }

  return(bData)
}

#' Caches b-values stored in CSV files
#'
#' Do not use this function in parallel.
#'
#' @details CSV files may contain data from multiple patients but the CSV
#'   parsing is done by patient [get_Bdata_from_CSVs()]. To avoid re-loading the
#'   same files multiple times, this function can be used first to load all
#'   files in memory and then parse them efficiently. This function must not be
#'   used in parallel since it stores the data in a package environment. All
#'   downstream analyses may be conducted in parallel.
#'
#' @param csv_dirs list of directories for which to read all CSV files
#'
#' @returns NULL
#' @seealso [get_Bdata_from_CSVs()] [delete_cached_CSV_data()]
#' @export

load_Bdata_from_CSVs <- function(csv_dirs){
  lapply(csv_dirs,function(csv_dir){
    if(is.null(.phyfumr_env$loadedCSVs[[csv_dir]])){
      theseCSVs <- list.files(csv_dir,pattern = ".csv$",full.names = T)
      theseDTs <- lapply(theseCSVs,function(x){data.table::setnames(data.table::fread(x),1,"CpG")})
      .phyfumr_env$loadedCSVs[[csv_dir]] <- Reduce(function(x,y) merge(x,y,by="CpG",all=T),theseDTs)
    } else {
      warning(paste0("Reusing previously-loaded CSV data from ",csv_dir))
    }})
}

#' Deletes CSV data cached by [load_Bdata_from_CSVs]
#'
#' @returns NULL
#' @seealso [get_Bdata_from_CSVs()] [load_Bdata_from_CSVs()]
#' @export

delete_cached_CSV_data <- function(){
  rm(list=c(loadedCSVs),envir=.phyfumr_env)
}

#' Obtains b-values from CSV files (cached or not)
#'
#' @details CSV files may contain data from multiple patients but the CSV
#'   parsing is done by patient. To avoid re-loading the same files multiple
#'   times, [load_Bdata_from_CSVs()] can be used first to load all files in memory
#'   and then parse them efficiently here. This function can be used in parallel
#'   if the data has been previously cached. This function assumes it is
#'   executed in parallel and stops when it is not safe to do so. If the user is
#'   certain the function is not run in parallel, isParallel should be set to F.
#'   The cached data can be deleted using [delete_cached_CSV_data()]
#'
#' @param csv_dir directory with CSV files to parse
#' @param patient_id patient ID to obtain the CSV column names to parse from
#'   sample_pairing
#' @param sample_pairing data.table with a patient column and a samplename
#'   column.
#' @inherit get_Bdata_from_XML params return
#' @param isParallel True by default, assuming this function is executed in
#'   parallel and warning when this is not safe to do.
#'
#' @seealso [load_Bdata_from_CSVs()] [delete_cached_CSV_data]
#' @export

get_Bdata_from_CSVs <- function(csv_dir,patient_id,sample_pairing,matrix=T,isParallel=T){
  if(is.null(.phyfumr_env$loadedCSVs[[csv_dir]]) && isParallel==F){
    load_Bdata_from_CSVs(csv_dir)
  } else if (isParallel==T) {
    stop(paste0("This function should not be run in parallel",
                " unless the data has been pre-loaded with load_Bdata_from_CSVs.",
                " If you are not executing this function in parallel, ",
                "feel free to use the isParallel=F argument and execute it again"))
  }

  if(is.null(.phyfumr_env$loadedCSVs[[csv_dir]])){
    warning(paste0("CSV files not found in ", csv_dir))
    return(NULL)
  }
  if(length(sample_pairing[patient==patient_id,samplename])==0){
    warning(paste0("No patient information for patient_id ",patient_id))
    return(NULL)
  }

  theseData <- .phyfumr_env$loadedCSVs[[csv_dir]][,.SD,.SDcols=sample_pairing[patient==patient_id,samplename]]
  cleanData <- theseData[rowSums(!is.na(theseData)) > 0,]
  if(matrix==T) return(t(as.matrix(cleanData)))
  return(data.table::melt(cleanData,variable.name = "sample",value.name = "b")[,`:=`(sample=as.character(sample))])
}


#' Parses a phyfum .log file and summarizes a list of parameters
#'
#' @param log_file name of the log file to parse
#' @param burnin proportion of samples to be removed as burnin
#' @param sumfun summary function. One of mean or median.
#' @param params list of parameters to parse
#'
#' @returns data.table with a a column per param with its point estimate
#' @export

estimate_params <- function(log_file,
                      burnin=0.1,
                      sumfun=c("mean","median"),
                      params=c("errorModel.deltaOffset",
                               "errorModel.etaOffset",
                               "errorModel.kappaScale",
                               "flipflop.mu",
                               "flipflop.gamma",
                               "flipflop.lambda",
                               "alignment.stemCells")) {

  sumfun <- match.arg(sumfun)

  if(!file.exists(log_file))
    stop(paste0(log_file," does not exits"))

  b_data <- data.table::fread(log_file)
  b_data[, valid := .I/.N>burnin]

  return(as.list(b_data[valid==T,lapply(.SD,sumfun),.SDcols=params]))
}

