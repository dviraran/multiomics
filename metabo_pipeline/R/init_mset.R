#' Initialize MetaboAnalystR mSet object for LC-MS
#' @param raw_files List of raw data file paths
#' @return mSet object
init_mset <- function(raw_files) {
  library(MetaboAnalystR)
  mSet <- InitDataObjects("pktable", "mspeak", FALSE)
  mSet <- Read.MSData(mSet, files = raw_files, format = "mzML", mode = "positive")
  return(mSet)
}
