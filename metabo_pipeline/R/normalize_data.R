#' Normalize and transform data
#' @param cleaned_data Cleaned feature table
#' @return Normalized data matrix
normalize_data <- function(cleaned_data) {
  library(MetaboAnalystR)
  mSet <- InitDataObjects("pktable", "stat", FALSE)
  mSet$dataSet$peakint <- cleaned_data
  mSet <- PerformNormalization(mSet, method = "SumNorm", trans = "Log", scale = "Auto")
  norm_data <- mSet$dataSet$norm
  return(norm_data)
}
