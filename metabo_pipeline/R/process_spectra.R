#' Process spectra: peak detection, alignment, quantification
#' @param mSet MetaboAnalystR object
#' @return Feature table (data.frame)
process_spectra <- function(mSet) {
  mSet <- PerformPeakAnnotation(mSet)
  mSet <- PerformPeakProfiling(mSet)
  feature_table <- mSet$dataSet$peakint
  return(feature_table)
}
