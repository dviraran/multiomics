#' Run univariate and multivariate statistics
#' @param norm_data Normalized data
#' @return List of results (fold changes, p-values, models)
run_stats <- function(norm_data) {
  library(MetaboAnalystR)
  mSet <- InitDataObjects("pktable", "stat", FALSE)
  mSet$dataSet$norm <- norm_data
  mSet <- Ttests.Anal(mSet, nonpar = FALSE, threshp = 0.05, paired = FALSE)
  mSet <- PCA.Anal(mSet)
  mSet <- PLSDA.Anal(mSet)
  results <- list(
    ttest = mSet$analSet$tt,
    pca = mSet$analSet$pca,
    plsda = mSet$analSet$plsda
  )
  return(results)
}
