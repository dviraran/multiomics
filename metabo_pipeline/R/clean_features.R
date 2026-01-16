#' Clean feature table: filter, impute, outlier detection
#' @param feature_table Raw feature table
#' @return Cleaned feature table
clean_features <- function(feature_table) {
  # Example: filter by prevalence, impute missing, PCA outlier removal
  library(dplyr)
  keep <- colMeans(!is.na(feature_table)) > 0.8
  ft <- feature_table[, keep]
  ft[is.na(ft)] <- apply(ft, 2, function(x) mean(x, na.rm = TRUE))
  # Outlier detection (PCA-based)
  pca <- prcomp(ft, scale. = TRUE)
  dists <- sqrt(rowSums(pca$x^2))
  outlier <- dists > (mean(dists) + 3 * sd(dists))
  ft <- ft[!outlier, ]
  return(ft)
}
