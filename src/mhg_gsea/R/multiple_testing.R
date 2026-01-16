#' Adjust P-values for Multiple Comparisons
#'
#' Applies Benjamini-Hochberg FDR correction to a vector of p-values.
#'
#' @param pvals A numeric vector of p-values.
#' @param method Correction method (default "BH").
#' @return A numeric vector of adjusted p-values.
#' @export
adjust_fdr <- function(pvals, method = "BH") {
    return(p.adjust(pvals, method = method))
}
