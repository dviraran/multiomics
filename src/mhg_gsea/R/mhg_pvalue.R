#' Calculate Exact mHG P-value
#'
#' Computes the exact p-value for a given mHG score using the dynamic programming algorithm.
#'
#' @param score The observed mHG statistic (minimum p-value).
#' @param N Total number of genes.
#' @param B Size of the gene set.
#' @return The exact p-value corrected for multiple cutoffs.
#' @export
mhg_pvalue <- function(score, N, B) {
    if (score >= 1) {
        return(1.0)
    }
    if (B == 0) {
        return(1.0)
    }

    # Calculate thresholds (limits) for the DP
    # We want to identify the "safe" region where HGT p-value > score.
    # For each cutoff n (1 to N), we need to find max hits k such that P(X >= k) > score.
    # P(X >= k) = P(X > k-1).
    # We use qhyper with lower.tail = FALSE:
    # qhyper(p, m, n, k, lower.tail = FALSE) returns the smallest x such that P(X > x) <= p.
    # Let x_crit = qhyper(score, ...).
    # Then P(X > x_crit) <= score. This implies that having x_crit+1 hits is significant (p <= score).
    # So we must have hits <= x_crit to remain non-significant.
    # Thus, the limit is exactly x_crit.

    n_vec <- 1:N
    limits <- qhyper(score, B, N - B, n_vec, lower.tail = FALSE)

    # Call C++ function
    # Ensure integer type for limits
    limits <- as.integer(limits)

    pval <- mhg_exact_pvalue(limits, N, B)

    return(pval)
}
