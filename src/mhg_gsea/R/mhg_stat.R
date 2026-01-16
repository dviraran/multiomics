#' Calculate minimum HyperGeometric (mHG) Statistic
#'
#' @param rank_vec Named numeric vector of gene scores, sorted descending.
#' @param gene_set Character vector of gene symbols in the set.
#' @param min_size Minimum intersection size to consider (default 5).
#' @return A list containing the mHG statistic, optimal cutoff index (n),
#' number of hits at that cutoff (b), total set size (B), and all p-values.
#' @export
mhg_stat <- function(rank_vec, gene_set, min_size = 5) {
    genes <- names(rank_vec)
    N <- length(rank_vec)

    # Logical vector: is gene in set?
    hits <- genes %in% gene_set
    B <- sum(hits)

    if (B < min_size) {
        return(list(mHG = 1.0, b = 0, n = 0, N = N, B = B, threshold = NULL, pvals = rep(1, N)))
    }

    # k[i] is the number of hits in top i genes
    k <- cumsum(hits)
    n <- 1:N

    # phyper(x, m, n, k, lower.tail=FALSE) calculates P(X > x)
    # We want P(X >= k_observed) = P(X > k_observed - 1)
    # m = B (size of gene set, white balls)
    # n = N - B (size of background, black balls)
    # k = n (sample size, balls drawn)

    pvals <- phyper(k - 1, B, N - B, n, lower.tail = FALSE)

    # Find minimum p-value
    min_p <- min(pvals)

    # Find index of minimum (use first occurrence for smallest cutoff)
    idx <- which.min(pvals)

    return(list(
        mHG = min_p,
        b = k[idx],
        n = n[idx],
        N = N,
        B = B,
        pvals = pvals,
        hits = hits
    ))
}
