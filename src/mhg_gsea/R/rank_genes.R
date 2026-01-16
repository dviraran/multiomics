#' Rank Genes for GSEA
#'
#' Ranks genes based on a signed log-p value metric: sign(log2FC) * (-log10(pvalue)).
#'
#' @param df A data frame containing differential expression results.
#' @param gene_col Character string specifying the column name for gene symbols (optional).
#' If NULL, rownames are used.
#' @param log2fc_col Character string specifying the column name for log2 fold changes.
#' @param pval_col Character string specifying the column name for p-values.
#' @return A named numeric vector of gene scores, sorted in descending order.
#' @export
rank_genes <- function(df, gene_col = NULL, log2fc_col = "log2FoldChange", pval_col = "pvalue") {
    if (!all(c(log2fc_col, pval_col) %in% colnames(df))) {
        stop(paste("Columns", log2fc_col, "and/or", pval_col, "not found in data frame."))
    }

    if (!is.null(gene_col)) {
        if (!(gene_col %in% colnames(df))) stop(paste("Gene column", gene_col, "not found."))
        gene_names <- df[[gene_col]]
    } else {
        gene_names <- rownames(df)
    }

    # Handle zero p-values
    pvals <- df[[pval_col]]
    min_nonzero <- min(pvals[pvals > 0], na.rm = TRUE)
    if (is.infinite(min_nonzero)) min_nonzero <- .Machine$double.xmin

    # Replace 0 with a value slightly smaller than min_nonzero
    pvals[pvals == 0] <- min_nonzero * 1e-1

    # Calculate metric
    scores <- sign(df[[log2fc_col]]) * (-log10(pvals))

    names(scores) <- gene_names

    # Return sorted
    return(sort(scores, decreasing = TRUE))
}
