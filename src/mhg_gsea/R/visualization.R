#' Plot Running Enrichment
#'
#' Generates a Running Enrichment Plot showing the HyperGeometric Tail (HGT)
#' profile across ranks and marking the optimal cutoff.
#'
#' @param stat_res The result list from the `mhg_stat` function.
#' @param title The title of the plot (optional).
#' @return A ggplot object.
#' @import ggplot2
#' @export
plot_running_enrichment <- function(stat_res, title = "Running Enrichment") {
    if (is.null(stat_res$pvals)) stop("No p-values in stat_res object. Ensure mhg_stat was run.")

    # Create data frame for plotting
    df <- data.frame(
        Rank = 1:length(stat_res$pvals),
        Score = -log10(stat_res$pvals)
    )

    # Optimal cutoff
    optimal_rank <- stat_res$n
    optimal_score <- -log10(stat_res$mHG)

    p <- ggplot(df, aes(x = Rank, y = Score)) +
        geom_line(color = "darkblue", linewidth = 1) +
        geom_vline(xintercept = optimal_rank, linetype = "dashed", color = "red") +
        geom_point(
            data = data.frame(Rank = optimal_rank, Score = optimal_score),
            aes(x = Rank, y = Score), color = "red", size = 3
        ) +
        labs(
            title = title,
            x = "Rank",
            y = "-log10(HyperGeometric Tail p-value)"
        ) +
        theme_minimal() +
        annotate("text",
            x = optimal_rank, y = optimal_score,
            label = paste0("  Cutoff: ", optimal_rank, "\n  b: ", stat_res$b),
            hjust = 0, vjust = 1, color = "red"
        )

    return(p)
}
