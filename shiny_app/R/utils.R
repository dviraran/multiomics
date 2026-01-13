# =============================================================================
# Utility Functions for Multi-Omics Viewer
# =============================================================================

#' Format p-value for display
#' @param p P-value
#' @param digits Number of significant digits
#' @return Formatted string
format_pvalue <- function(p, digits = 3) {
    if (is.na(p)) return("NA")
    if (p < 0.001) {
        return(format(p, digits = 2, scientific = TRUE))
    }
    return(format(round(p, digits), nsmall = digits))
}

#' Format numbers with commas
#' @param x Numeric value
#' @return Formatted string
format_number <- function(x) {
    if (is.na(x)) return("NA")
    format(x, big.mark = ",", scientific = FALSE)
}

#' Extract variance explained from PCA results
#' @param pca_data Data frame with PCA results
#' @return Named numeric vector of variance explained
extract_variance_explained <- function(pca_data) {
    # Look for variance columns
    var_cols <- grep("^var|variance", colnames(pca_data), value = TRUE, ignore.case = TRUE)

    if (length(var_cols) > 0) {
        # If variance stored in columns
        return(as.numeric(pca_data[1, var_cols]))
    }

    # Otherwise, try to extract from column names (e.g., "PC1_12.5%")
    pc_cols <- grep("^PC[0-9]", colnames(pca_data), value = TRUE)
    var_exp <- sapply(pc_cols, function(col) {
        match <- regmatches(col, regexpr("[0-9.]+%", col))
        if (length(match) > 0) {
            return(as.numeric(gsub("%", "", match)))
        }
        NA
    })

    if (all(is.na(var_exp))) {
        # Return dummy values
        return(rep(NA, length(pc_cols)))
    }

    var_exp
}

#' Get metadata columns suitable for coloring/grouping
#' @param metadata Data frame with sample metadata
#' @param exclude Columns to exclude
#' @return Character vector of column names
get_metadata_columns <- function(metadata, exclude = c("sample", "sample_id")) {
    if (is.null(metadata)) return(character(0))

    cols <- colnames(metadata)
    cols <- setdiff(cols, exclude)

    # Filter to reasonable columns (not too many unique values)
    good_cols <- sapply(cols, function(col) {
        n_unique <- length(unique(metadata[[col]]))
        n_unique >= 2 && n_unique <= 50 && !all(is.na(metadata[[col]]))
    })

    cols[good_cols]
}

#' Prepare expression data for boxplot (wide to long format)
#' @param expr_matrix Expression matrix (features x samples)
#' @param feature_id Feature to extract
#' @param metadata Sample metadata
#' @param sample_col Column in metadata with sample IDs
#' @return Long-format data frame
prepare_expression_data <- function(expr_matrix, feature_id, metadata,
                                    sample_col = "sample") {

    # Find the feature
    if (feature_id %in% rownames(expr_matrix)) {
        expr_values <- as.numeric(expr_matrix[feature_id, ])
    } else if (feature_id %in% expr_matrix[[1]]) {
        # First column might be feature IDs
        idx <- which(expr_matrix[[1]] == feature_id)
        expr_values <- as.numeric(expr_matrix[idx, -1])
    } else {
        return(NULL)
    }

    # Get sample names
    sample_names <- colnames(expr_matrix)
    if (sample_names[1] %in% c("feature_id", "gene_id", "protein_id")) {
        sample_names <- sample_names[-1]
    }

    # Create data frame
    df <- data.frame(
        sample = sample_names,
        expression = expr_values,
        stringsAsFactors = FALSE
    )

    # Merge with metadata if available
    if (!is.null(metadata)) {
        df <- merge(df, metadata, by.x = "sample", by.y = sample_col, all.x = TRUE)
    }

    df
}

#' Create summary statistics for a data table
#' @param de_results DE results data frame
#' @param pval_col P-value column name
#' @param log2fc_col Log2FC column name
#' @param pval_thresh P-value threshold
#' @param log2fc_thresh Log2FC threshold
#' @return List with summary statistics
summarize_de_results <- function(de_results, pval_col = "padj",
                                 log2fc_col = "log2FoldChange",
                                 pval_thresh = 0.05, log2fc_thresh = 1) {

    n_total <- nrow(de_results)
    sig_mask <- de_results[[pval_col]] < pval_thresh
    up_mask <- sig_mask & de_results[[log2fc_col]] > log2fc_thresh
    down_mask <- sig_mask & de_results[[log2fc_col]] < -log2fc_thresh

    list(
        n_total = n_total,
        n_significant = sum(sig_mask, na.rm = TRUE),
        n_up = sum(up_mask, na.rm = TRUE),
        n_down = sum(down_mask, na.rm = TRUE),
        pct_significant = round(100 * sum(sig_mask, na.rm = TRUE) / n_total, 1)
    )
}

#' Get available contrasts from DE results
#' @param de_results List of DE result data frames
#' @return Character vector of contrast names
get_contrasts <- function(de_results) {
    if (is.null(de_results)) return(character(0))
    if (is.data.frame(de_results)) {
        if ("contrast" %in% colnames(de_results)) {
            return(unique(de_results$contrast))
        }
        return("default")
    }
    if (is.list(de_results)) {
        return(names(de_results))
    }
    character(0)
}

#' Filter DE results by contrast
#' @param de_results List or data frame of DE results
#' @param contrast Contrast name to filter
#' @return Data frame for the specified contrast
get_de_for_contrast <- function(de_results, contrast) {
    if (is.null(de_results)) return(NULL)

    if (is.data.frame(de_results)) {
        if ("contrast" %in% colnames(de_results)) {
            return(de_results[de_results$contrast == contrast, ])
        }
        return(de_results)
    }

    if (is.list(de_results) && contrast %in% names(de_results)) {
        return(de_results[[contrast]])
    }

    NULL
}

#' Create info box content
#' @param title Box title
#' @param value Main value
#' @param subtitle Optional subtitle
#' @param color Box color
#' @return HTML for info box
create_info_box <- function(title, value, subtitle = NULL, color = "#3498db") {
    subtitle_html <- ""
    if (!is.null(subtitle)) {
        subtitle_html <- paste0('<div style="font-size: 12px; color: #666;">',
                                subtitle, '</div>')
    }

    HTML(paste0(
        '<div style="background: ', color, '; color: white; padding: 15px; ',
        'border-radius: 5px; text-align: center; margin: 5px;">',
        '<div style="font-size: 14px; font-weight: bold;">', title, '</div>',
        '<div style="font-size: 24px; font-weight: bold;">', value, '</div>',
        subtitle_html,
        '</div>'
    ))
}

#' Truncate text to max length
#' @param text Text to truncate
#' @param max_length Maximum length
#' @return Truncated text with "..." if needed
truncate_text <- function(text, max_length = 50) {
    if (nchar(text) > max_length) {
        return(paste0(substr(text, 1, max_length - 3), "..."))
    }
    text
}
