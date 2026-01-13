# =============================================================================
# Pathway and Set Enrichment Analysis
# =============================================================================
#
# This file uses shared utilities from ../shared/R/ when available for:
# - GMT file reading and validation
# - Support for custom pathway definitions

# Source shared utilities if available
.source_shared_utils <- function() {
    possible_paths <- c(
        file.path(dirname(dirname(getwd())), "shared", "R"),
        file.path(dirname(getwd()), "shared", "R"),
        file.path(getwd(), "..", "shared", "R"),
        file.path(getwd(), "..", "..", "shared", "R")
    )

    for (shared_dir in possible_paths) {
        if (dir.exists(shared_dir)) {
            gmt_file <- file.path(shared_dir, "gmt_utils.R")
            if (file.exists(gmt_file)) {
                source(gmt_file)
                message("Loaded shared GMT utilities")
                return(TRUE)
            }
        }
    }
    return(FALSE)
}

.source_shared_utils()

#' Run enrichment analysis
#'
#' @param de_results Differential analysis results
#' @param imputed_data Imputed data
#' @param ingested_data Ingested data (for annotations)
#' @param config Configuration list
#' @return Enrichment results
run_enrichment_analysis <- function(de_results, imputed_data, ingested_data, config) {
  log_message("=== Starting Enrichment Analysis ===")

  if (!config$enrichment$run_enrichment %||% TRUE) {
    log_message("Enrichment analysis disabled. Skipping.")
    return(NULL)
  }

  if (is.null(de_results)) {
    log_message("No differential analysis results. Skipping enrichment.")
    return(NULL)
  }

  annotations <- ingested_data$annotations
  all_results <- list()

  # Check what enrichment resources are available
  has_pathway_mapping <- !is.null(annotations$pathway_mapping)
  has_annotation_table <- !is.null(annotations$annotation_table)
  has_gmt <- !is.null(annotations$metabolite_sets)
  has_class_annotation <- has_annotation_table && "class" %in% colnames(annotations$annotation_table)

  if (!has_pathway_mapping && !has_gmt && !has_class_annotation) {
    log_message("No pathway mapping, GMT file, or class annotations available.")
    log_message("Skipping enrichment analysis. Provide pathway_mapping, gmt_file, or class annotations in config.")
    return(NULL)
  }

  for (contrast_name in names(de_results$results)) {
    log_message("Running enrichment for: ", contrast_name)

    de_table <- de_results$results[[contrast_name]]$table
    contrast_results <- list()

    # Pathway enrichment (if pathway mapping available)
    if (has_pathway_mapping) {
      pathway_results <- run_pathway_enrichment(de_table, annotations$pathway_mapping, config)
      if (!is.null(pathway_results)) {
        contrast_results$pathway <- pathway_results
      }
    }

    # GMT-based enrichment
    if (has_gmt) {
      gmt_results <- run_gmt_enrichment(de_table, annotations$metabolite_sets, config)
      if (!is.null(gmt_results)) {
        contrast_results$gmt <- gmt_results
      }
    }

    # Chemical class enrichment
    if (has_class_annotation && config$enrichment$run_class_enrichment %||% TRUE) {
      class_results <- run_class_enrichment(de_table, annotations$annotation_table, config)
      if (!is.null(class_results)) {
        contrast_results$class <- class_results
      }
    }

    all_results[[contrast_name]] <- contrast_results

    # Save results
    clean_name <- gsub(" ", "_", gsub("-", "_vs_", contrast_name))
    if (!is.null(contrast_results$pathway)) {
      save_table(contrast_results$pathway, paste0("enrichment_pathway_", clean_name, ".csv"), config, "tables")
    }
    if (!is.null(contrast_results$gmt)) {
      save_table(contrast_results$gmt, paste0("enrichment_gmt_", clean_name, ".csv"), config, "tables")
    }
    if (!is.null(contrast_results$class)) {
      save_table(contrast_results$class, paste0("enrichment_class_", clean_name, ".csv"), config, "tables")
    }
  }

  # Create plots
  plots <- create_enrichment_plots(all_results, config)

  log_message("=== Enrichment Analysis Complete ===")

  list(
    results = all_results,
    plots = plots
  )
}

#' Run pathway enrichment using custom mapping
#'
#' @param de_table DE results table
#' @param pathway_mapping Pathway mapping data frame
#' @param config Configuration
#' @return Enrichment results
run_pathway_enrichment <- function(de_table, pathway_mapping, config) {
  log_message("Running pathway enrichment...")

  # Identify the ID column to use for matching
  id_cols <- c("metabolite_id", "hmdb_id", "kegg_id", "feature_id")
  match_col <- NULL

  for (col in id_cols) {
    if (col %in% colnames(pathway_mapping) && col %in% colnames(de_table)) {
      match_col <- col
      break
    }
  }

  if (is.null(match_col)) {
    # Try matching by feature_id
    match_col <- "feature_id"
    if (!"metabolite_id" %in% colnames(pathway_mapping)) {
      log_message("Cannot match DE results to pathway mapping. Skipping.")
      return(NULL)
    }
    pathway_mapping$feature_id <- pathway_mapping$metabolite_id
  }

  # Get significant and all tested metabolites
  sig_metabolites <- de_table[[match_col]][de_table$significant & !is.na(de_table[[match_col]])]
  all_metabolites <- de_table[[match_col]][!is.na(de_table[[match_col]])]

  if (length(sig_metabolites) < 2) {
    log_message("Too few significant metabolites for enrichment")
    return(NULL)
  }

  # Get unique pathways
  pathways <- unique(pathway_mapping$pathway_id)

  # Run ORA for each pathway
  results <- lapply(pathways, function(pw_id) {
    pw_metabolites <- pathway_mapping[[match_col]][pathway_mapping$pathway_id == pw_id]
    pw_metabolites <- unique(pw_metabolites[!is.na(pw_metabolites)])

    # Overlap with significant
    sig_in_pw <- length(intersect(sig_metabolites, pw_metabolites))
    sig_not_in_pw <- length(setdiff(sig_metabolites, pw_metabolites))

    # Background
    bg_in_pw <- length(intersect(all_metabolites, pw_metabolites)) - sig_in_pw
    bg_not_in_pw <- length(all_metabolites) - sig_in_pw - sig_not_in_pw - bg_in_pw

    # Contingency table
    mat <- matrix(c(sig_in_pw, sig_not_in_pw, bg_in_pw, bg_not_in_pw), nrow = 2)

    if (any(mat < 0)) return(NULL)
    if (sig_in_pw == 0) return(NULL)

    test <- fisher.test(mat, alternative = "greater")

    # Get pathway name
    pw_name <- pathway_mapping$pathway_name[pathway_mapping$pathway_id == pw_id][1]

    data.frame(
      pathway_id = pw_id,
      pathway_name = pw_name,
      n_pathway = length(pw_metabolites),
      n_overlap = sig_in_pw,
      p.value = test$p.value,
      odds_ratio = test$estimate,
      metabolites = paste(intersect(sig_metabolites, pw_metabolites), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })

  results <- do.call(rbind, results[!sapply(results, is.null)])

  if (is.null(results) || nrow(results) == 0) {
    return(NULL)
  }

  results$p.adjust <- p.adjust(results$p.value, method = "BH")
  results <- results[order(results$p.adjust), ]

  log_message("Found ", sum(results$p.adjust < 0.1), " pathways with FDR < 0.1")

  return(results)
}

#' Run GMT-based enrichment
#'
#' @param de_table DE results table
#' @param metabolite_sets Named list of metabolite sets
#' @param config Configuration
#' @return Enrichment results
run_gmt_enrichment <- function(de_table, metabolite_sets, config) {
  log_message("Running GMT-based enrichment...")

  sig_features <- de_table$feature_id[de_table$significant]
  all_features <- de_table$feature_id

  if (length(sig_features) < 2) {
    log_message("Too few significant features for GMT enrichment")
    return(NULL)
  }

  min_size <- config$enrichment$min_set_size %||% 3
  max_size <- config$enrichment$max_set_size %||% 500

  # Filter sets by size
  set_sizes <- sapply(metabolite_sets, length)
  valid_sets <- names(metabolite_sets)[set_sizes >= min_size & set_sizes <= max_size]

  if (length(valid_sets) == 0) {
    log_message("No metabolite sets with valid size")
    return(NULL)
  }

  results <- lapply(valid_sets, function(set_name) {
    set_members <- metabolite_sets[[set_name]]

    sig_in_set <- length(intersect(sig_features, set_members))
    sig_not_in_set <- length(setdiff(sig_features, set_members))
    bg_in_set <- length(intersect(all_features, set_members)) - sig_in_set
    bg_not_in_set <- length(all_features) - sig_in_set - sig_not_in_set - bg_in_set

    mat <- matrix(c(sig_in_set, sig_not_in_set, bg_in_set, bg_not_in_set), nrow = 2)

    if (any(mat < 0) || sig_in_set == 0) return(NULL)

    test <- fisher.test(mat, alternative = "greater")

    data.frame(
      set_name = set_name,
      n_set = length(set_members),
      n_overlap = sig_in_set,
      p.value = test$p.value,
      odds_ratio = test$estimate,
      features = paste(intersect(sig_features, set_members), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })

  results <- do.call(rbind, results[!sapply(results, is.null)])

  if (is.null(results) || nrow(results) == 0) {
    return(NULL)
  }

  results$p.adjust <- p.adjust(results$p.value, method = "BH")
  results <- results[order(results$p.adjust), ]

  return(results)
}

#' Run chemical class enrichment
#'
#' @param de_table DE results table
#' @param annotation_table Annotation table with class column
#' @param config Configuration
#' @return Class enrichment results
run_class_enrichment <- function(de_table, annotation_table, config) {
  log_message("Running chemical class enrichment...")

  if (!"class" %in% colnames(annotation_table)) {
    return(NULL)
  }

  # Match annotations to DE table
  match_idx <- match(de_table$feature_id, annotation_table$feature_id)
  de_table$class <- annotation_table$class[match_idx]

  # Remove entries without class annotation
  de_with_class <- de_table[!is.na(de_table$class), ]

  if (nrow(de_with_class) < 10) {
    log_message("Too few features with class annotation")
    return(NULL)
  }

  sig_features <- de_with_class[de_with_class$significant, ]

  if (nrow(sig_features) < 2) {
    log_message("Too few significant features with class annotation")
    return(NULL)
  }

  # Get unique classes
  classes <- unique(de_with_class$class)

  results <- lapply(classes, function(cls) {
    n_sig_in_class <- sum(sig_features$class == cls)
    n_sig_not_in_class <- nrow(sig_features) - n_sig_in_class
    n_bg_in_class <- sum(de_with_class$class == cls) - n_sig_in_class
    n_bg_not_in_class <- nrow(de_with_class) - n_sig_in_class - n_sig_not_in_class - n_bg_in_class

    mat <- matrix(c(n_sig_in_class, n_sig_not_in_class, n_bg_in_class, n_bg_not_in_class), nrow = 2)

    if (any(mat < 0) || n_sig_in_class == 0) return(NULL)

    test <- fisher.test(mat, alternative = "greater")

    # Direction (up or down enriched)
    sig_in_class <- sig_features[sig_features$class == cls, ]
    pct_up <- mean(sig_in_class$direction == "up") * 100

    data.frame(
      class = cls,
      n_class = sum(de_with_class$class == cls),
      n_significant = n_sig_in_class,
      pct_up = round(pct_up, 1),
      p.value = test$p.value,
      odds_ratio = test$estimate,
      stringsAsFactors = FALSE
    )
  })

  results <- do.call(rbind, results[!sapply(results, is.null)])

  if (is.null(results) || nrow(results) == 0) {
    return(NULL)
  }

  results$p.adjust <- p.adjust(results$p.value, method = "BH")
  results <- results[order(results$p.adjust), ]

  log_message("Found ", sum(results$p.adjust < 0.1), " classes with FDR < 0.1")

  return(results)
}

#' Create enrichment plots
#'
#' @param enrichment_results Enrichment results
#' @param config Configuration
#' @return List of plots
create_enrichment_plots <- function(enrichment_results, config) {
  all_plots <- list()

  for (contrast_name in names(enrichment_results)) {
    clean_name <- gsub(" ", "_", gsub("-", "_vs_", contrast_name))
    results <- enrichment_results[[contrast_name]]
    plots <- list()

    # Pathway enrichment dotplot
    if (!is.null(results$pathway) && nrow(results$pathway) > 0) {
      top_pathways <- head(results$pathway[results$pathway$p.adjust < 0.2, ], 15)

      if (nrow(top_pathways) > 0) {
        top_pathways$short_name <- substr(top_pathways$pathway_name, 1, 40)

        plots$pathway <- ggplot2::ggplot(top_pathways,
          ggplot2::aes(x = -log10(p.adjust), y = reorder(short_name, -log10(p.adjust)),
                       size = n_overlap)) +
          ggplot2::geom_point(color = "#0072B2") +
          ggplot2::labs(
            title = paste("Pathway Enrichment:", contrast_name),
            x = "-Log10(Adjusted P-value)",
            y = "Pathway",
            size = "Overlap"
          ) +
          ggplot2::theme_minimal() +
          ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))

        save_plot(plots$pathway, paste0("enrichment_pathway_", clean_name, ".png"),
                  config, height = 8, subdir = "plots")
      }
    }

    # Class enrichment barplot
    if (!is.null(results$class) && nrow(results$class) > 0) {
      top_classes <- head(results$class[results$class$p.adjust < 0.2, ], 15)

      if (nrow(top_classes) > 0) {
        plots$class <- ggplot2::ggplot(top_classes,
          ggplot2::aes(x = -log10(p.adjust), y = reorder(class, -log10(p.adjust)),
                       fill = pct_up > 50)) +
          ggplot2::geom_bar(stat = "identity") +
          ggplot2::scale_fill_manual(values = c("TRUE" = "#D55E00", "FALSE" = "#0072B2"),
                                      labels = c("Mostly down", "Mostly up"),
                                      name = "Direction") +
          ggplot2::labs(
            title = paste("Chemical Class Enrichment:", contrast_name),
            x = "-Log10(Adjusted P-value)",
            y = "Chemical Class"
          ) +
          ggplot2::theme_minimal()

        save_plot(plots$class, paste0("enrichment_class_", clean_name, ".png"),
                  config, height = 8, subdir = "plots")
      }
    }

    all_plots[[contrast_name]] <- plots
  }

  return(all_plots)
}

#' Create annotation summary
#'
#' @param ingested_data Ingested data
#' @param config Configuration
#' @return Annotation summary
create_annotation_summary <- function(ingested_data, config) {
  log_message("Creating annotation summary...")

  annotations <- ingested_data$annotations
  feature_ids <- rownames(ingested_data$matrix)
  n_features <- length(feature_ids)

  summary_rows <- list()

  summary_rows[[1]] <- data.frame(
    category = "Total features",
    count = n_features,
    percentage = 100,
    stringsAsFactors = FALSE
  )

  if (!is.null(annotations$annotation_table)) {
    ann <- annotations$annotation_table

    # Coverage by ID type
    id_cols <- c("hmdb_id", "kegg_id", "pubchem_id", "chebi_id", "metabolite_name", "class")

    for (col in id_cols) {
      if (col %in% colnames(ann)) {
        matched <- sum(feature_ids %in% ann$feature_id[!is.na(ann[[col]])])
        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          category = paste("With", col),
          count = matched,
          percentage = round(matched / n_features * 100, 1),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  summary_df <- do.call(rbind, summary_rows)
  save_table(summary_df, "annotation_summary.csv", config, "qc")

  return(summary_df)
}
