# =============================================================================
# DIABLO (mixOmics) Supervised Integration
# =============================================================================

#' Run DIABLO integration
run_diablo_integration <- function(feature_data, config) {
  log_message("=== Running DIABLO Integration ===")

  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    log_message("mixOmics package not available. Skipping DIABLO integration.")
    return(NULL)
  }

  # Extract parameters
  diablo_config <- config$integration$diablo
  ncomp <- diablo_config$ncomp %||% 3
  design_type <- diablo_config$design_matrix %||% "full"
  cv_folds <- diablo_config$cv_folds %||% 5
  cv_repeats <- diablo_config$cv_repeats %||% 10
  min_samples <- diablo_config$min_samples_per_group %||% 5

  # Get data
  filtered_mae <- feature_data$filtered_mae
  matrices <- extract_matrices_for_integration(filtered_mae)
  metadata <- get_sample_metadata(filtered_mae)

  condition_col <- config$design$condition_column

  # Check requirements
  if (length(matrices) < 2) {
    log_message("DIABLO requires at least 2 omics layers. Found: ", length(matrices))
    return(NULL)
  }

  if (!condition_col %in% colnames(metadata)) {
    log_message("Condition column '", condition_col, "' not found in metadata")
    return(NULL)
  }

  # Get outcome variable
  Y <- factor(metadata[[condition_col]])
  sample_ids <- rownames(metadata)

  # Check sample sizes per group
  group_sizes <- table(Y)
  if (any(group_sizes < min_samples)) {
    log_message("Warning: Some groups have fewer than ", min_samples, " samples")
    log_message("Group sizes: ", paste(names(group_sizes), "=", group_sizes, collapse = ", "))
  }

  # Prepare data list (samples x features for mixOmics)
  X <- list()
  for (omic in names(matrices)) {
    mat <- matrices[[omic]]
    # Transpose: DIABLO expects samples in rows
    mat_t <- t(mat)
    # Remove features with zero variance
    var_check <- apply(mat_t, 2, var, na.rm = TRUE)
    mat_t <- mat_t[, var_check > 0, drop = FALSE]
    # Impute remaining NAs with column means
    for (j in seq_len(ncol(mat_t))) {
      na_idx <- is.na(mat_t[, j])
      if (any(na_idx)) {
        mat_t[na_idx, j] <- mean(mat_t[, j], na.rm = TRUE)
      }
    }
    X[[omic]] <- mat_t
    log_message("  ", omic, ": ", nrow(mat_t), " samples x ", ncol(mat_t), " features")
  }

  # Create design matrix
  design <- create_diablo_design(names(X), design_type)

  # Limit ncomp to reasonable value
  ncomp <- min(ncomp, min(sapply(X, ncol)), length(levels(Y)) - 1 + 1)
  ncomp <- max(ncomp, 1)

  log_message("Running DIABLO with ", ncomp, " components, design: ", design_type)

  # Run DIABLO
  diablo_result <- tryCatch(
    {
      mixOmics::block.splsda(
        X = X,
        Y = Y,
        ncomp = ncomp,
        design = design
      )
    },
    error = function(e) {
      log_message("Error in DIABLO: ", e$message)
      return(NULL)
    }
  )

  if (is.null(diablo_result)) {
    return(NULL)
  }

  log_message("DIABLO training complete")

  # Extract and save results
  diablo_results <- extract_diablo_results(diablo_result, config)

  # Performance evaluation via cross-validation
  if (cv_folds > 1 && cv_repeats > 0) {
    perf_result <- evaluate_diablo_performance(
      diablo_result, X, Y, design, ncomp, cv_folds, cv_repeats, config
    )
    diablo_results$performance <- perf_result
  }

  # Create visualizations
  create_diablo_plots(diablo_result, diablo_results, metadata, condition_col, config)

  list(
    model = diablo_result,
    results = diablo_results,
    config = diablo_config
  )
}

#' Create DIABLO design matrix
create_diablo_design <- function(block_names, design_type = "full") {
  n_blocks <- length(block_names)
  design <- matrix(0,
    nrow = n_blocks, ncol = n_blocks,
    dimnames = list(block_names, block_names)
  )

  if (design_type == "full") {
    # All blocks connected (correlation = 1)
    design[, ] <- 1
    diag(design) <- 0
  } else if (design_type == "null") {
    # No connections between blocks
    design[, ] <- 0
  } else if (is.numeric(design_type)) {
    # Custom correlation value
    design[, ] <- design_type
    diag(design) <- 0
  }

  design
}

#' Extract DIABLO results
extract_diablo_results <- function(diablo_model, config) {
  log_message("Extracting DIABLO results...")

  # Sample scores per block
  variates <- diablo_model$variates

  # Feature loadings per block
  loadings <- diablo_model$loadings

  # Selected features per component per block
  selected_vars <- list()
  for (block in names(loadings)) {
    if (block == "Y") next
    load_mat <- loadings[[block]]
    selected_vars[[block]] <- rownames(load_mat)[rowSums(load_mat != 0) > 0]
  }

  # Save sample scores
  for (block in names(variates)) {
    scores <- variates[[block]]
    scores_df <- as.data.frame(scores)
    scores_df$sample_id <- rownames(scores)
    scores_df <- scores_df[, c("sample_id", setdiff(colnames(scores_df), "sample_id"))]
    save_table(scores_df, paste0("diablo_scores_", block, ".csv"), config)
  }

  # Save loadings
  for (block in names(loadings)) {
    if (block == "Y") next
    load_df <- as.data.frame(loadings[[block]])
    load_df$feature_id <- rownames(load_df)
    load_df <- load_df[, c("feature_id", setdiff(colnames(load_df), "feature_id"))]
    save_table(load_df, paste0("diablo_loadings_", block, ".csv"), config)
  }

  # Save selected features summary
  selected_summary <- data.frame(
    block = names(selected_vars),
    n_selected = sapply(selected_vars, length),
    stringsAsFactors = FALSE
  )
  save_table(selected_summary, "diablo_selected_features_summary.csv", config)

  list(
    variates = variates,
    loadings = loadings,
    selected_vars = selected_vars,
    Y = diablo_model$Y
  )
}

#' Evaluate DIABLO performance via CV
evaluate_diablo_performance <- function(diablo_model, X, Y, design, ncomp,
                                        cv_folds, cv_repeats, config) {
  log_message(
    "Evaluating DIABLO performance (", cv_folds, "-fold CV x ",
    cv_repeats, " repeats)..."
  )

  perf_result <- tryCatch(
    {
      mixOmics::perf(
        diablo_model,
        validation = "Mfold",
        folds = cv_folds,
        nrepeat = cv_repeats,
        progressBar = FALSE
      )
    },
    error = function(e) {
      log_message("CV performance evaluation failed: ", e$message)
      return(NULL)
    }
  )

  if (is.null(perf_result)) {
    return(NULL)
  }

  # Extract error rates
  error_rates <- perf_result$WeightedVote.error.rate

  if (!is.null(error_rates)) {
    # Convert to data frame
    if (is.list(error_rates)) {
      error_df <- do.call(rbind, lapply(names(error_rates), function(dist) {
        err <- error_rates[[dist]]
        data.frame(
          distance = dist,
          component = seq_len(nrow(err)),
          error_rate = err[, 1],
          stringsAsFactors = FALSE
        )
      }))
    } else {
      error_df <- data.frame(
        component = seq_len(length(error_rates)),
        error_rate = error_rates,
        stringsAsFactors = FALSE
      )
    }

    save_table(error_df, "diablo_cv_error_rates.csv", config)
    log_message("CV error rates saved")
  }

  # Optimal number of components
  optimal_ncomp <- tryCatch(
    {
      perf_result$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]
    },
    error = function(e) ncomp
  )

  log_message("Suggested optimal ncomp: ", optimal_ncomp)

  list(
    perf = perf_result,
    error_rates = error_rates,
    optimal_ncomp = optimal_ncomp
  )
}

#' Create DIABLO visualizations
create_diablo_plots <- function(diablo_model, diablo_results, metadata,
                                condition_col, config) {
  log_message("Creating DIABLO plots...")

  Y <- diablo_results$Y

  # 1. Sample plot (first 2 components)
  tryCatch(
    {
      p <- plot_diablo_samples(diablo_model, Y, comp = c(1, 2))
      save_plot(p, "diablo_sample_plot_comp1_2.png", config, width = 10, height = 8)
    },
    error = function(e) log_message("Failed to create sample plot: ", e$message)
  )

  # 2. Correlation circle plot
  tryCatch(
    {
      p <- plot_diablo_correlation_circle(diablo_model)
      save_plot(p, "diablo_correlation_circle.png", config, width = 10, height = 10)
    },
    error = function(e) log_message("Failed to create correlation plot: ", e$message)
  )

  # 3. Loadings plot per block
  tryCatch(
    {
      for (block in names(diablo_results$loadings)) {
        if (block == "Y") next
        p <- plot_diablo_loadings(diablo_model, block)
        save_plot(p, paste0("diablo_loadings_", block, ".png"), config, width = 10, height = 8)
      }
    },
    error = function(e) log_message("Failed to create loadings plots: ", e$message)
  )



  # 4. Circos plot for feature correlations
  tryCatch(
    {
      circos_file <- file.path(config$output$output_dir, "plots", "diablo_circos.png")

      # Use png device directly as circosPlot plots to device
      png(circos_file, width = 1000, height = 1000, res = 150)
      mixOmics::circosPlot(diablo_model,
        cutoff = 0.7, line = TRUE,
        color.blocks = c("darkorchid", "brown1", "lightgreen"),
        color.cor = c("chocolate3", "grey20"), size.labels = 1.5
      )
      dev.off()

      log_message("Saved circos plot to ", circos_file)
    },
    error = function(e) log_message("Failed to create circos plot: ", e$message)
  )

  # 5. CIM (Clustered Image Map) plot
  tryCatch(
    {
      cim_file <- file.path(config$output$output_dir, "plots", "diablo_cim.png")

      # Use png device directly as cim plots to device
      png(cim_file, width = 1200, height = 1000, res = 150)
      mixOmics::cim(diablo_model,
        comp = 1:2,  # Use first 2 components
        margins = c(10, 5),
        row.names = FALSE,  # Hide row names for clarity with many features
        col.names = TRUE,   # Show sample names
        title = "DIABLO: Clustered Image Map (CIM)",
        legend = TRUE,
        color = colorRampPalette(c("blue", "white", "red"))(100),
        cutoff = 0.5  # Only show features with loading > 0.5
      )
      dev.off()

      log_message("Saved CIM plot to ", cim_file)
    },
    error = function(e) log_message("Failed to create CIM plot: ", e$message)
  )

  # 6. Enrichment on Loadings
  tryCatch(
    {
      run_diablo_enrichment(diablo_results, config)
    },
    error = function(e) log_message("Failed to run DIABLO enrichment: ", e$message)
  )


  log_message("DIABLO plots complete")
}

#' Plot DIABLO sample scores
plot_diablo_samples <- function(diablo_model, Y, comp = c(1, 2)) {
  variates <- diablo_model$variates
  blocks <- setdiff(names(variates), "Y")

  # Combine scores from all blocks
  all_scores <- do.call(cbind, lapply(blocks, function(b) {
    v <- variates[[b]][, comp]
    colnames(v) <- paste0(b, "_comp", comp)
    v
  }))

  # Average across blocks for visualization
  avg_scores <- data.frame(
    comp1 = rowMeans(sapply(blocks, function(b) variates[[b]][, comp[1]])),
    comp2 = rowMeans(sapply(blocks, function(b) variates[[b]][, comp[2]])),
    condition = Y,
    stringsAsFactors = FALSE
  )

  ggplot2::ggplot(avg_scores, ggplot2::aes(x = comp1, y = comp2, color = condition)) +
    ggplot2::geom_point(size = 3, alpha = 0.8) +
    ggplot2::stat_ellipse(level = 0.95) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "DIABLO: Sample Scores",
      x = paste0("Component ", comp[1]),
      y = paste0("Component ", comp[2]),
      color = "Condition"
    )
}

#' Plot correlation circle
plot_diablo_correlation_circle <- function(diablo_model) {
  # Use mixOmics plotVar internally, but create ggplot version
  loadings <- diablo_model$loadings
  blocks <- setdiff(names(loadings), "Y")

  all_loadings <- do.call(rbind, lapply(blocks, function(b) {
    load_mat <- loadings[[b]]
    if (ncol(load_mat) >= 2) {
      df <- data.frame(
        comp1 = load_mat[, 1],
        comp2 = load_mat[, 2],
        block = b,
        feature = rownames(load_mat),
        stringsAsFactors = FALSE
      )
      # Filter to non-zero loadings
      df <- df[df$comp1 != 0 | df$comp2 != 0, ]
      df
    } else {
      NULL
    }
  }))

  if (is.null(all_loadings) || nrow(all_loadings) == 0) {
    return(ggplot2::ggplot() +
      ggplot2::ggtitle("No loadings to display"))
  }

  # Draw correlation circle
  circle <- data.frame(
    x = cos(seq(0, 2 * pi, length.out = 100)),
    y = sin(seq(0, 2 * pi, length.out = 100))
  )

  ggplot2::ggplot(all_loadings) +
    ggplot2::geom_path(data = circle, ggplot2::aes(x = x, y = y), color = "gray50") +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, y = 0, xend = comp1, yend = comp2, color = block),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.1, "cm")),
      alpha = 0.5
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "DIABLO: Correlation Circle",
      x = "Component 1",
      y = "Component 2",
      color = "Block"
    )
}

#' Plot loadings for a block
plot_diablo_loadings <- function(diablo_model, block, comp = 1, n_top = 20) {
  loadings <- diablo_model$loadings[[block]]

  if (ncol(loadings) < comp) {
    return(ggplot2::ggplot() +
      ggplot2::ggtitle("Component not available"))
  }

  load_vec <- loadings[, comp]
  load_vec <- load_vec[load_vec != 0]

  if (length(load_vec) == 0) {
    return(ggplot2::ggplot() +
      ggplot2::ggtitle("No selected features"))
  }

  # Top features by absolute loading
  ord <- order(abs(load_vec), decreasing = TRUE)
  top_idx <- head(ord, n_top)

  df <- data.frame(
    feature = names(load_vec)[top_idx],
    loading = load_vec[top_idx],
    stringsAsFactors = FALSE
  )
  df$feature <- factor(df$feature, levels = df$feature[order(abs(df$loading))])
  df$direction <- ifelse(df$loading > 0, "positive", "negative")

  ggplot2::ggplot(df, ggplot2::aes(x = loading, y = feature, fill = direction)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = c("positive" = "firebrick", "negative" = "steelblue")) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste0("DIABLO: ", block, " Loadings (Component ", comp, ")"),
      x = "Loading",
      y = NULL
    ) +
    ggplot2::theme(legend.position = "none")
}

#' Get discriminating features from DIABLO
get_diablo_discriminating_features <- function(diablo_results, min_loading = 0.1) {
  selected <- diablo_results$selected_vars
  loadings <- diablo_results$loadings

  all_features <- list()

  for (block in names(selected)) {
    load_mat <- loadings[[block]]
    feat_ids <- selected[[block]]

    # Get max absolute loading across components for each feature
    max_loadings <- apply(abs(load_mat[feat_ids, , drop = FALSE]), 1, max)

    # Filter by minimum loading
    keep <- max_loadings >= min_loading

    all_features[[block]] <- data.frame(
      feature_id = feat_ids[keep],
      max_loading = max_loadings[keep],
      block = block,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, all_features)
}

#' Run DIABLO enrichment on loadings
run_diablo_enrichment <- function(diablo_results, config) {
  log_message("Running enrichment on DIABLO loadings...")

  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    log_message("clusterProfiler not installed, skipping enrichment")
    return(NULL)
  }

  loadings <- diablo_results$loadings
  enrich_results <- list()

  # Identify organism-specific settings
  organism_name <- config$global$organism %||% "human"
  org_db <- NULL
  kegg_code <- NULL

  if (organism_name == "c_elegans") {
    if (requireNamespace("org.Ce.eg.db", quietly = TRUE)) {
      org_db <- org.Ce.eg.db::org.Ce.eg.db
      kegg_code <- "cel"
    }
  } else {
    if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      org_db <- org.Hs.eg.db::org.Hs.eg.db
      kegg_code <- "hsa"
    }
  }

  if (is.null(org_db)) {
    log_message("OrgDb not found for ", organism_name, ", skipping enrichment")
    return(NULL)
  }

  # Iterate through blocks (omics)
  for (block in names(loadings)) {
    if (block == "Y") next

    # Skip metabolomics for gene enrichment
    if (block == "metabolomics") {
      log_message("  Skipping gene enrichment for metabolomics block")
      next
    }

    log_message("  Processing block: ", block)
    load_mat <- loadings[[block]]

    # Iterate through component columns
    for (k in seq_len(ncol(load_mat))) {
      comp_name <- paste0("comp", k)

      # Get top 100 features by absolute weight
      w <- load_mat[, k]
      top_features <- names(sort(abs(w), decreasing = TRUE))[1:100]

      log_message("    Top features example: ", paste(head(top_features, 3), collapse = ", "))

      # Convert ID type based on block
      key_type <- "SYMBOL" # Default

      # Heuristic: Check ID structure
      if (all(grepl("^WBGene", top_features[1:5], ignore.case = TRUE))) {
        key_type <- "WORMBASE"
      } else if (all(grepl("^ENS", top_features[1:5], ignore.case = TRUE))) {
        key_type <- "ENSEMBL"
      } else if (any(grepl("^[OPQ][0-9]", top_features[1:5], ignore.case = TRUE)) ||
        any(grepl("^A0A", top_features[1:5], ignore.case = TRUE))) {
        key_type <- "UNIPROT"
      }

      log_message("    Detected keyType: ", key_type)

      tryCatch(
        {
          gene_map <- clusterProfiler::bitr(top_features,
            fromType = key_type,
            toType = "ENTREZID", OrgDb = org_db
          )

          if (!is.null(gene_map) && nrow(gene_map) > 5) {
            genes <- gene_map$ENTREZID

            # 1. KEGG Enrichment
            tryCatch(
              {
                c_kegg <- clusterProfiler::enrichKEGG(gene = genes, organism = kegg_code, pvalueCutoff = 0.05)

                if (!is.null(c_kegg) && nrow(as.data.frame(c_kegg)) > 0) {
                  res_kegg <- as.data.frame(c_kegg)
                  res_kegg$block <- block
                  res_kegg$component <- comp_name
                  res_kegg$type <- "KEGG"

                  out_name <- paste0("diablo_enrichment_kegg_", block, "_", comp_name)
                  save_table(res_kegg, paste0(out_name, ".csv"), config)

                  p <- clusterProfiler::dotplot(c_kegg, showCategory = 15) +
                    ggplot2::ggtitle(paste0("KEGG: ", block, " - ", comp_name))
                  save_plot(p, paste0(out_name, ".png"), config)

                  enrich_results[[out_name]] <- res_kegg
                  log_message("    KEGG success: ", nrow(res_kegg), " pathways")
                }
              },
              error = function(e) log_message("    KEGG failed: ", e$message)
            )

            # 2. GO Enrichment (BP) - Robust Fallback
            tryCatch(
              {
                c_go <- clusterProfiler::enrichGO(gene = genes, OrgDb = org_db, ont = "BP", pvalueCutoff = 0.05)

                if (!is.null(c_go) && nrow(as.data.frame(c_go)) > 0) {
                  res_go <- as.data.frame(c_go)
                  res_go$block <- block
                  res_go$component <- comp_name
                  res_go$type <- "GO"

                  out_name <- paste0("diablo_enrichment_go_", block, "_", comp_name)
                  save_table(res_go, paste0(out_name, ".csv"), config)

                  p <- clusterProfiler::dotplot(c_go, showCategory = 15) +
                    ggplot2::ggtitle(paste0("GO BP: ", block, " - ", comp_name))
                  save_plot(p, paste0(out_name, ".png"), config)

                  enrich_results[[out_name]] <- res_go
                  log_message("    GO success: ", nrow(res_go), " terms")
                }
              },
              error = function(e) log_message("    GO failed: ", e$message)
            )

            # 3. Reactome Enrichment (if applicable)
            if (requireNamespace("ReactomePA", quietly = TRUE) && organism_name == "c_elegans") {
              tryCatch(
                {
                  c_react <- ReactomePA::enrichPathway(gene = genes, organism = "celegans", pvalueCutoff = 0.05)
                  if (!is.null(c_react) && nrow(as.data.frame(c_react)) > 0) {
                    res_react <- as.data.frame(c_react)
                    res_react$block <- block
                    res_react$component <- comp_name
                    res_react$type <- "Reactome"

                    out_name <- paste0("diablo_enrichment_reactome_", block, "_", comp_name)
                    save_table(res_react, paste0(out_name, ".csv"), config)

                    p <- clusterProfiler::dotplot(c_react, showCategory = 15) +
                      ggplot2::ggtitle(paste0("Reactome: ", block, " - ", comp_name))
                    save_plot(p, paste0(out_name, ".png"), config)

                    enrich_results[[out_name]] <- res_react
                    log_message("    Reactome success: ", nrow(res_react), " pathways")
                  }
                },
                error = function(e) log_message("    Reactome failed: ", e$message)
              )
            }
          } else {
            log_message("    Could not map enough genes to ENTREZID")
          }
        },
        error = function(e) {
          log_message("    Enrichment error for ", block, " ", comp_name, ": ", e$message)
        }
      )
    }
  }

  return(enrich_results)
}
