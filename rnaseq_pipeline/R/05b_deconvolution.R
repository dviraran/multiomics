# =============================================================================
# Cell Type Deconvolution Analysis
# =============================================================================
# Bulk RNA-seq cell type deconvolution using xCell 2.0
# Supports configurable references and non-human organisms
#
# Key analyses:
# - xCell 2.0 cell type enrichment analysis
# - Differential composition analysis across conditions
# - Cell type correlation with phenotypes
# - Visualization of cell type proportions
# =============================================================================

# -----------------------------------------------------------------------------
# Main Orchestrator Function
# -----------------------------------------------------------------------------

#' Run Cell Type Deconvolution Analysis
#'
#' Comprehensive cell type deconvolution using xCell 2.0
#'
#' @param dds DESeqDataSet object with normalized counts
#' @param de_results Data frame with DE results
#' @param config Pipeline configuration list
#' @return List containing deconvolution results and plots
#' @export
run_deconvolution_analysis <- function(dds, de_results, config) {
  log_message("=== Running Cell Type Deconvolution Analysis ===")

  # Get config settings
  deconv_config <- config$deconvolution %||% list()

  if (!(deconv_config$run_deconvolution %||% TRUE)) {
    log_message("Deconvolution disabled in config. Skipping.")
    return(NULL)
  }

  # Validate inputs
  if (is.null(dds) || length(dds) == 0) {
    log_message("WARNING: Invalid dds object. Skipping deconvolution.")
    return(NULL)
  }

  # Set up output directories
  output_dir <- file.path(config$output$output_dir %||% "outputs", "deconvolution")
  plots_dir <- file.path(output_dir, "plots")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

  # Check for xCell2 package
  if (!requireNamespace("xCell2", quietly = TRUE)) {
    log_message("WARNING: xCell2 package not installed. Install from GitHub: dviraran/xCell2")
    log_message("Attempting to use fallback xCell if available...")

    if (!requireNamespace("xCell", quietly = TRUE)) {
      log_message("WARNING: Neither xCell2 nor xCell available. Skipping deconvolution.")
      return(NULL)
    }

    # Use legacy xCell
    return(tryCatch(
      run_xcell_legacy(dds, de_results, config, output_dir, plots_dir),
      error = function(e) {
        log_message("WARNING: xCell legacy analysis failed: ", conditionMessage(e))
        NULL
      }
    ))
  }

  # Extract expression matrix
  expr_matrix <- get_expression_matrix(dds, config)

  # Get metadata
  metadata <- as.data.frame(SummarizedExperiment::colData(dds))

  # Run xCell 2.0
  xcell2_results <- run_xcell2_analysis(
    expr_matrix = expr_matrix,
    metadata = metadata,
    config = deconv_config,
    output_dir = output_dir
  )

  if (is.null(xcell2_results)) {
    log_message("xCell2 analysis failed. Returning NULL.")
    return(NULL)
  }

  # Differential composition analysis
  diff_comp_results <- NULL
  if (deconv_config$differential_composition %||% TRUE) {
    diff_comp_results <- analyze_differential_composition(
      xcell2_scores = xcell2_results$scores,
      metadata = metadata,
      config = config,
      output_dir = output_dir
    )
  }

  # Cell type-phenotype correlations
  phenotype_corr <- NULL
  if (deconv_config$phenotype_correlations %||% TRUE) {
    phenotype_corr <- correlate_celltype_phenotypes(
      xcell2_scores = xcell2_results$scores,
      metadata = metadata,
      config = config,
      output_dir = output_dir
    )
  }

  # Generate visualizations
  figures <- generate_deconvolution_plots(
    xcell2_results = xcell2_results,
    diff_comp_results = diff_comp_results,
    phenotype_corr = phenotype_corr,
    metadata = metadata,
    config = config,
    plots_dir = plots_dir
  )

  # Compile results
  results <- list(
    scores = xcell2_results$scores,
    scores_raw = xcell2_results$scores_raw,
    reference_used = xcell2_results$reference,
    organism = deconv_config$organism %||% "human",
    differential_composition = diff_comp_results,
    phenotype_correlations = phenotype_corr,
    figures = figures,
    summary = create_deconvolution_summary(xcell2_results, diff_comp_results, phenotype_corr)
  )

  # Save outputs
  save_deconvolution_outputs(results, output_dir)

  log_message("Cell type deconvolution analysis complete!")
  return(results)
}

# -----------------------------------------------------------------------------
# xCell 2.0 Analysis
# -----------------------------------------------------------------------------

#' Run xCell 2.0 Analysis
#'
#' @param expr_matrix Expression matrix (genes x samples)
#' @param metadata Sample metadata
#' @param config Deconvolution config
#' @param output_dir Output directory
#' @return List with xCell2 scores and metadata
run_xcell2_analysis <- function(expr_matrix, metadata, config, output_dir) {
  log_message("Running xCell 2.0 analysis...")

  # Get configuration parameters
  organism_raw <- config$organism %||% "human"
  organism <- normalize_organism(organism_raw)
  reference <- config$reference %||% get_default_reference(organism)
  custom_reference <- config$custom_reference %||% NULL

  log_message("  Organism: ", organism)
  log_message("  Reference: ", reference)

  tryCatch({
    # Load xCell2
    library(xCell2)

    # Determine reference to use
    ref_object <- NULL

    if (!is.null(custom_reference)) {
      # Load custom reference
      log_message("  Loading custom reference from: ", custom_reference)
      if (file.exists(custom_reference)) {
        ref_object <- readRDS(custom_reference)
      } else {
        log_message("WARNING: Custom reference file not found. Using default.")
      }
    }

    # Run xCell2 with appropriate settings
    if (!is.null(ref_object)) {
      # Use custom reference
      scores <- xCell2::xCell2Analysis(
        expr = expr_matrix,
        xcell2_ref = ref_object,
        parallel = config$parallel %||% TRUE,
        BPPARAM = BiocParallel::MulticoreParam(workers = config$n_cores %||% 4)
      )
    } else if (organism == "human") {
      # Use built-in human reference
      scores <- run_xcell2_human(expr_matrix, reference, config)
    } else if (organism == "mouse") {
      # Use mouse reference
      scores <- run_xcell2_mouse(expr_matrix, reference, config)
    } else {
      # Attempt with custom organism settings
      log_message("WARNING: Organism '", organism, "' may require custom reference")
      scores <- run_xcell2_human(expr_matrix, reference, config)
    }

    # Process scores
    scores_df <- as.data.frame(t(scores))
    scores_df$sample_id <- rownames(scores_df)

    # Filter low-confidence cell types
    min_score <- config$min_score %||% 0.01
    scores_filtered <- filter_low_confidence_celltypes(scores_df, min_score)

    log_message("  Detected ", ncol(scores_filtered) - 1, " cell types")

    return(list(
      scores = scores_filtered,
      scores_raw = scores_df,
      reference = reference,
      organism = organism
    ))

  }, error = function(e) {
    log_message("ERROR in xCell2 analysis: ", e$message)
    return(NULL)
  })
}

#' Run xCell2 with Human Reference
run_xcell2_human <- function(expr_matrix, reference, config) {
  log_message("  Using human reference: ", reference)

  # Map reference name to xCell2 data object name (new API)
  ref_map <- c(
    "BlueprintEncode" = "BlueprintEncode.xCell2Ref",
    "LM22" = "LM22.xCell2Ref",
    "DICE" = "DICE_demo.xCell2Ref",
    "ImmuneCompendium" = "ImmuneCompendium.xCell2Ref",
    "PanCancer" = "PanCancer.xCell2Ref",
    "TMECompendium" = "TMECompendium.xCell2Ref",
    "TabulaSapiensBlood" = "TabulaSapiensBlood.xCell2Ref"
  )

  # Get data object name
  ref_data_name <- ref_map[[reference]]
  if (is.null(ref_data_name)) {
    log_message("WARNING: Reference '", reference, "' not recognized. Using BlueprintEncode.")
    ref_data_name <- "BlueprintEncode.xCell2Ref"
  }

  # Get reference object
  ref_obj <- get_xcell2_reference(ref_data_name)

  if (is.null(ref_obj)) {
    log_message("  Could not load xCell2 reference. Skipping deconvolution.")
    return(NULL)
  }

  # Run analysis with new API: mix, xcell2object
  scores <- tryCatch({
    xCell2::xCell2Analysis(
      mix = expr_matrix,
      xcell2object = ref_obj,
      minSharedGenes = 0.8,
      rawScores = FALSE,
      spillover = TRUE
    )
  }, error = function(e) {
    log_message("  xCell2 analysis failed: ", e$message)
    NULL
  })

  return(scores)
}

#' Run xCell2 with Mouse Reference
run_xcell2_mouse <- function(expr_matrix, reference, config) {
  log_message("  Using mouse reference: ", reference)

  # Map mouse reference names to xCell2 data object names
  mouse_ref_map <- c(
    "ImmGen" = "ImmGenData.xCell2Ref",
    "mouse_default" = "ImmGenData.xCell2Ref",
    "MouseRNAseq" = "MouseRNAseqData.xCell2Ref",
    "TabulaMurisBlood" = "TabulaMurisBlood.xCell2Ref"
  )

  ref_data_name <- mouse_ref_map[[reference]]
  if (!is.null(ref_data_name)) {
    ref_obj <- get_xcell2_reference(ref_data_name)

    if (!is.null(ref_obj)) {
      scores <- tryCatch({
        xCell2::xCell2Analysis(
          mix = expr_matrix,
          xcell2object = ref_obj,
          minSharedGenes = 0.8,
          rawScores = FALSE,
          spillover = TRUE
        )
      }, error = function(e) {
        log_message("  xCell2 mouse analysis failed: ", e$message)
        NULL
      })
      return(scores)
    }
  }

  # Fallback: use homolog mapping to human reference
  log_message("  Mapping mouse genes to human orthologs for analysis...")
  expr_human <- convert_mouse_to_human(expr_matrix)

  ref_obj <- get_xcell2_reference("BlueprintEncode.xCell2Ref")

  if (is.null(ref_obj)) {
    log_message("  Could not load human reference for mouse fallback.")
    return(NULL)
  }

  scores <- tryCatch({
    xCell2::xCell2Analysis(
      mix = expr_human,
      xcell2object = ref_obj,
      minSharedGenes = 0.8,
      rawScores = FALSE,
      spillover = TRUE
    )
  }, error = function(e) {
    log_message("  xCell2 mouse fallback analysis failed: ", e$message)
    NULL
  })

  return(scores)
}

#' Get xCell2 Reference Object
get_xcell2_reference <- function(reference_name) {
  # New xCell2 API uses data objects like "BlueprintEncode.xCell2Ref"
  ref_obj <- NULL

  # Approach 1: Load from package data (new API - xCell2 >= 2.0)
  ref_obj <- tryCatch({
    data_env <- new.env()
    # Reference names are like "BlueprintEncode.xCell2Ref"
    data(list = reference_name, package = "xCell2", envir = data_env)
    get(reference_name, envir = data_env)
  }, error = function(e) NULL)

  # Approach 2: Try older naming convention
  if (is.null(ref_obj)) {
    ref_obj <- tryCatch({
      data_env <- new.env()
      # Try without .xCell2Ref suffix
      base_name <- sub("\\.xCell2Ref$", "", reference_name)
      data(list = base_name, package = "xCell2", envir = data_env)
      get(base_name, envir = data_env)
    }, error = function(e) NULL)
  }

  # Approach 3: Try getXCell2Ref function (older API)
  if (is.null(ref_obj)) {
    ref_obj <- tryCatch({
      if (exists("getXCell2Ref", where = asNamespace("xCell2"), mode = "function")) {
        base_name <- sub("\\.xCell2Ref$", "", reference_name)
        xCell2::getXCell2Ref(base_name)
      } else {
        NULL
      }
    }, error = function(e) NULL)
  }

  if (is.null(ref_obj)) {
    log_message("  Could not load reference '", reference_name, "'.")
  }

  return(ref_obj)
}

#' Convert Mouse Genes to Human Orthologs
convert_mouse_to_human <- function(expr_matrix) {
  log_message("  Converting mouse genes to human orthologs...")

  # Use biomaRt for ortholog mapping
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    log_message("WARNING: biomaRt not available. Attempting simple capitalization...")
    # Simple fallback: capitalize gene names (works for many genes)
    rownames(expr_matrix) <- toupper(rownames(expr_matrix))
    return(expr_matrix)
  }

  tryCatch({
    # Get ortholog mapping
    mouse_mart <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")

    orthologs <- biomaRt::getBM(
      attributes = c("mgi_symbol", "hsapiens_homolog_associated_gene_name"),
      filters = "mgi_symbol",
      values = rownames(expr_matrix),
      mart = mouse_mart
    )

    # Filter valid mappings
    orthologs <- orthologs[orthologs$hsapiens_homolog_associated_gene_name != "", ]

    # Map and aggregate (take mean for duplicates)
    expr_mapped <- expr_matrix[rownames(expr_matrix) %in% orthologs$mgi_symbol, , drop = FALSE]

    # Create mapping
    gene_map <- setNames(
      orthologs$hsapiens_homolog_associated_gene_name,
      orthologs$mgi_symbol
    )

    # Rename rows
    new_names <- gene_map[rownames(expr_mapped)]
    rownames(expr_mapped) <- new_names

    # Aggregate duplicates
    expr_agg <- aggregate_duplicate_genes(expr_mapped)

    log_message("  Mapped ", nrow(expr_agg), " genes to human orthologs")
    return(expr_agg)

  }, error = function(e) {
    log_message("WARNING: Ortholog mapping failed: ", e$message)
    log_message("  Using simple capitalization fallback...")
    rownames(expr_matrix) <- toupper(rownames(expr_matrix))
    return(expr_matrix)
  })
}

#' Aggregate Duplicate Genes
aggregate_duplicate_genes <- function(expr_matrix) {
  if (!any(duplicated(rownames(expr_matrix)))) {
    return(expr_matrix)
  }

  # Take mean of duplicates
  expr_df <- as.data.frame(expr_matrix)
  expr_df$gene <- rownames(expr_df)

  expr_agg <- expr_df %>%
    group_by(gene) %>%
    summarise(across(everything(), mean, na.rm = TRUE)) %>%
    as.data.frame()

  rownames(expr_agg) <- expr_agg$gene
  expr_agg$gene <- NULL

  as.matrix(expr_agg)
}

#' Filter Low-Confidence Cell Types
filter_low_confidence_celltypes <- function(scores_df, min_score = 0.01) {
  sample_col <- which(colnames(scores_df) == "sample_id")
  score_cols <- setdiff(1:ncol(scores_df), sample_col)

  # Calculate max score per cell type
  max_scores <- apply(scores_df[, score_cols, drop = FALSE], 2, max, na.rm = TRUE)

  # Keep cell types with at least one sample above threshold
  keep_types <- names(max_scores)[max_scores >= min_score]

  # Return filtered data frame
  scores_df[, c(keep_types, "sample_id")]
}

# -----------------------------------------------------------------------------
# Differential Composition Analysis
# -----------------------------------------------------------------------------

#' Analyze Differential Composition Across Conditions
#'
#' @param xcell2_scores xCell2 scores data frame
#' @param metadata Sample metadata
#' @param config Pipeline configuration
#' @param output_dir Output directory
#' @return Data frame with differential composition results
analyze_differential_composition <- function(xcell2_scores, metadata, config, output_dir) {
  log_message("Analyzing differential cell type composition...")

  # Get condition column
  condition_col <- config$design$condition_column %||% "condition"

  if (!condition_col %in% colnames(metadata)) {
    log_message("WARNING: Condition column '", condition_col, "' not found. Skipping differential analysis.")
    return(NULL)
  }

  # Merge scores with metadata
  scores_meta <- merge(
    xcell2_scores,
    metadata[, c("sample_id", condition_col), drop = FALSE],
    by.x = "sample_id",
    by.y = "sample_id"
  )

  # Get cell type columns
  celltype_cols <- setdiff(colnames(xcell2_scores), "sample_id")

  # Perform tests for each cell type
  diff_results <- lapply(celltype_cols, function(ct) {
    test_celltype_differential(scores_meta, ct, condition_col, config)
  })

  diff_df <- do.call(rbind, diff_results)

  # Adjust p-values
  diff_df$padj <- p.adjust(diff_df$pvalue, method = "BH")

  # Sort by significance
  diff_df <- diff_df[order(diff_df$padj), ]

  # Save results
  write.csv(
    diff_df,
    file.path(output_dir, "differential_composition.csv"),
    row.names = FALSE
  )

  log_message("  Found ", sum(diff_df$padj < 0.05, na.rm = TRUE),
              " significantly different cell types (FDR < 0.05)")

  return(diff_df)
}

#' Test Single Cell Type for Differential Abundance
test_celltype_differential <- function(scores_meta, celltype, condition_col, config) {
  # Get values
  values <- scores_meta[[celltype]]
  conditions <- scores_meta[[condition_col]]

  # Get reference level
  ref_level <- config$design$reference_level %||% levels(factor(conditions))[1]

  # Get unique conditions
  cond_levels <- unique(conditions)

  if (length(cond_levels) == 2) {
    # Two-group comparison (Wilcoxon)
    group1 <- values[conditions == ref_level]
    group2 <- values[conditions != ref_level]

    test_result <- tryCatch({
      wilcox.test(group2, group1)
    }, error = function(e) {
      list(p.value = NA, statistic = NA)
    })

    # Calculate effect size (Cliff's delta)
    effect_size <- tryCatch({
      calculate_cliffs_delta(group2, group1)
    }, error = function(e) NA)

    # Calculate log2 fold change of means
    mean1 <- mean(group1, na.rm = TRUE)
    mean2 <- mean(group2, na.rm = TRUE)
    log2fc <- log2((mean2 + 0.001) / (mean1 + 0.001))

    return(data.frame(
      cell_type = celltype,
      comparison = paste0(setdiff(cond_levels, ref_level), "_vs_", ref_level),
      mean_ref = mean1,
      mean_alt = mean2,
      log2fc = log2fc,
      effect_size = effect_size,
      statistic = test_result$statistic,
      pvalue = test_result$p.value,
      stringsAsFactors = FALSE
    ))

  } else {
    # Multi-group comparison (Kruskal-Wallis)
    test_result <- tryCatch({
      kruskal.test(values ~ factor(conditions))
    }, error = function(e) {
      list(p.value = NA, statistic = NA)
    })

    # Calculate eta-squared as effect size
    effect_size <- tryCatch({
      calculate_eta_squared_kw(values, conditions)
    }, error = function(e) NA)

    return(data.frame(
      cell_type = celltype,
      comparison = "multi_group",
      mean_ref = mean(values[conditions == ref_level], na.rm = TRUE),
      mean_alt = mean(values[conditions != ref_level], na.rm = TRUE),
      log2fc = NA,
      effect_size = effect_size,
      statistic = test_result$statistic,
      pvalue = test_result$p.value,
      stringsAsFactors = FALSE
    ))
  }
}

#' Calculate Cliff's Delta Effect Size
calculate_cliffs_delta <- function(x, y) {
  nx <- length(x)
  ny <- length(y)

  dominance <- outer(x, y, function(a, b) sign(a - b))

  sum(dominance) / (nx * ny)
}

#' Calculate Eta-Squared for Kruskal-Wallis
calculate_eta_squared_kw <- function(values, groups) {
  n <- length(values)
  k <- length(unique(groups))

  H <- kruskal.test(values ~ factor(groups))$statistic

  # Eta-squared approximation
  (H - k + 1) / (n - k)
}

# -----------------------------------------------------------------------------
# Cell Type-Phenotype Correlations
# -----------------------------------------------------------------------------

#' Correlate Cell Types with Phenotypes
#'
#' @param xcell2_scores xCell2 scores data frame
#' @param metadata Sample metadata
#' @param config Pipeline configuration
#' @param output_dir Output directory
#' @return Data frame with correlation results
correlate_celltype_phenotypes <- function(xcell2_scores, metadata, config, output_dir) {
  log_message("Correlating cell types with phenotypes...")

  # Get phenotype columns
  phenotype_cols <- config$deconvolution$phenotype_columns %||% NULL

  if (is.null(phenotype_cols)) {
    # Auto-detect numeric columns in metadata
    numeric_cols <- sapply(metadata, is.numeric)
    phenotype_cols <- names(metadata)[numeric_cols]
    phenotype_cols <- setdiff(phenotype_cols, c("sample_id", "batch"))
  }

  if (length(phenotype_cols) == 0) {
    log_message("  No phenotype columns found for correlation analysis.")
    return(NULL)
  }

  log_message("  Phenotype columns: ", paste(phenotype_cols, collapse = ", "))

  # Get cell type columns
  celltype_cols <- setdiff(colnames(xcell2_scores), "sample_id")

  # Merge scores with metadata
  merged <- merge(
    xcell2_scores,
    metadata,
    by.x = "sample_id",
    by.y = "sample_id"
  )

  # Compute correlations
  corr_results <- list()

  for (ct in celltype_cols) {
    for (pheno in phenotype_cols) {
      if (!pheno %in% colnames(merged)) next

      ct_vals <- merged[[ct]]
      pheno_vals <- merged[[pheno]]

      # Skip if not enough valid values
      valid_idx <- complete.cases(ct_vals, pheno_vals)
      if (sum(valid_idx) < 5) next

      # Spearman correlation
      cor_test <- tryCatch({
        cor.test(ct_vals[valid_idx], pheno_vals[valid_idx], method = "spearman")
      }, error = function(e) NULL)

      if (!is.null(cor_test)) {
        corr_results[[length(corr_results) + 1]] <- data.frame(
          cell_type = ct,
          phenotype = pheno,
          correlation = cor_test$estimate,
          pvalue = cor_test$p.value,
          n_samples = sum(valid_idx),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(corr_results) == 0) {
    log_message("  No valid correlations computed.")
    return(NULL)
  }

  corr_df <- do.call(rbind, corr_results)

  # Adjust p-values
  corr_df$padj <- p.adjust(corr_df$pvalue, method = "BH")

  # Sort by significance
  corr_df <- corr_df[order(corr_df$padj), ]

  # Save results
  write.csv(
    corr_df,
    file.path(output_dir, "celltype_phenotype_correlations.csv"),
    row.names = FALSE
  )

  log_message("  Found ", sum(corr_df$padj < 0.05, na.rm = TRUE),
              " significant correlations (FDR < 0.05)")

  return(corr_df)
}

# -----------------------------------------------------------------------------
# Visualization Functions
# -----------------------------------------------------------------------------

#' Generate Deconvolution Plots
#'
#' @param xcell2_results xCell2 analysis results
#' @param diff_comp_results Differential composition results
#' @param phenotype_corr Phenotype correlation results
#' @param metadata Sample metadata
#' @param config Pipeline configuration
#' @param plots_dir Plots output directory
#' @return List of figure paths
generate_deconvolution_plots <- function(xcell2_results, diff_comp_results, phenotype_corr,
                                          metadata, config, plots_dir) {
  log_message("Generating deconvolution plots...")

  figures <- list()

  # 1. Cell type heatmap
  figures$heatmap <- plot_xcell2_heatmap(
    xcell2_results$scores,
    metadata,
    config,
    plots_dir
  )

  # 2. Stacked bar plot by condition
  figures$stacked_bar <- plot_celltype_stacked_bar(
    xcell2_results$scores,
    metadata,
    config,
    plots_dir
  )

  # 3. Box plots for significant cell types
  if (!is.null(diff_comp_results)) {
    figures$diff_boxplots <- plot_differential_celltypes(
      xcell2_results$scores,
      diff_comp_results,
      metadata,
      config,
      plots_dir
    )
  }

  # 4. Phenotype correlation heatmap
  if (!is.null(phenotype_corr)) {
    figures$phenotype_corr <- plot_phenotype_correlations(
      phenotype_corr,
      plots_dir
    )
  }

  # 5. Cell type composition PCA
  figures$pca <- plot_celltype_pca(
    xcell2_results$scores,
    metadata,
    config,
    plots_dir
  )

  return(figures)
}

#' Plot xCell2 Heatmap
plot_xcell2_heatmap <- function(scores, metadata, config, plots_dir) {
  log_message("  Creating cell type heatmap...")

  # Prepare matrix
  score_cols <- setdiff(colnames(scores), "sample_id")
  score_mat <- as.matrix(scores[, score_cols])
  rownames(score_mat) <- scores$sample_id

  # Get condition for annotation
  condition_col <- config$design$condition_column %||% "condition"

  if (condition_col %in% colnames(metadata)) {
    # Create annotation
    sample_order <- match(scores$sample_id, metadata$sample_id)
    ha <- ComplexHeatmap::HeatmapAnnotation(
      Condition = metadata[[condition_col]][sample_order],
      col = list(Condition = get_condition_colors(metadata[[condition_col]]))
    )
  } else {
    ha <- NULL
  }

  # Scale for visualization
  score_mat_scaled <- t(scale(t(score_mat)))

  # Create heatmap
  ht <- ComplexHeatmap::Heatmap(
    t(score_mat_scaled),
    name = "z-score",
    top_annotation = ha,
    show_column_names = ncol(score_mat_scaled) <= 50,
    show_row_names = TRUE,
    column_title = "Cell Type Enrichment Scores (xCell 2.0)",
    row_names_gp = grid::gpar(fontsize = 8),
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = "euclidean",
    clustering_method_rows = "ward.D2",
    clustering_method_columns = "ward.D2",
    col = circlize::colorRamp2(
      c(-2, 0, 2),
      c("blue", "white", "red")
    )
  )

  # Save
  fig_path <- file.path(plots_dir, "xcell2_heatmap.png")
  png(fig_path, width = 12, height = 10, units = "in", res = 150)
  ComplexHeatmap::draw(ht)
  dev.off()

  return(fig_path)
}

#' Plot Cell Type Stacked Bar
plot_celltype_stacked_bar <- function(scores, metadata, config, plots_dir) {
  log_message("  Creating stacked bar plot...")

  condition_col <- config$design$condition_column %||% "condition"

  # Prepare data
  score_cols <- setdiff(colnames(scores), "sample_id")

  scores_long <- scores %>%
    tidyr::pivot_longer(
      cols = all_of(score_cols),
      names_to = "cell_type",
      values_to = "score"
    )

  # Add condition
  if (condition_col %in% colnames(metadata)) {
    scores_long <- merge(
      scores_long,
      metadata[, c("sample_id", condition_col)],
      by.x = "sample_id",
      by.y = "sample_id"
    )

    # Calculate mean per condition and cell type
    summary_df <- scores_long %>%
      group_by(!!sym(condition_col), cell_type) %>%
      summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop")

    # Normalize to proportions within each condition
    summary_df <- summary_df %>%
      group_by(!!sym(condition_col)) %>%
      mutate(proportion = mean_score / sum(mean_score, na.rm = TRUE)) %>%
      ungroup()

    # Plot
    p <- ggplot(summary_df, aes(x = !!sym(condition_col), y = proportion, fill = cell_type)) +
      geom_bar(stat = "identity", position = "stack") +
      labs(
        title = "Cell Type Composition by Condition",
        x = "Condition",
        y = "Proportion",
        fill = "Cell Type"
      ) +
      theme_minimal() +
      theme(
        legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1)
      )

  } else {
    # Single bar for all samples
    summary_df <- scores_long %>%
      group_by(cell_type) %>%
      summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
      mutate(proportion = mean_score / sum(mean_score, na.rm = TRUE))

    p <- ggplot(summary_df, aes(x = "All Samples", y = proportion, fill = cell_type)) +
      geom_bar(stat = "identity", position = "stack") +
      labs(
        title = "Cell Type Composition",
        x = "",
        y = "Proportion",
        fill = "Cell Type"
      ) +
      theme_minimal()
  }

  fig_path <- file.path(plots_dir, "celltype_stacked_bar.png")
  ggsave(fig_path, p, width = 10, height = 8)

  return(fig_path)
}

#' Plot Differential Cell Types
plot_differential_celltypes <- function(scores, diff_results, metadata, config, plots_dir) {
  log_message("  Creating differential cell type plots...")

  condition_col <- config$design$condition_column %||% "condition"

  # Get significant cell types
  sig_types <- diff_results$cell_type[diff_results$padj < 0.05]

  if (length(sig_types) == 0) {
    log_message("  No significant cell types to plot.")
    return(NULL)
  }

  # Limit to top 12
  if (length(sig_types) > 12) {
    sig_types <- sig_types[1:12]
  }

  # Prepare data
  scores_meta <- merge(
    scores,
    metadata[, c("sample_id", condition_col)],
    by.x = "sample_id",
    by.y = "sample_id"
  )

  # Create box plots
  plots <- lapply(sig_types, function(ct) {
    padj <- diff_results$padj[diff_results$cell_type == ct]

    ggplot(scores_meta, aes(x = !!sym(condition_col), y = !!sym(ct), fill = !!sym(condition_col))) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.5) +
      labs(
        title = paste0(ct, "\n(adj.p = ", format(padj, digits = 2), ")"),
        x = "",
        y = "xCell2 Score"
      ) +
      theme_minimal() +
      theme(legend.position = "none")
  })

  # Combine plots
  combined <- patchwork::wrap_plots(plots, ncol = 3)

  fig_path <- file.path(plots_dir, "differential_celltypes.png")
  ggsave(fig_path, combined, width = 12, height = 4 * ceiling(length(sig_types) / 3))

  return(fig_path)
}

#' Plot Phenotype Correlations
plot_phenotype_correlations <- function(phenotype_corr, plots_dir) {
  log_message("  Creating phenotype correlation heatmap...")

  # Pivot to matrix
  corr_mat <- phenotype_corr %>%
    select(cell_type, phenotype, correlation) %>%
    tidyr::pivot_wider(
      names_from = phenotype,
      values_from = correlation
    ) %>%
    tibble::column_to_rownames("cell_type") %>%
    as.matrix()

  # Get significance matrix
  sig_mat <- phenotype_corr %>%
    select(cell_type, phenotype, padj) %>%
    tidyr::pivot_wider(
      names_from = phenotype,
      values_from = padj
    ) %>%
    tibble::column_to_rownames("cell_type") %>%
    as.matrix()

  # Create annotation for significance
  sig_fun <- function(j, i, x, y, w, h, fill) {
    if (!is.na(sig_mat[i, j]) && sig_mat[i, j] < 0.05) {
      grid::grid.text("*", x, y, gp = grid::gpar(fontsize = 14))
    }
  }

  # Create heatmap
  ht <- ComplexHeatmap::Heatmap(
    corr_mat,
    name = "Correlation",
    cell_fun = sig_fun,
    column_title = "Cell Type - Phenotype Correlations",
    row_names_gp = grid::gpar(fontsize = 9),
    column_names_gp = grid::gpar(fontsize = 10),
    col = circlize::colorRamp2(
      c(-1, 0, 1),
      c("blue", "white", "red")
    ),
    na_col = "grey90"
  )

  fig_path <- file.path(plots_dir, "phenotype_correlations.png")
  png(fig_path, width = 10, height = 8, units = "in", res = 150)
  ComplexHeatmap::draw(ht)
  dev.off()

  return(fig_path)
}

#' Plot Cell Type PCA
plot_celltype_pca <- function(scores, metadata, config, plots_dir) {
  log_message("  Creating cell type PCA...")

  condition_col <- config$design$condition_column %||% "condition"

  # Prepare matrix
  score_cols <- setdiff(colnames(scores), "sample_id")
  score_mat <- as.matrix(scores[, score_cols])
  rownames(score_mat) <- scores$sample_id

  # Run PCA
  pca_result <- prcomp(score_mat, scale. = TRUE, center = TRUE)

  # Create plot data
  pca_df <- as.data.frame(pca_result$x[, 1:min(2, ncol(pca_result$x))])
  pca_df$sample_id <- rownames(pca_df)

  # Add metadata
  if (condition_col %in% colnames(metadata)) {
    pca_df <- merge(pca_df, metadata[, c("sample_id", condition_col)], by = "sample_id")

    # Variance explained
    var_exp <- summary(pca_result)$importance[2, 1:2] * 100

    p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = !!sym(condition_col))) +
      geom_point(size = 3) +
      labs(
        title = "PCA of Cell Type Composition",
        x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
        y = paste0("PC2 (", round(var_exp[2], 1), "%)")
      ) +
      theme_minimal()
  } else {
    var_exp <- summary(pca_result)$importance[2, 1:2] * 100

    p <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
      geom_point(size = 3) +
      labs(
        title = "PCA of Cell Type Composition",
        x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
        y = paste0("PC2 (", round(var_exp[2], 1), "%)")
      ) +
      theme_minimal()
  }

  fig_path <- file.path(plots_dir, "celltype_pca.png")
  ggsave(fig_path, p, width = 8, height = 6)

  return(fig_path)
}

# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------

#' Get Expression Matrix from DESeqDataSet
get_expression_matrix <- function(dds, config) {
  # Use normalized counts or TPM if available
  if ("tpm" %in% SummarizedExperiment::assayNames(dds)) {
    expr <- SummarizedExperiment::assay(dds, "tpm")
  } else if ("normalized" %in% SummarizedExperiment::assayNames(dds)) {
    expr <- SummarizedExperiment::assay(dds, "normalized")
  } else {
    # Calculate normalized counts
    expr <- DESeq2::counts(dds, normalized = TRUE)
  }

  # Convert to matrix with gene symbols as rownames
  expr_mat <- as.matrix(expr)

  # If rownames are Ensembl IDs, try to convert to symbols
  if (grepl("^ENS", rownames(expr_mat)[1])) {
    expr_mat <- convert_ensembl_to_symbol(expr_mat, config)
  }

  return(expr_mat)
}

#' Convert Ensembl IDs to Gene Symbols
convert_ensembl_to_symbol <- function(expr_mat, config) {
  log_message("  Converting Ensembl IDs to gene symbols...")

  organism <- config$deconvolution$organism %||% "human"

  org_db <- switch(
    organism,
    "human" = "org.Hs.eg.db",
    "mouse" = "org.Mm.eg.db",
    "org.Hs.eg.db"  # default
  )

  if (!requireNamespace(org_db, quietly = TRUE)) {
    log_message("WARNING: ", org_db, " not available. Keeping Ensembl IDs.")
    return(expr_mat)
  }

  tryCatch({
    # Strip version numbers
    ensembl_ids <- sub("\\..*", "", rownames(expr_mat))

    # Get symbols
    symbols <- AnnotationDbi::mapIds(
      get(org_db),
      keys = ensembl_ids,
      keytype = "ENSEMBL",
      column = "SYMBOL"
    )

    # Replace rownames
    valid_symbols <- !is.na(symbols)
    rownames(expr_mat)[valid_symbols] <- symbols[valid_symbols]

    # Remove duplicates
    expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), ]

    return(expr_mat)

  }, error = function(e) {
    log_message("WARNING: Symbol conversion failed: ", e$message)
    return(expr_mat)
  })
}

#' Get Condition Colors
get_condition_colors <- function(conditions) {
  unique_conds <- unique(conditions)
  n_conds <- length(unique_conds)

  colors <- RColorBrewer::brewer.pal(max(3, n_conds), "Set1")[1:n_conds]
  names(colors) <- unique_conds

  return(colors)
}

#' Create Deconvolution Summary
create_deconvolution_summary <- function(xcell2_results, diff_comp_results, phenotype_corr) {
  summary <- list(
    n_cell_types = ncol(xcell2_results$scores) - 1,
    reference = xcell2_results$reference,
    organism = xcell2_results$organism
  )

  if (!is.null(diff_comp_results)) {
    summary$n_differential <- sum(diff_comp_results$padj < 0.05, na.rm = TRUE)
    summary$top_differential <- head(diff_comp_results$cell_type[diff_comp_results$padj < 0.05], 5)
  }

  if (!is.null(phenotype_corr)) {
    summary$n_phenotype_correlations <- sum(phenotype_corr$padj < 0.05, na.rm = TRUE)
  }

  return(summary)
}

#' Save Deconvolution Outputs
save_deconvolution_outputs <- function(results, output_dir) {
  log_message("Saving deconvolution outputs...")

  # Save scores
  write.csv(
    results$scores,
    file.path(output_dir, "xcell2_scores.csv"),
    row.names = FALSE
  )

  # Save raw scores (unfiltered)
  write.csv(
    results$scores_raw,
    file.path(output_dir, "xcell2_scores_raw.csv"),
    row.names = FALSE
  )

  # Save full results as RDS
  saveRDS(results, file.path(output_dir, "deconvolution_results.rds"))

  log_message("  Outputs saved to: ", output_dir)
}

# -----------------------------------------------------------------------------
# Legacy xCell Fallback
# -----------------------------------------------------------------------------

#' Run Legacy xCell Analysis (Fallback)
run_xcell_legacy <- function(dds, de_results, config, output_dir, plots_dir) {
  log_message("Running legacy xCell analysis (xCell 2.0 not available)...")

  if (!requireNamespace("xCell", quietly = TRUE)) {
    log_message("ERROR: xCell package not available.")
    return(NULL)
  }

  # Get expression matrix
  expr_matrix <- get_expression_matrix(dds, config)

  tryCatch({
    # Run xCell
    scores <- xCell::xCellAnalysis(expr_matrix)

    # Convert to data frame
    scores_df <- as.data.frame(t(scores))
    scores_df$sample_id <- rownames(scores_df)

    # Save
    write.csv(scores_df, file.path(output_dir, "xcell_scores.csv"), row.names = FALSE)

    return(list(
      scores = scores_df,
      scores_raw = scores_df,
      reference = "xCell_default",
      organism = "human",
      differential_composition = NULL,
      phenotype_correlations = NULL,
      figures = list(),
      summary = list(n_cell_types = ncol(scores_df) - 1, note = "Legacy xCell used")
    ))

  }, error = function(e) {
    log_message("ERROR in legacy xCell: ", e$message)
    return(NULL)
  })
}

# -----------------------------------------------------------------------------
# Organism Helper Functions
# -----------------------------------------------------------------------------

#' Normalize organism name to standard format
#' @param organism Character string with organism name
#' @return Standardized organism name ("human", "mouse", or original)
normalize_organism <- function(organism) {
  if (is.null(organism)) return("human")

  organism_lower <- tolower(organism)

  # Map various organism names to standardized forms
  if (organism_lower %in% c("human", "homo sapiens", "homo_sapiens", "hs", "hsapiens")) {
    return("human")
  }

  if (organism_lower %in% c("mouse", "mus musculus", "mus_musculus", "mm", "mmusculus")) {
    return("mouse")
  }

  if (organism_lower %in% c("rat", "rattus norvegicus", "rattus_norvegicus", "rn", "rnorvegicus")) {
    return("rat")
  }

  # Return original if no match
  return(organism)
}

#' Get default xCell2 reference for organism
#' @param organism Normalized organism name
#' @return Default reference name
get_default_reference <- function(organism) {
  switch(organism,
    "human" = "BlueprintEncode",
    "mouse" = "ImmGen",
    "BlueprintEncode"  # Default fallback
  )
}
