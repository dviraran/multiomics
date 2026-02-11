# =============================================================================
# 07_enrichment.R — MetaboAnalystR Quantitative Enrichment Analysis (QEA)
# Uses globalTest on full data matrix (no need for significant features)
# Supports KEGG compound and SMPDB pathway libraries
# =============================================================================

run_enrichment <- function(de_results, normalized_data, ingested_data, config) {
  log_message("=== Starting Enrichment Analysis (QEA) ===")

  if (!isTRUE(config$enrichment$run_enrichment)) {
    log_message("Enrichment analysis disabled in config")
    return(NULL)
  }

  if (!requireNamespace("MetaboAnalystR", quietly = TRUE)) {
    warning("MetaboAnalystR not installed. Skipping enrichment.")
    return(NULL)
  }

  # Libraries to run (default: both SMPDB and KEGG)
  libraries <- config$enrichment$libraries
  if (is.null(libraries)) {
    lib_single <- config$enrichment$library
    if (!is.null(lib_single)) {
      libraries <- list(lib_single)
    } else {
      libraries <- list("smpdb_pathway", "kegg_pathway")
    }
  }
  gmt_files <- config$enrichment$gmt_file
  if (!is.null(gmt_files)) {
    # Normalise to a character vector (YAML may parse as list or vector)
    gmt_files <- unlist(gmt_files)
    for (gf in gmt_files) {
      if (!file.exists(gf)) {
        log_message("Warning: GMT file not found: ", gf)
      } else {
        # Derive a library label from filename (e.g. "kegg_legionella" / "smpdb")
        gmt_label <- tools::file_path_sans_ext(basename(gf))
        log_message("Using GMT file for enrichment: ", gf, " (label: ", gmt_label, ")")
        libraries <- c(libraries, gmt_label)
      }
    }
    # Store as named vector for lookup later
    names(gmt_files) <- tools::file_path_sans_ext(basename(gmt_files))
  } else {
    gmt_files <- character(0)
  }

  # Prepare compound data file for MetaboAnalystR
  data_file <- prepare_qea_data(ingested_data, config)
  if (is.null(data_file)) {
    return(NULL)
  }

  # Run QEA for each library
  all_results <- list()
  for (lib in libraries) {
    if (lib %in% names(gmt_files)) {
      # Custom GMT file
      gf <- gmt_files[[lib]]
      log_message("Running globalTest with GMT: ", lib, " (", gf, ")")
      result <- tryCatch(
        run_qea_gmt(data_file, gf, config),
        error = function(e) {
          log_message("QEA failed for GMT ", lib, ": ", e$message)
          NULL
        }
      )
    } else {
      log_message("Running QEA globalTest with library: ", lib)
      result <- tryCatch(
        run_qea_globaltest(data_file, lib, config),
        error = function(e) {
          log_message("QEA failed for ", lib, ": ", e$message)
          NULL
        }
      )
    }
    if (!is.null(result)) {
      all_results[[lib]] <- result
    }
  }

  # Clean up temp dir
  unlink(dirname(data_file), recursive = TRUE)

  if (length(all_results) == 0) {
    log_message("No enrichment results from any library")
    return(NULL)
  }

  # Combine results from all libraries
  combined_df <- do.call(rbind, lapply(names(all_results), function(lib) {
    df <- all_results[[lib]]
    df$library <- lib
    df
  }))

  # FDR correction across all pathways
  combined_df$FDR <- p.adjust(combined_df$`Raw p`, method = "fdr")
  combined_df <- combined_df[order(combined_df$FDR), ]

  save_table(combined_df, "enrichment_results.csv", config, "tables")

  # Combined plot (all libraries together)
  plot <- plot_enrichment_barplot(combined_df, config)

  # Per-library plots (separate KEGG and SMPDB bar charts)
  per_lib_plots <- list()
  for (lib_name in unique(combined_df$library)) {
    lib_df <- combined_df[combined_df$library == lib_name, ]
    lib_df <- lib_df[order(lib_df$FDR), ]
    lib_label <- gsub("_pathway$", "", lib_name)
    lib_label <- toupper(lib_label)
    p_lib <- plot_enrichment_barplot_single(lib_df, config, lib_label)
    per_lib_plots[[lib_name]] <- p_lib
  }

  n_sig <- sum(combined_df$FDR < 0.05, na.rm = TRUE)
  log_message(
    "Enrichment complete: ", nrow(combined_df), " pathways tested, ",
    n_sig, " with FDR < 0.05"
  )
  log_message("=== Enrichment Analysis Complete ===")

  list(
    table = combined_df,
    plot = plot,
    per_library_plots = per_lib_plots,
    method = "globaltest_qea"
  )
}

# ---------- Prepare data file for MetaboAnalystR QEA -------------------------

prepare_qea_data <- function(ingested_data, config) {
  mat_raw <- ingested_data$matrix
  metadata <- ingested_data$metadata
  fmeta <- ingested_data$feature_metadata

  sample_col <- config$input$sample_id_column
  condition_col <- config$design$condition_column

  # Extract compound names for MetaboAnalystR mapping
  if ("Molecule" %in% colnames(fmeta)) {
    compound_names <- fmeta$Molecule
  } else {
    # Try to extract from feature_id (HMDB|Molecule format)
    compound_names <- sub("^[^|]+\\|", "", rownames(mat_raw))
  }

  # --- HMDB to KEGG Mapping (if provided) ---
  mapping_file <- config$enrichment$mapping_file
  if (!is.null(mapping_file) && file.exists(mapping_file)) {
    log_message("Using external mapping file: ", mapping_file)
    mapping_df <- readr::read_delim(mapping_file, delim = "\t", col_types = readr::cols())

    # Expected columns: Molecule, KEGG (adjust based on actual file)
    # If the file has no headers, we might need to adjust. Assuming headers for now.
    # User's file likely matches Molecule or HMDB to KEGG.

    # Try to map based on row names or a specific ID column
    # The input matrix rownames are likely HMDB IDs or Molecule names
    # Let's assume rownames are HMDB IDs or can be mapped from mapping_file

    # Check if we have HMDB IDs in feature_metadata
    hmdb_ids <- NULL
    if ("HMDB" %in% colnames(fmeta)) {
      hmdb_ids <- fmeta$HMDB
    } else {
      # Fallback: assume rownames are HMDB IDs
      hmdb_ids <- rownames(mat_raw)
    }

    # Clean HMDB IDs (remove 'HMDB' prefix if needed, though usually standard)
    # Mapping file format check:
    # If mapping_file is 'HMDB_to_KEGG.level1.txt', it likely has HMDB and KEGG columns.

    if ("HMDB" %in% colnames(mapping_df) && "KEGG" %in% colnames(mapping_df)) {
      # simple mapping
      map_vec <- setNames(mapping_df$KEGG, mapping_df$HMDB)
      mapped_kegg <- map_vec[hmdb_ids]

      # Replace compound_names with KEGG IDs where mapping exists
      has_map <- !is.na(mapped_kegg) & mapped_kegg != ""
      compound_names[has_map] <- mapped_kegg[has_map]

      log_message("Mapped ", sum(has_map), " features to KEGG IDs using external file")
    } else {
      log_message("Mapping file columns not recognized. Expected 'HMDB' and 'KEGG'.")
      # Attempt flexible column matching if needed, or log warning
      print(head(mapping_df))
    }
  }

  # Filter: keep only features with valid compound names
  valid <- !is.na(compound_names) & compound_names != "" & compound_names != "NA"
  mat_use <- mat_raw[valid, , drop = FALSE]
  names_use <- compound_names[valid]

  # Deduplicate compound names (keep first occurrence)
  dup <- duplicated(names_use)
  mat_use <- mat_use[!dup, , drop = FALSE]
  rownames(mat_use) <- names_use[!dup]

  log_message(
    "Prepared ", nrow(mat_use), " compounds for enrichment (",
    sum(valid) - sum(!dup[valid]), " duplicates removed)"
  )

  # Build MetaboAnalystR format: rows = samples, cols = compounds
  # First two columns must be Sample and Group
  df_t <- as.data.frame(t(mat_use), check.names = FALSE)
  conditions <- metadata[[condition_col]][
    match(rownames(df_t), metadata[[sample_col]])
  ]
  df_t <- cbind(
    data.frame(
      Sample = rownames(df_t),
      Group = as.character(conditions),
      stringsAsFactors = FALSE
    ),
    df_t
  )

  # Write to temp file
  tmp_dir <- tempfile("metabo_qea_")
  dir.create(tmp_dir, recursive = TRUE)
  out_path <- file.path(tmp_dir, "combined_data.txt")
  write.table(df_t, out_path, sep = "\t", row.names = FALSE, quote = FALSE)

  log_message(
    "Wrote QEA data file: ", nrow(df_t), " samples x ",
    ncol(df_t) - 2, " compounds"
  )
  out_path
}

# ---------- Run QEA globalTest for a single library ---------------------------

run_qea_globaltest <- function(data_file, library_name, config) {
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)

  row_norm <- config$normalization$row_norm %||% "NULL"
  trans_norm <- config$normalization$trans_norm %||% "LogNorm"
  scale_norm <- config$normalization$scale_norm %||% "MeanCenter"

  mSet <- MetaboAnalystR::InitDataObjects("conc", "msetqea", FALSE)
  mSet <- MetaboAnalystR::Read.TextData(mSet, data_file, "rowu", "disc")
  mSet <- MetaboAnalystR::SanityCheckData(mSet)
  mSet <- MetaboAnalystR::ReplaceMin(mSet)
  mSet <- MetaboAnalystR::PreparePrenormData(mSet)
  mSet <- MetaboAnalystR::Normalization(mSet, row_norm, trans_norm, scale_norm,
    "S10T0",
    ratio = FALSE, ratioNum = 20
  )

  # Cross-reference compound names
  mSet <- MetaboAnalystR::CrossReferencing(mSet, "name")
  mSet <- MetaboAnalystR::CreateMappingResultTable(mSet)

  # Log mapping stats
  map_tbl <- mSet$dataSet$map.table
  if (!is.null(map_tbl)) {
    n_mapped <- sum(map_tbl[, "Match"] != "", na.rm = TRUE)
    log_message("Compound mapping: ", n_mapped, " / ", nrow(map_tbl), " mapped")
  }

  # Set pathway library and run globalTest
  mSet <- MetaboAnalystR::SetMetabolomeFilter(mSet, FALSE)
  mSet <- MetaboAnalystR::SetCurrentMsetLib(mSet, library_name, 2)
  mSet <- MetaboAnalystR::CalculateGlobalTestScore(mSet)

  # Extract QEA results
  qea_mat <- mSet$analSet$qea.mat
  if (is.null(qea_mat) || nrow(qea_mat) == 0) {
    log_message("No results for library: ", library_name)
    return(NULL)
  }

  enrich_df <- as.data.frame(qea_mat)
  enrich_df$pathway <- rownames(enrich_df)
  enrich_df <- enrich_df[order(enrich_df$`Raw p`), ]

  n_sig <- sum(enrich_df$`Raw p` < 0.05, na.rm = TRUE)
  log_message(
    library_name, ": ", nrow(enrich_df), " pathways, ",
    n_sig, " significant (p < 0.05)"
  )
  if (nrow(enrich_df) > 0) {
    log_message(
      "  Top: ", enrich_df$pathway[1],
      " (p = ", signif(enrich_df$`Raw p`[1], 3), ")"
    )
  }

  enrich_df
}

# ---------- Enrichment barplot ------------------------------------------------

plot_enrichment_barplot <- function(enrich_df, config, top_n = 20) {
  log_message("Generating enrichment barplot")

  top_df <- head(enrich_df, top_n)
  top_df$neg_log10_fdr <- -log10(pmax(top_df$FDR, 1e-20))
  top_df$pathway_short <- ifelse(
    nchar(top_df$pathway) > 50,
    paste0(substr(top_df$pathway, 1, 47), "..."),
    top_df$pathway
  )

  # Color by library if present
  if ("library" %in% colnames(top_df)) {
    top_df$lib_label <- gsub("_pathway", "", top_df$library)
    top_df$lib_label <- toupper(top_df$lib_label)
  }

  top_df$pathway_short <- factor(
    top_df$pathway_short,
    levels = rev(top_df$pathway_short)
  )

  if ("library" %in% colnames(top_df)) {
    p <- ggplot2::ggplot(
      top_df,
      ggplot2::aes(x = pathway_short, y = neg_log10_fdr, fill = lib_label)
    ) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::labs(fill = "Database")
  } else {
    p <- ggplot2::ggplot(
      top_df,
      ggplot2::aes(x = pathway_short, y = neg_log10_fdr, fill = neg_log10_fdr)
    ) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_gradient(
        low = "steelblue", high = "firebrick",
        name = "-log10(FDR)"
      )
  }

  p <- p +
    ggplot2::geom_hline(
      yintercept = -log10(0.05), linetype = "dashed",
      color = "grey40"
    ) +
    ggplot2::labs(
      title = "Pathway Enrichment (GlobalTest QEA)",
      subtitle = paste0("Top ", nrow(top_df), " pathways"),
      x = NULL,
      y = "-log10(FDR)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 9),
      plot.title = ggplot2::element_text(face = "bold")
    )

  save_plot(p, "enrichment_barplot.png", config,
    width = 12, height = 8, subdir = "plots"
  )

  p
}

# ---------- Single-library enrichment barplot ---------------------------------

plot_enrichment_barplot_single <- function(enrich_df, config, lib_label,
                                           top_n = 20) {
  top_df <- head(enrich_df, top_n)
  if (nrow(top_df) == 0) return(NULL)

  top_df$neg_log10_fdr <- -log10(pmax(top_df$FDR, 1e-20))
  top_df$pathway_short <- ifelse(
    nchar(top_df$pathway) > 50,
    paste0(substr(top_df$pathway, 1, 47), "..."),
    top_df$pathway
  )
  top_df$pathway_short <- factor(
    top_df$pathway_short,
    levels = rev(top_df$pathway_short)
  )

  p <- ggplot2::ggplot(
    top_df,
    ggplot2::aes(x = pathway_short, y = neg_log10_fdr, fill = neg_log10_fdr)
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_gradient(
      low = "steelblue", high = "firebrick",
      name = "-log10(FDR)"
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(0.05), linetype = "dashed", color = "grey40"
    ) +
    ggplot2::labs(
      title = paste0(lib_label, " Pathway Enrichment"),
      subtitle = paste0("Top ", nrow(top_df), " pathways (GlobalTest QEA)"),
      x = NULL,
      y = "-log10(FDR)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 9),
      plot.title = ggplot2::element_text(face = "bold")
    )

  fname <- paste0("enrichment_barplot_", tolower(lib_label), ".png")
  save_plot(p, fname, config, width = 12, height = 8, subdir = "plots")

  p
}

# ---------- Run QEA with GMT file ---------------------------------------------
run_qea_gmt <- function(data_file, gmt_file, config) {
  # 1. Parse GMT (with descriptions for labeling)
  gmt_parsed <- read_gmt_list(gmt_file, include_descriptions = TRUE)
  gmt_list <- gmt_parsed$sets
  desc_map <- gmt_parsed$descriptions
  if (length(gmt_list) == 0) {
    return(NULL)
  }

  log_message("Loaded ", length(gmt_list), " pathways from GMT")

  # Translate HMDB IDs in GMT to KEGG IDs using mapping file (if available),
  # so that HMDB-based GMTs (e.g. SMPDB) match the KEGG-mapped data columns.
  mapping_file <- config$enrichment$mapping_file
  if (!is.null(mapping_file) && file.exists(mapping_file)) {
    mapping_df <- readr::read_delim(mapping_file, delim = "\t",
                                    col_types = readr::cols(),
                                    show_col_types = FALSE)
    if ("HMDB" %in% colnames(mapping_df) && "KEGG" %in% colnames(mapping_df)) {
      hmdb2kegg <- setNames(mapping_df$KEGG, mapping_df$HMDB)
      gmt_list <- lapply(gmt_list, function(cpds) {
        mapped <- hmdb2kegg[cpds]
        ifelse(!is.na(mapped) & mapped != "", mapped, cpds)
      })
      log_message("Translated HMDB IDs in GMT to KEGG IDs using mapping file")
    }
  }

  # Rename gene set keys to "ID - Description" where descriptions exist
  labeled_names <- make_pathway_labels(names(gmt_list), desc_map)
  names(gmt_list) <- labeled_names

  if (!requireNamespace("globaltest", quietly = TRUE)) {
    log_message("globaltest package required for GMT enrichment but not found.")
    return(NULL)
  }

  # Load data
  df <- read.table(data_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

  # MetaboAnalystR format: Sample, Group, Compound1, Compound2...
  # globaltest expects: response variable (Group), and data matrix (Compounds)

  response <- factor(df$Group)
  # Data matrix: rows = samples, cols = compounds
  # df has 2 metadata columns
  X <- as.matrix(df[, -c(1, 2), drop = FALSE])

  # GlobalTest needs features in columns?
  # gt(response, alternative)
  # If alternative is matrix, rows should be samples. Matches X structure.

  # 3. Run globaltest for each pathway in GMT
  # We need to match GMT compound IDs to X colnames
  # X colnames are mapped KEGG IDs or Names from prepare_qea_data

  # Filter GMT subsets to only include compounds present in the data matrix.
  # globaltest crashes with "subscript out of bounds" if subset elements
  # reference columns that don't exist in X.
  available <- colnames(X)
  subsets <- lapply(gmt_list, function(cpds) cpds[cpds %in% available])
  # Remove pathways with fewer than 2 matching compounds
  keep <- vapply(subsets, length, integer(1)) >= 2L
  subsets <- subsets[keep]

  if (length(subsets) == 0) {
    log_message("No GMT pathways have >= 2 matching compounds in the data")
    return(NULL)
  }
  log_message(
    "Retained ", length(subsets), " / ", length(gmt_list),
    " pathways with >= 2 matching compounds"
  )

  res_gt <- tryCatch(
    globaltest::gt(response, X, subsets = subsets),
    error = function(e) {
      log_message("globaltest error: ", e$message)
      NULL
    }
  )

  if (is.null(res_gt)) {
    return(NULL)
  }

  # Extract results table
  # globaltest result object printing gives summary, result(res_gt) gives table
  res_tbl <- globaltest::result(res_gt)

  # Columns: p-value, Statistic, Expected, Std.Dev, #Cov
  # We need to format like MetaboAnalystR output for consistency
  # MetaboAnalystR QEA out: "P-value", "Holm adjusted", "FDR", "Hits", ... (variable names differ)
  # Our 'combine' step expects 'Raw p' column.

  out_df <- data.frame(
    pathway = rownames(res_tbl),
    `Raw p` = res_tbl[, "p-value"],
    Hits = res_tbl[, "#Cov"], # approximate match to hits
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  colnames(out_df)[colnames(out_df) == "Raw p"] <- "Raw p"

  # Add library column if needed (done in main loop)

  return(out_df)
}

# ---------- Make pathway labels from IDs + descriptions -----------------------
# Produces "ID - Description" when a description exists, otherwise keeps ID.

make_pathway_labels <- function(ids, desc_map) {
  vapply(ids, function(id) {
    desc <- desc_map[[id]]
    if (!is.null(desc) && nzchar(desc)) {
      paste0(id, " - ", desc)
    } else {
      id
    }
  }, character(1), USE.NAMES = FALSE)
}

read_gmt_list <- function(gmt_file, include_descriptions = FALSE) {
  if (!file.exists(gmt_file)) {
    if (include_descriptions) return(list(sets = list(), descriptions = character(0)))
    return(list())
  }

  lines <- readLines(gmt_file)
  gmt_list <- list()
  desc_map <- character(0)

  for (line in lines) {
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) < 3) next # ID, Desc, Gene1...

    pathway_id <- parts[1]
    pathway_desc <- parts[2]
    genes <- parts[3:length(parts)]

    # Remove empty strings
    genes <- genes[genes != ""]

    if (length(genes) > 0) {
      gmt_list[[pathway_id]] <- genes
      if (nzchar(pathway_desc)) {
        desc_map[[pathway_id]] <- pathway_desc
      }
    }
  }

  if (include_descriptions) {
    return(list(sets = gmt_list, descriptions = desc_map))
  }
  return(gmt_list)
}

# =============================================================================
# ssGSEA Pathway Enrichment (GSVA + Wilcoxon rank-sum)
# =============================================================================

run_ssgsea_enrichment <- function(normalized_data, ingested_data, config) {
  log_message("=== Starting ssGSEA Pathway Enrichment ===")

  if (!isTRUE(config$enrichment$run_enrichment)) {
    log_message("Enrichment analysis disabled in config")
    return(NULL)
  }

  if (!requireNamespace("GSVA", quietly = TRUE)) {
    warning("GSVA not installed. Skipping ssGSEA enrichment.")
    return(NULL)
  }

  # --- Build expression matrix (features x samples) with compound IDs --------
  expr_mat <- prepare_ssgsea_data(ingested_data, config)
  if (is.null(expr_mat) || nrow(expr_mat) < 2) {
    log_message("Too few features for ssGSEA")
    return(NULL)
  }

  # --- Build gene sets from GMT files and/or built-in libraries ---------------
  gene_sets <- build_ssgsea_gene_sets(config)
  if (length(gene_sets) == 0) {
    log_message("No gene sets available for ssGSEA")
    return(NULL)
  }

  # Filter gene sets: keep only those with >= 2 members present in expr_mat
  available_ids <- rownames(expr_mat)
  gene_sets <- lapply(gene_sets, function(cpds) cpds[cpds %in% available_ids])
  set_sizes <- vapply(gene_sets, length, integer(1))
  gene_sets <- gene_sets[set_sizes >= 2L]

  if (length(gene_sets) == 0) {
    log_message("No gene sets have >= 2 matching compounds in data")
    return(NULL)
  }
  log_message(
    "ssGSEA: ", length(gene_sets), " gene sets with >= 2 members in data"
  )

  # --- Run ssGSEA via GSVA ---------------------------------------------------
  scores <- tryCatch({
    # Try new API (GSVA >= 1.46)
    param <- GSVA::ssgseaParam(expr_mat, gene_sets)
    GSVA::gsva(param)
  }, error = function(e1) {
    tryCatch({
      # Fallback: old API
      GSVA::gsva(expr_mat, gene_sets, method = "ssgsea", verbose = FALSE)
    }, error = function(e2) {
      log_message("ssGSEA failed: ", e2$message)
      NULL
    })
  })

  if (is.null(scores) || nrow(scores) == 0) {
    log_message("ssGSEA produced no results")
    return(NULL)
  }
  log_message("ssGSEA computed scores for ", nrow(scores), " pathways x ",
              ncol(scores), " samples")

  # --- Wilcoxon test per pathway ----------------------------------------------
  metadata <- ingested_data$metadata
  sample_col <- config$input$sample_id_column
  condition_col <- config$design$condition_column

  conditions <- metadata[[condition_col]][
    match(colnames(scores), metadata[[sample_col]])
  ]
  conditions <- factor(conditions)
  cond_levels <- levels(conditions)

  if (length(cond_levels) != 2) {
    log_message("ssGSEA Wilcoxon test requires exactly 2 conditions, found: ",
                length(cond_levels))
    return(NULL)
  }

  grp1_idx <- which(conditions == cond_levels[1])
  grp2_idx <- which(conditions == cond_levels[2])

  pvalues <- vapply(seq_len(nrow(scores)), function(i) {
    tryCatch(
      wilcox.test(scores[i, grp1_idx], scores[i, grp2_idx])$p.value,
      error = function(e) NA_real_
    )
  }, numeric(1))

  fdr <- p.adjust(pvalues, method = "fdr")

  mean1 <- rowMeans(scores[, grp1_idx, drop = FALSE], na.rm = TRUE)
  mean2 <- rowMeans(scores[, grp2_idx, drop = FALSE], na.rm = TRUE)

  results_df <- data.frame(
    pathway = rownames(scores),
    p_value = pvalues,
    FDR = fdr,
    mean_score_cond1 = mean1,
    mean_score_cond2 = mean2,
    score_diff = mean2 - mean1,
    significant = !is.na(fdr) & fdr < 0.05,
    stringsAsFactors = FALSE
  )
  colnames(results_df)[colnames(results_df) == "mean_score_cond1"] <-
    paste0("mean_", cond_levels[1])
  colnames(results_df)[colnames(results_df) == "mean_score_cond2"] <-
    paste0("mean_", cond_levels[2])

  results_df <- results_df[order(results_df$p_value), ]

  save_table(results_df, "ssgsea_results.csv", config, "tables")

  n_sig <- sum(results_df$significant, na.rm = TRUE)
  log_message("ssGSEA: ", nrow(results_df), " pathways tested, ",
              n_sig, " significant (FDR < 0.05)")

  # --- Plots ------------------------------------------------------------------
  barplot_gg <- plot_ssgsea_barplot(results_df, config)
  boxplot_gg <- plot_ssgsea_boxplots(scores, conditions, results_df, config)

  log_message("=== ssGSEA Pathway Enrichment Complete ===")

  list(
    table = results_df,
    scores = scores,
    barplot = barplot_gg,
    boxplots = boxplot_gg,
    method = "ssgsea_wilcoxon"
  )
}

# ---------- Prepare expression matrix for ssGSEA ------------------------------

prepare_ssgsea_data <- function(ingested_data, config) {
  mat_raw <- ingested_data$matrix
  fmeta <- ingested_data$feature_metadata

  # Extract compound names
  if ("Molecule" %in% colnames(fmeta)) {
    compound_names <- fmeta$Molecule
  } else {
    compound_names <- sub("^[^|]+\\|", "", rownames(mat_raw))
  }

  # HMDB to KEGG mapping (if provided)
  mapping_file <- config$enrichment$mapping_file
  if (!is.null(mapping_file) && file.exists(mapping_file)) {
    mapping_df <- readr::read_delim(mapping_file, delim = "\t",
                                    col_types = readr::cols(),
                                    show_col_types = FALSE)

    hmdb_ids <- if ("HMDB" %in% colnames(fmeta)) fmeta$HMDB else rownames(mat_raw)

    if ("HMDB" %in% colnames(mapping_df) && "KEGG" %in% colnames(mapping_df)) {
      map_vec <- setNames(mapping_df$KEGG, mapping_df$HMDB)
      mapped_kegg <- map_vec[hmdb_ids]
      has_map <- !is.na(mapped_kegg) & mapped_kegg != ""
      compound_names[has_map] <- mapped_kegg[has_map]
      log_message("ssGSEA: Mapped ", sum(has_map), " features to KEGG IDs")
    }
  }

  # Filter valid names
  valid <- !is.na(compound_names) & compound_names != "" & compound_names != "NA"
  mat_use <- mat_raw[valid, , drop = FALSE]
  names_use <- compound_names[valid]

  # Deduplicate
  dup <- duplicated(names_use)
  mat_use <- mat_use[!dup, , drop = FALSE]
  rownames(mat_use) <- names_use[!dup]

  log_message("ssGSEA: Prepared ", nrow(mat_use), " compounds x ",
              ncol(mat_use), " samples")

  as.matrix(mat_use)
}

# ---------- Build gene sets for ssGSEA ----------------------------------------

build_ssgsea_gene_sets <- function(config) {
  gene_sets <- list()

  # From custom GMT files
  gmt_files <- config$enrichment$gmt_file
  if (!is.null(gmt_files)) {
    gmt_files <- unlist(gmt_files)
    for (gf in gmt_files) {
      if (file.exists(gf)) {
        gmt_parsed <- read_gmt_list(gf, include_descriptions = TRUE)
        gmt <- gmt_parsed$sets
        desc_map <- gmt_parsed$descriptions

        # Translate HMDB IDs to KEGG if mapping file available
        mapping_file <- config$enrichment$mapping_file
        if (!is.null(mapping_file) && file.exists(mapping_file)) {
          mapping_df <- readr::read_delim(mapping_file, delim = "\t",
                                          col_types = readr::cols(),
                                          show_col_types = FALSE)
          if ("HMDB" %in% colnames(mapping_df) && "KEGG" %in% colnames(mapping_df)) {
            hmdb2kegg <- setNames(mapping_df$KEGG, mapping_df$HMDB)
            gmt <- lapply(gmt, function(cpds) {
              mapped <- hmdb2kegg[cpds]
              ifelse(!is.na(mapped) & mapped != "", mapped, cpds)
            })
          }
        }

        # Rename keys to "ID - Description" where descriptions exist
        labeled_names <- make_pathway_labels(names(gmt), desc_map)
        names(gmt) <- labeled_names

        gene_sets <- c(gene_sets, gmt)
        log_message("ssGSEA: Loaded ", length(gmt), " sets from ",
                    basename(gf))
      }
    }
  }

  # From built-in MetaboAnalystR libraries (try to find cached GMT)
  # If no GMT files provided and no sets loaded, warn
  if (length(gene_sets) == 0) {
    log_message("ssGSEA: No GMT files found in config$enrichment$gmt_file. ",
                "Provide GMT pathway files for ssGSEA analysis.")
  }

  gene_sets
}

# ---------- ssGSEA barplot (significant or top pathways) ----------------------

plot_ssgsea_barplot <- function(results_df, config, top_n = 20) {
  sig_df <- results_df[results_df$significant == TRUE, ]

  if (nrow(sig_df) > 0) {
    plot_df <- head(sig_df[order(sig_df$FDR), ], top_n)
    subtitle <- paste0(nrow(sig_df), " significant pathways (FDR < 0.05)")
  } else {
    plot_df <- head(results_df, top_n)
    subtitle <- paste0("Top ", min(top_n, nrow(plot_df)),
                       " pathways (none significant at FDR < 0.05)")
  }

  plot_df$neg_log10_fdr <- -log10(pmax(plot_df$FDR, 1e-20))
  plot_df$pathway_short <- ifelse(
    nchar(plot_df$pathway) > 50,
    paste0(substr(plot_df$pathway, 1, 47), "..."),
    plot_df$pathway
  )
  plot_df$pathway_short <- factor(
    plot_df$pathway_short,
    levels = rev(plot_df$pathway_short)
  )

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = pathway_short, y = neg_log10_fdr, fill = neg_log10_fdr)
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_gradient(
      low = "steelblue", high = "firebrick",
      name = "-log10(FDR)"
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(0.05), linetype = "dashed", color = "grey40"
    ) +
    ggplot2::labs(
      title = "ssGSEA Pathway Enrichment",
      subtitle = subtitle,
      x = NULL,
      y = "-log10(FDR)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 9),
      plot.title = ggplot2::element_text(face = "bold")
    )

  save_plot(p, "ssgsea_barplot.png", config,
            width = 12, height = 8, subdir = "plots")
  p
}

# ---------- ssGSEA boxplots per significant pathway ---------------------------

plot_ssgsea_boxplots <- function(scores, conditions, results_df, config,
                                 max_pathways = 12) {
  sig_pathways <- results_df$pathway[results_df$significant == TRUE]
  if (length(sig_pathways) == 0) {
    log_message("ssGSEA: No significant pathways for boxplots")
    return(NULL)
  }

  sig_pathways <- head(sig_pathways, max_pathways)
  sig_scores <- scores[sig_pathways, , drop = FALSE]

  plot_data <- do.call(rbind, lapply(sig_pathways, function(pw) {
    data.frame(
      pathway = pw,
      score = sig_scores[pw, ],
      condition = as.character(conditions),
      stringsAsFactors = FALSE
    )
  }))

  # Shorten long pathway names for facet labels
  plot_data$pathway_short <- ifelse(
    nchar(plot_data$pathway) > 40,
    paste0(substr(plot_data$pathway, 1, 37), "..."),
    plot_data$pathway
  )

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = condition, y = score, fill = condition)
  ) +
    ggplot2::geom_boxplot(outlier.shape = 16, outlier.size = 1.5) +
    ggplot2::geom_jitter(width = 0.15, size = 1, alpha = 0.5) +
    ggplot2::facet_wrap(~pathway_short, scales = "free_y") +
    ggplot2::labs(
      title = "ssGSEA Pathway Scores by Condition",
      subtitle = paste0("Top ", length(sig_pathways),
                        " significant pathways (FDR < 0.05)"),
      x = NULL,
      y = "ssGSEA Enrichment Score",
      fill = "Condition"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(face = "bold")
    )

  n_facets <- length(sig_pathways)
  plot_height <- max(6, ceiling(n_facets / 3) * 3.5)

  save_plot(p, "ssgsea_boxplots.png", config,
            width = 12, height = plot_height, subdir = "plots")
  p
}
