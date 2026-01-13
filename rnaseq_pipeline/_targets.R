# _targets.R
# Main targets pipeline for RNA-seq differential expression and pathway analysis
# Enhanced with batch correction, WGCNA, GSVA, and cell type deconvolution

library(targets)
library(tarchetypes)

# Source all functions
tar_source("R/")

# Note: New analysis scripts added:
# - R/04b_batch_correction.R - Batch effect detection and correction
# - R/07b_advanced_pathway.R - GSVA and extended pathway databases
# - R/07c_coexpression_networks.R - WGCNA co-expression analysis
# - R/05b_deconvolution.R - xCell 2.0 cell type deconvolution

# Load configuration
config <- yaml::read_yaml("config.yml")

# Set target options
tar_option_set(
  packages = c(
    "tidyverse", "DESeq2", "edgeR", "limma", "AnnotationHub", "biomaRt",
    "clusterProfiler", "fgsea", "enrichplot", "ggplot2", "pheatmap",
    "ggrepel", "RColorBrewer", "here", "yaml", "janitor", "matrixStats",
    "jsonlite", "htmltools"
  ),
  format = "qs",
  error = "continue"
)

# Define the pipeline
list(
  # ==========================================================================
  # Stage 0: Load configuration
  # ==========================================================================
  tar_target(
    pipeline_config,
    load_params("config.yml")
  ),

  # ==========================================================================
  # Stage 1: Data ingestion and validation
  # ==========================================================================
  tar_target(
    counts_raw,
    read_counts_matrix(pipeline_config$counts_file)
  ),

  tar_target(
    metadata_raw,
    read_metadata(pipeline_config$metadata_file, sample_id_col = pipeline_config$sample_id_col)
  ),

  tar_target(
    validated_data,
    validate_inputs(
      counts = counts_raw,
      metadata = metadata_raw,
      sample_id_col = pipeline_config$sample_id_col
    )
  ),

  tar_target(
    counts_validated,
    validated_data$counts
  ),

  tar_target(
    metadata_validated,
    validated_data$metadata
  ),

  # ==========================================================================
  # Stage 2: Gene ID processing and annotation
  # ==========================================================================
  tar_target(
    gene_ids_processed,
    process_gene_ids(
      counts = counts_validated,
      strip_version = pipeline_config$strip_ensembl_version
    )
  ),

  tar_target(
    counts_dedup,
    gene_ids_processed$counts
  ),

  tar_target(
    gene_id_mapping_raw,
    gene_ids_processed$id_mapping
  ),

  tar_target(
    annotation_result,
    annotate_genes(
      gene_ids = rownames(counts_dedup),
      organism = pipeline_config$organism,
      gene_id_type = pipeline_config$gene_id_type,
      mapping_file = pipeline_config$mapping_file,
      annotation_source = pipeline_config$annotation_source
    )
  ),

  tar_target(
    gene_annotation,
    annotation_result$annotation
  ),

  tar_target(
    annotation_stats,
    annotation_result$stats
  ),

  tar_target(
    annotation_csv,
    save_csv(gene_annotation, here::here("outputs", "gene_annotation.csv")),
    format = "file"
  ),

  # ==========================================================================
  # Stage 3: Filtering low-expression genes
  # ==========================================================================
  tar_target(
    filtering_result,
    filter_low_expression(
      counts = counts_dedup,
      metadata = metadata_validated,
      group_col = pipeline_config$group_col,
      min_count = pipeline_config$min_count,
      min_samples = pipeline_config$min_samples
    )
  ),

  tar_target(
    counts_filtered,
    filtering_result$counts
  ),

  tar_target(
    filtering_stats,
    filtering_result$stats
  ),

  # ==========================================================================
  # Stage 4: Normalization and transformation
  # ==========================================================================
  tar_target(
    dds_object,
    create_deseq_object(
      counts = counts_filtered,
      metadata = metadata_validated,
      design_formula = as.formula(pipeline_config$design_formula)
    )
  ),

  tar_target(
    dds_normalized,
    normalize_deseq(dds_object)
  ),

  tar_target(
    normalized_counts,
    get_normalized_counts(dds_normalized)
  ),

  tar_target(
    vst_counts,
    get_vst_counts(dds_normalized)
  ),

  tar_target(
    normalized_counts_csv,
    save_csv(
      as.data.frame(normalized_counts) %>%
        rownames_to_column("gene_id"),
      here::here("outputs", "normalized_counts.csv")
    ),
    format = "file"
  ),

  tar_target(
    vst_counts_csv,
    save_csv(
      as.data.frame(vst_counts) %>%
        rownames_to_column("gene_id"),
      here::here("outputs", "vst_counts.csv")
    ),
    format = "file"
  ),

  # ==========================================================================
  # Stage 4b: Batch Effect Detection and Correction
  # ==========================================================================
  tar_target(
    batch_analysis,
    {
      bc <- pipeline_config$batch_correction %||% list()
      if (bc$detect_batch %||% TRUE) {
        run_batch_analysis(
          vst_counts = vst_counts,
          counts_filtered = counts_filtered,
          metadata = metadata_validated,
          config = pipeline_config
        )
      } else {
        NULL
      }
    }
  ),

  tar_target(
    batch_corrected_dds,
    {
      if (!is.null(batch_analysis) && !is.null(batch_analysis$corrected_dds)) {
        batch_analysis$corrected_dds
      } else {
        dds_normalized
      }
    }
  ),

  # ==========================================================================
  # Stage 5: QC and exploratory analysis
  # ==========================================================================
  tar_target(
    qc_metrics,
    compute_qc_metrics(
      counts = counts_filtered,
      normalized_counts = normalized_counts,
      metadata = metadata_validated
    )
  ),

  tar_target(
    qc_metrics_csv,
    save_csv(qc_metrics, here::here("outputs", "qc_metrics.csv")),
    format = "file"
  ),

  tar_target(
    sample_correlation,
    compute_sample_correlation(vst_counts)
  ),

  tar_target(
    sample_correlation_csv,
    save_csv(
      as.data.frame(sample_correlation) %>%
        rownames_to_column("sample"),
      here::here("outputs", "sample_correlation.csv")
    ),
    format = "file"
  ),

  tar_target(
    pca_result,
    compute_pca(
      vst_counts = vst_counts,
      metadata = metadata_validated,
      sample_id_col = pipeline_config$sample_id_col
    )
  ),

  tar_target(
    pca_csv,
    save_csv(pca_result$pca_data, here::here("outputs", "pca_results.csv")),
    format = "file"
  ),

  tar_target(
    outlier_detection,
    detect_outliers(
      vst_counts = vst_counts,
      metadata = metadata_validated,
      sample_id_col = pipeline_config$sample_id_col,
      sd_threshold = pipeline_config$outlier_sd_threshold
    )
  ),

  tar_target(
    outlier_csv,
    save_csv(outlier_detection, here::here("outputs", "outlier_flags.csv")),
    format = "file"
  ),

  tar_target(
    qc_plots,
    generate_qc_plots(
      qc_metrics = qc_metrics,
      sample_correlation = sample_correlation,
      pca_result = pca_result,
      metadata = metadata_validated,
      group_col = pipeline_config$group_col,
      output_dir = here::here("outputs", "plots")
    )
  ),

  # ==========================================================================
  # Stage 6: Differential expression analysis
  # ==========================================================================
  tar_target(
    de_results_list,
    run_differential_expression(
      dds = dds_normalized,
      contrasts = pipeline_config$contrasts,
      group_col = pipeline_config$group_col,
      alpha = pipeline_config$alpha,
      lfc_threshold = pipeline_config$lfc_threshold
    )
  ),

  tar_target(
    de_results_annotated,
    annotate_de_results(
      de_results_list = de_results_list,
      annotation = gene_annotation
    )
  ),

  tar_target(
    de_results_csv,
    save_de_results(
      de_results_annotated,
      output_dir = here::here("outputs", "de_results")
    ),
    format = "file"
  ),

  tar_target(
    de_summary,
    summarize_de_results(
      de_results_annotated,
      alpha = pipeline_config$alpha,
      lfc_threshold = pipeline_config$lfc_threshold
    )
  ),

  tar_target(
    de_plots,
    generate_de_plots(
      de_results_annotated = de_results_annotated,
      vst_counts = vst_counts,
      metadata = metadata_validated,
      group_col = pipeline_config$group_col,
      output_dir = here::here("outputs", "plots")
    )
  ),

  # ==========================================================================
  # Stage 7: Pathway and gene set analysis
  # ==========================================================================
  tar_target(
    gene_sets,
    load_gene_sets(
      organism = pipeline_config$organism,
      pathway_database = pipeline_config$pathway_database,
      gmt_file = pipeline_config$gmt_file,
      annotation = gene_annotation,
      gene_id_type = pipeline_config$gene_id_type
    )
  ),

  tar_target(
    pathway_results,
    run_pathway_analysis(
      de_results = de_results_annotated,
      gene_sets = gene_sets,
      annotation = gene_annotation,
      method = pipeline_config$pathway_method,
      min_size = pipeline_config$pathway_min_size,
      max_size = pipeline_config$pathway_max_size
    )
  ),

  tar_target(
    pathway_csv,
    save_pathway_results(
      pathway_results,
      output_dir = here::here("outputs", "pathway_results")
    ),
    format = "file"
  ),

  tar_target(
    pathway_plots,
    generate_pathway_plots(
      pathway_results = pathway_results,
      output_dir = here::here("outputs", "plots")
    )
  ),

  # ==========================================================================
  # Stage 7b: Advanced Pathway Analysis (GSVA + Extended Databases)
  # ==========================================================================
  tar_target(
    gsva_results,
    {
      ap <- pipeline_config$advanced_pathway %||% list()
      if (ap$run_gsva %||% TRUE) {
        run_advanced_pathway_analysis(
          normalized_counts = vst_counts,
          de_results_annotated = de_results_annotated,
          gene_annotation = annotation_stats,
          metadata = metadata_validated,
          config = pipeline_config
        )
      } else {
        NULL
      }
    }
  ),

  tar_target(
    gsva_scores_csv,
    {
      if (!is.null(gsva_results) && !is.null(gsva_results$gsva_scores)) {
        out_path <- here::here("outputs", "gsva_scores.csv")
        write.csv(gsva_results$gsva_scores, out_path, row.names = TRUE)
        out_path
      } else {
        NULL
      }
    },
    format = "file"
  ),

  # ==========================================================================
  # Stage 7c: Co-Expression Network Analysis (WGCNA)
  # ==========================================================================
  tar_target(
    wgcna_results,
    {
      wc <- pipeline_config$coexpression %||% list()
      if (wc$run_wgcna %||% TRUE) {
        run_wgcna_analysis(
          vst_counts = vst_counts,
          metadata = metadata_validated,
          gene_annotation = annotation_stats,
          config = pipeline_config
        )
      } else {
        NULL
      }
    }
  ),

  tar_target(
    wgcna_modules_csv,
    {
      if (!is.null(wgcna_results) && !is.null(wgcna_results$module_membership)) {
        out_path <- here::here("outputs", "wgcna_modules.csv")
        write.csv(wgcna_results$module_membership, out_path, row.names = FALSE)
        out_path
      } else {
        NULL
      }
    },
    format = "file"
  ),

  # ==========================================================================
  # Stage 5b: Cell Type Deconvolution (xCell 2.0)
  # ==========================================================================
  tar_target(
    deconvolution_results,
    {
      dc <- pipeline_config$deconvolution %||% list()
      if (dc$run_deconvolution %||% TRUE) {
        run_deconvolution_analysis(
          dds = batch_corrected_dds,
          de_results = de_results_annotated,
          config = pipeline_config
        )
      } else {
        NULL
      }
    }
  ),

  tar_target(
    xcell2_scores_csv,
    {
      if (!is.null(deconvolution_results) && !is.null(deconvolution_results$scores)) {
        out_path <- here::here("outputs", "xcell2_scores.csv")
        write.csv(deconvolution_results$scores, out_path, row.names = FALSE)
        out_path
      } else {
        NULL
      }
    },
    format = "file"
  ),

  # ==========================================================================
  # Stage 8: Figure Commentary Generation
  # ==========================================================================

  # Build table of all figures with metadata
  tar_target(
    figures_tbl,
    build_figures_table(
      qc_plots = qc_plots,
      de_plots = de_plots,
      pathway_plots = pathway_plots,
      batch_analysis = batch_analysis,
      gsva_results = gsva_results,
      wgcna_results = wgcna_results,
      deconvolution_results = deconvolution_results,
      config = pipeline_config
    )
  ),

  # Generate commentary for all figures
  tar_target(
    commentary_tbl,
    generate_all_commentary(
      figures_tbl = figures_tbl,
      config = pipeline_config,
      qc_metrics = qc_metrics,
      pca_result = pca_result,
      sample_correlation = sample_correlation,
      de_summary = de_summary,
      outlier_detection = outlier_detection,
      batch_analysis = batch_analysis,
      gsva_results = gsva_results,
      wgcna_results = wgcna_results,
      deconvolution_results = deconvolution_results,
      output_dir = here::here("outputs", "commentary")
    )
  ),

  # Save figures table
  tar_target(
    figures_tbl_csv,
    save_csv(figures_tbl, here::here("outputs", "commentary", "figures_metadata.csv")),
    format = "file"
  ),

  # ==========================================================================
  # Stage 9: Report generation
  # ==========================================================================
  tar_render(
    report,
    path = here::here("reports", "analysis_report.Rmd"),
    output_dir = here::here("outputs"),
    params = list(
      config = pipeline_config,
      qc_metrics = qc_metrics,
      filtering_stats = filtering_stats,
      annotation_stats = annotation_stats,
      pca_result = pca_result,
      sample_correlation = sample_correlation,
      outlier_detection = outlier_detection,
      de_results = de_results_annotated,
      de_summary = de_summary,
      pathway_results = pathway_results,
      commentary_tbl = commentary_tbl
    )
  ),

  # ==========================================================================
  # Stage 10: Pipeline summary
  # ==========================================================================
  tar_target(
    pipeline_summary,
    create_pipeline_summary(
      filtering_stats = filtering_stats,
      annotation_stats = annotation_stats,
      de_summary = de_summary,
      pathway_results = pathway_results,
      params = pipeline_config
    )
  ),

  tar_target(
    pipeline_summary_csv,
    save_csv(pipeline_summary, here::here("outputs", "pipeline_summary.csv")),
    format = "file"
  )
)
