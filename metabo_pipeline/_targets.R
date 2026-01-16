# _targets.R: Main pipeline definition for MetaboAnalystR LC-MS workflow
source("packages.R")
lapply(list.files("R", full.names = TRUE, pattern = "\\.R$"), source)
library(targets)
library(tarchetypes)

# Set target options
tar_option_set(
  packages = c("MetaboAnalystR", "ggplot2", "plotly", "dplyr", "readr", "tibble", "stringr", "purrr", "tidyr"),
  format = "rds"
)

list(
  # 1. Raw Data Input & Spectra Processing
  tar_target(
    raw_files,
    list.files("data", pattern = "\\.(mzML|mzXML|cdf)$", full.names = TRUE)
  ),
  tar_target(
    mset_raw,
    init_mset(raw_files)
  ),
  tar_target(
    feature_table,
    process_spectra(mset_raw)
  ),

  # 2. Data Integrity Check & Cleaning
  tar_target(
    cleaned_data,
    clean_features(feature_table)
  ),

  # 3. Normalization & Transformation
  tar_target(
    norm_data,
    normalize_data(cleaned_data)
  ),

  # 4. Statistical & Biomarker Analysis
  tar_target(
    stats_results,
    run_stats(norm_data)
  ),

  # 5. Compound Identification & Annotation
  tar_target(
    annotated_table,
    annotate_compounds(norm_data)
  ),

  # 6. Functional & Pathway Interpretation
  tar_target(
    pathway_results,
    pathway_analysis(annotated_table)
  ),

  # Quarto Reports
  tar_quarto(
    report_raw,
    path = "reports/raw_data_report.qmd",
    execute_params = list(feature_table = feature_table)
  ),
  tar_quarto(
    report_clean,
    path = "reports/cleaning_report.qmd",
    execute_params = list(cleaned_data = cleaned_data)
  ),
  tar_quarto(
    report_norm,
    path = "reports/normalization_report.qmd",
    execute_params = list(norm_data = norm_data)
  ),
  tar_quarto(
    report_stats,
    path = "reports/statistics_report.qmd",
    execute_params = list(stats_results = stats_results)
  ),
  tar_quarto(
    report_annot,
    path = "reports/annotation_report.qmd",
    execute_params = list(annotated_table = annotated_table)
  ),
  tar_quarto(
    report_pathway,
    path = "reports/pathway_report.qmd",
    execute_params = list(pathway_results = pathway_results)
  )
)
