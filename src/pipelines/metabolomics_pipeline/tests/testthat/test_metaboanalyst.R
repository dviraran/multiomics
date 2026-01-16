library(testthat)

context("MetaboAnalyst integration helpers")

# Helper: find project root containing config.yml
find_project_root <- function() {
  candidates <- c(".", "..", "../..", "../../..", "../../../..")
  for (p in candidates) {
    if (file.exists(file.path(p, "config.yml"))) return(normalizePath(p))
  }
  stop("Could not find config.yml in parent directories")
}

root <- find_project_root()
source(file.path(root, "R/00_utils.R"))
source(file.path(root, "R/10_metaboanalyst_integration.R"))

cfg <- load_config(file.path(root, "config.yml"))


test_that("init_mset_from_table handles missing files gracefully", {
  # Provide a non-existent file
  mset <- init_mset_from_table("nonexistent_file.tsv", cfg)
  expect_null(mset)
})

test_that("export_normalized_samples_tsv writes TSV when matrix present", {
  tmpdir <- tempdir()
  cfg$output$output_dir <- tmpdir
  create_output_dirs(cfg)
  # Create a small normalized matrix
  mat <- matrix(runif(20), nrow = 4)
  rownames(mat) <- paste0("F", seq_len(nrow(mat)))
  colnames(mat) <- paste0("S", seq_len(ncol(mat)))
  norm_list <- list(matrix = mat)

  out <- export_normalized_samples_tsv(norm_list, cfg)
  expect_true(file.exists(out))
  df <- readr::read_tsv(out)
  expect_true(nrow(df) == ncol(mat))
})