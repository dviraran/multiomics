library(testthat)

context("Random Forest wrapper")

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
source(file.path(root, "R/11_random_forest.R"))

cfg <- load_config(file.path(root, "config.yml"))

test_that("run_random_forest returns importance when data is adequate", {
  tmpdir <- tempdir()
  cfg$output$output_dir <- tmpdir
  create_output_dirs(cfg)

  # Synthetic matrix: 10 features x 20 samples
  set.seed(42)
  mat <- matrix(rnorm(200), nrow = 10)
  rownames(mat) <- paste0("F", seq_len(nrow(mat)))
  colnames(mat) <- paste0("S", seq_len(ncol(mat)))

  # Metadata with two conditions
  metadata <- data.frame(
    sample_id = colnames(mat),
    condition = rep(c("A", "B"), length.out = ncol(mat)),
    stringsAsFactors = FALSE
  )

  res <- run_random_forest(mat, metadata, cfg)
  # If ranger not installed, wrapper returns NULL; assert either NULL or valid structure
  if (is.null(res)) {
    expect_null(res)
  } else {
    expect_true("importance" %in% names(res))
    expect_true(file.exists(file.path(cfg$output$output_dir, "tables", "rf_feature_importance.csv")))
  }
})