# Run unit tests with: Rscript -e 'testthat::test_dir("tests/testthat")'
if (requireNamespace("testthat", quietly = TRUE)) {
  testthat::test_dir("tests/testthat")
} else {
  message("Install testthat to run unit tests: install.packages('testthat')")
}