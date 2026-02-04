# Use Posit Package Manager for pre-compiled binaries
options(repos = c(
  CRAN = "https://packagemanager.posit.co/cran/__linux__/jammy/latest"
))

# Install essential packages for report rendering
packages <- c("tidyverse", "rmarkdown", "knitr", "yaml")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing ", pkg, "...")
    install.packages(pkg, quiet = FALSE)
  } else {
    message(pkg, " already installed")
  }
}

message("\nPackages installed successfully!")
