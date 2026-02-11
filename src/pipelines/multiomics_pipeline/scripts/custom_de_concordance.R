library(dplyr)
library(ggplot2)
library(readr)

# Define paths
protein_file <- "data/de_protein_PFOS_c_elegants.csv"
rna_file <- "data/de_rna_PFOS_c_elegants.csv"
output_plot <- "outputs/custom_rna_protein_concordance.png"

# Ensure output directory exists
dir.create(dirname(output_plot), recursive = TRUE, showWarnings = FALSE)

# Load data
message("Loading data...")
protein_data <- read_csv(protein_file, show_col_types = FALSE)
rna_data <- read_csv(rna_file, show_col_types = FALSE)

# Check columns
if (!"log2FC" %in% colnames(protein_data)) {
    stop("Protein data missing 'log2FC' column. Available columns: ", paste(colnames(protein_data), collapse = ", "))
}
if (!"log2FC" %in% colnames(rna_data)) {
    stop("RNA data missing 'log2FC' column. Available columns: ", paste(colnames(rna_data), collapse = ", "))
}

# Merge data (ensure we only keep common genes)
# Protein has 'Wormbase_id', RNA has 'gene_id'
message("Merging data...")
merged_data <- inner_join(protein_data, rna_data, by = c("Wormbase_id" = "gene_id"), suffix = c("_protein", "_rna"))

message(paste("Found", nrow(merged_data), "common genes."))

if (nrow(merged_data) < 3) {
    stop("Not enough common genes to calculate correlation.")
}

# Calculate correlation
correlation <- cor(merged_data$log2FC_protein, merged_data$log2FC_rna, use = "complete.obs", method = "pearson")
message(paste("Correlation:", round(correlation, 4)))

# Create plot
p <- ggplot(merged_data, aes(x = log2FC_rna, y = log2FC_protein)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", color = "red") +
    theme_minimal() +
    labs(
        title = paste("RNA vs Protein log2FC (Pearson r =", round(correlation, 3), ")"),
        x = "RNA log2FC",
        y = "Protein log2FC"
    )

# Save plot
ggsave(output_plot, plot = p, width = 8, height = 6)
message(paste("Plot saved to:", output_plot))
