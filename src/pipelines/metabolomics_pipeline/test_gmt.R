source("R/07_enrichment.R")

# Create a dummy GMT file
cat("Pathway1\tDesc\tC00001\tC00002\nPathway2\tDesc\tC00003\n", file = "test.gmt")

# Test read_gmt_list
gmt <- read_gmt_list("test.gmt")
print(gmt)

# Clean up
unlink("test.gmt")
