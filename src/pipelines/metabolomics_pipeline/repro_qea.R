tryCatch(
    {
        library(MetaboAnalystR)
        print("MetaboAnalystR loaded successfully")
    },
    error = function(e) {
        print(paste("Error loading MetaboAnalystR:", e$message))
    }
)

# Mock data creation similar to prepare_qea_data
# Create a dummy data file
df <- data.frame(
    Sample = c("S1", "S2", "S3", "S4"),
    Group = c("A", "A", "B", "B"),
    C00022 = c(1.1, 1.2, 0.9, 0.8), # Pyruvate
    C00024 = c(0.5, 0.6, 1.5, 1.4) # Acetyl-CoA
)

tmp_file <- "test_qea_data.txt"
write.table(df, tmp_file, sep = "\t", row.names = FALSE, quote = FALSE)

.libPaths(c("/home/ozsol/miniconda3/pkgs/r-metaboanalyst-2.0.0-r40hdfd78af_2/lib/R/library", .libPaths()))
tryCatch(
    {
        print("Initializing mSet...")
        mSet <- InitDataObjects("conc", "msetqea", FALSE)

        print("Reading data...")
        mSet <- Read.TextData(mSet, tmp_file, "rowu", "disc")

        print("Sanity check...")
        mSet <- SanityCheckData(mSet)
        mSet <- ReplaceMin(mSet)
        mSet <- PreparePrenormData(mSet)
        mSet <- Normalization(mSet, "NULL", "LogNorm", "MeanCenter", "S10T0", ratio = FALSE, ratioNum = 20)

        print("Cross referencing...")
        mSet <- CrossReferencing(mSet, "name")
        mSet <- CreateMappingResultTable(mSet)

        print("Setting Metabolome Filter...")
        mSet <- SetMetabolomeFilter(mSet, FALSE)

        print("Setting Library: kegg_pathway")
        mSet <- SetCurrentMsetLib(mSet, "kegg_pathway", 2)

        print("Calculating Global Test Score...")
        mSet <- CalculateGlobalTestScore(mSet)

        print("Success!")
        print(head(mSet$analSet$qea.mat))
    },
    error = function(e) {
        print(paste("FAILED:", e$message))
        # traceback()
    }
)

file.remove(tmp_file)
