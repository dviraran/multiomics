#!/bin/bash
# Example script to run MOFA+ analysis

# This is an example of how to run the MOFA+ script from the command line
# Adjust paths and parameters as needed for your analysis

# Set paths (adjust these to your actual data locations)
RNA_FILE="data/abundance_RNA.csv"
PROT_FILE="data/abundance_protein.csv"
METABO_FILE="data/abundance_metabo.csv"

OUTPUT_MODEL="reports/mofa_model.hdf5"
OUTPUT_DIR="reports/mofa_outputs"

# Create output directory if it doesn't exist
mkdir -p reports/mofa_outputs

# Run MOFA with all three omics layers
python3 scripts/run_mofa.py \
  --views \
    rna="${RNA_FILE}" \
    protein="${PROT_FILE}" \
    metabolomics="${METABO_FILE}" \
  --outfile "${OUTPUT_MODEL}" \
  --outdir "${OUTPUT_DIR}" \
  --factors 15 \
  --seed 42 \
  --max_iter 1000 \
  --convergence_mode fast \
  --verbose

echo ""
echo "MOFA analysis complete!"
echo "Results saved to: ${OUTPUT_DIR}"
echo ""
echo "Output files:"
echo "  - ${OUTPUT_MODEL} (trained model)"
echo "  - ${OUTPUT_DIR}/factors.csv (latent factors)"
echo "  - ${OUTPUT_DIR}/weights_rna.csv (RNA loadings)"
echo "  - ${OUTPUT_DIR}/weights_protein.csv (Protein loadings)"
echo "  - ${OUTPUT_DIR}/weights_metabolomics.csv (Metabolomics loadings)"
echo "  - ${OUTPUT_DIR}/variance_explained_*.csv (R2 statistics)"
