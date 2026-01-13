# Shared Utilities for Multi-Omics Pipelines

This directory contains shared R utilities used across all pipelines. These modules provide organism-agnostic functionality that prioritizes support for non-model organisms.

## Modules

### `R/annotation_utils.R` - Unified Annotation System

Provides a flexible annotation system with a fallback chain: custom file → biomaRt → OrgDb → none.

**Key Functions:**
- `annotate_features()` - Main annotation function with graceful fallback
- `get_annotation_coverage()` - Check what percentage of features were annotated
- `check_enrichment_ready()` - Verify sufficient annotation for pathway analysis

**Usage:**
```r
source("shared/R/annotation_utils.R")

# Annotate with fallback chain
annotations <- annotate_features(
    feature_ids = gene_ids,
    config = list(
        organism = "Sus scrofa",  # pig
        annotation = list(
            custom_mapping_file = "my_mapping.csv",  # checked first
            fallback_chain = c("custom", "biomart", "orgdb")
        )
    )
)

# Check coverage
coverage <- get_annotation_coverage(annotations)
if (coverage < 0.3) {
    message("Low annotation coverage - consider providing custom mapping")
}
```

### `R/organism_detection.R` - Auto-Detection Utilities

Automatically detects organism and ID type from gene/protein IDs.

**Key Functions:**
- `detect_id_type()` - Detect Ensembl, Entrez, UniProt, symbols, etc.
- `detect_organism_from_ids()` - Guess organism from ID patterns
- `validate_ids()` - Check for duplicates, NAs, mixed types
- `normalize_organism_name()` - Standardize organism names
- `is_model_organism()` - Check if good annotation available

**Usage:**
```r
source("shared/R/organism_detection.R")

# Detect ID type
id_type <- detect_id_type(gene_ids)
# Returns: "ensembl_gene_id", "symbol", "entrez", etc.

# Detect organism
org <- detect_organism_from_ids(gene_ids)
# Returns: list(organism = "Mus musculus", confidence = "high", evidence = "...")

# Validate IDs
validation <- validate_ids(gene_ids)
if (length(validation$warnings) > 0) {
    message("ID issues found: ", paste(validation$warnings, collapse = "; "))
}
```

### `R/gmt_utils.R` - GMT File Handling

Read, write, validate, and generate GMT files for pathway analysis.

**Key Functions:**
- `read_gmt()` / `write_gmt()` - I/O for GMT files
- `validate_gmt()` - Check GMT coverage against your data
- `generate_gmt_from_go()` - Create GO GMT for any organism via biomaRt
- `generate_gmt_from_kegg()` - Create KEGG GMT from KEGG API
- `merge_gmt()` - Combine multiple GMT files
- `convert_gmt_ids()` - Convert gene IDs in GMT using a mapping table

**Usage:**
```r
source("shared/R/gmt_utils.R")

# Read and validate a GMT file
pathways <- read_gmt("pathways.gmt")
validation <- validate_gmt(pathways, my_gene_ids, min_coverage = 0.1)

if (!validation$valid) {
    message("GMT coverage too low - ", validation$coverage * 100, "%")
}

# Use filtered pathways for enrichment
filtered_pathways <- validation$filtered_pathways
```

### `R/config_validation.R` - Configuration Validation

Validate pipeline configurations before running.

**Key Functions:**
- `validate_config()` - Full validation with detailed report
- `is_valid_config()` - Quick TRUE/FALSE check

**Usage:**
```r
source("shared/R/config_validation.R")

config <- yaml::read_yaml("config.yml")
result <- validate_config(config, pipeline_type = "rnaseq", base_path = ".")

if (!result$valid) {
    stop("Configuration errors:\n", paste(result$errors, collapse = "\n"))
}
```

### `R/preflight_checks.R` - Data Validation

Validate data files before running pipelines.

**Key Functions:**
- `run_preflight_checks()` - Comprehensive data validation

**Checks performed:**
- Matrix format (numeric columns, no all-NA rows)
- Sample ID matching between data and metadata
- Group sizes for statistical validity
- ID type detection and validation
- Missing value patterns

**Usage:**
```r
source("shared/R/preflight_checks.R")

config <- yaml::read_yaml("config.yml")
result <- run_preflight_checks(config, pipeline_type = "rnaseq")

if (!result$passed) {
    stop("Data issues found:\n", paste(result$errors, collapse = "\n"))
}

# View statistics
print(result$stats)
# $n_genes, $n_samples, $group_sizes, etc.
```

## Scripts

### `scripts/generate_gmt.R` - GMT Generation Tool

Command-line tool and R functions to generate GMT files for any organism.

**Command Line:**
```bash
# Generate GO GMT for zebrafish
Rscript scripts/generate_gmt.R --organism "Danio rerio" --ontology BP

# Generate KEGG GMT for pig
Rscript scripts/generate_gmt.R --kegg-org ssc

# List available organisms
Rscript scripts/generate_gmt.R --list-organisms --search "cattle"
Rscript scripts/generate_gmt.R --list-kegg --search "pig"
```

**From R:**
```r
source("scripts/generate_gmt.R")

# Generate GO GMT
generate_gmt_for_organism(
    organism = "Bos taurus",
    ontology = "BP",
    id_type = "symbol",
    output_file = "cow_go_bp.gmt"
)

# Generate KEGG GMT
generate_kegg_gmt(
    kegg_org = "bta",
    output_file = "cow_kegg.gmt"
)

# List available organisms
list_available_organisms(pattern = "sheep")
list_kegg_organisms(pattern = "ovis")
```

## Non-Model Organism Support

All utilities are designed with non-model organisms in mind:

1. **Custom files take priority** - User-provided mapping/GMT files are checked first
2. **Graceful degradation** - If annotation fails, pipelines continue with original IDs
3. **Flexible ID types** - Support for any gene/protein ID system
4. **biomaRt fallback** - Works with any organism in Ensembl (100+ species)
5. **GMT generation** - Create pathway files for organisms lacking pre-built databases

### Recommended Workflow for Non-Model Organisms

1. **Check ID type**:
   ```r
   source("shared/R/organism_detection.R")
   detect_id_type(my_gene_ids)
   ```

2. **Generate GMT file**:
   ```bash
   Rscript scripts/generate_gmt.R --organism "My organism" --output my_pathways.gmt
   ```

3. **Validate GMT coverage**:
   ```r
   source("shared/R/gmt_utils.R")
   validate_gmt("my_pathways.gmt", my_gene_ids)
   ```

4. **Configure pipeline**:
   ```yaml
   organism: "My organism"
   annotation:
     custom_gmt_file: "my_pathways.gmt"
     skip_annotation: false  # or true if no annotation needed
   ```

5. **Run preflight checks**:
   ```r
   source("shared/R/preflight_checks.R")
   run_preflight_checks(config, "rnaseq")
   ```

## Supported Organisms

### Ensembl (100+ species via biomaRt)
Common examples with OrgDb packages:
- Homo sapiens, Mus musculus, Rattus norvegicus
- Danio rerio, Drosophila melanogaster, C. elegans
- Arabidopsis thaliana, Saccharomyces cerevisiae

Livestock and other species (biomaRt only):
- Sus scrofa (pig), Bos taurus (cow), Ovis aries (sheep)
- Gallus gallus (chicken), Equus caballus (horse)
- Canis familiaris (dog), Felis catus (cat)
- And many more...

### KEGG (8000+ organisms)
Use organism codes like:
- `hsa` (human), `mmu` (mouse), `rno` (rat)
- `ssc` (pig), `bta` (cow), `oar` (sheep)
- See full list: https://www.genome.jp/kegg/catalog/org_list.html

## Adding Utilities to Pipelines

To use these utilities in a pipeline's R code:

```r
# At the top of the R file
shared_dir <- file.path(dirname(getwd()), "shared", "R")
source(file.path(shared_dir, "annotation_utils.R"))
source(file.path(shared_dir, "organism_detection.R"))
source(file.path(shared_dir, "gmt_utils.R"))
```

Or in `_targets.R`:
```r
# Source shared utilities
tar_source("../shared/R/")
```
