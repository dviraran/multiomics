# Multi-Omics Pipeline Visualization Enhancements

## Summary of Changes

This document summarizes the enhancements made to the Multi-Omics Pipeline based on inspiration from `BulkRNAseq_ATGC/render.Rmd`.

---

## 1. Interactive Top N Terms Selection

### Backend Changes ([R/10_enrichment.R](R/10_enrichment.R))

- **Modified `run_gene_enrichment()`**: Now supports multiple gene set collections and generates plots for different top N values
- **Added `run_enrichment_for_collection()`**: Runs enrichment for a specific collection (GO, KEGG, MSigDB, or custom GMT)
- **Added `run_gene_enrichment_legacy()`**: Fallback for backward compatibility when no collections are defined
- **Added `read_gmt_file()`**: Helper function to read custom GMT files

**Key Features:**
- Automatically generates multiple plots for each collection: `_dotplot_top10.png`, `_dotplot_top20.png`, `_dotplot_top50.png`
- Saves all plots to `output/plots/` directory
- Configurable via `config.yml` `top_n_terms` parameter

### Frontend Changes ([reports/analysis_report.qmd](reports/analysis_report.qmd))

- **Interactive Image Switcher**: JavaScript-based dropdown menu to toggle between Top 10, 20, and 50 terms
- **Base64 Embedding**: All images embedded directly in HTML for self-contained reports
- **Dynamic Plot IDs**: Unique IDs for each omics-collection pair prevent conflicts

**Usage in Report:**
- Each enrichment result displays a dropdown menu above the plot
- Select "Top 10 terms", "Top 20 terms", or "Top 50 terms" to switch views
- Plot updates instantly without page reload

---

## 2. Downloadable Tables with DT::datatable

### Changes

**All static `kable()` tables replaced with interactive `DT::datatable()`:**

- ✅ MAE Summary
- ✅ Sample Alignment
- ✅ Sample Metadata
- ✅ Feature Selection Summary
- ✅ MOFA Variance Explained
- ✅ MOFA Factor Associations
- ✅ DIABLO Selected Features
- ✅ DIABLO CV Error Rates
- ✅ SNF Cluster Assignments
- ✅ RNA-Protein Concordance Summary
- ✅ Enrichment Results (all collections)
- ✅ Combined Multi-Omics Enrichment
- ✅ Methods Configuration
- ✅ Feature Selection Settings

**Features:**
- **Download Buttons**: Copy, CSV, Excel, PDF, Print
- **Scrollable**: Horizontal scroll for wide tables
- **Searchable**: Built-in search functionality
- **Sortable**: Click column headers to sort
- **Paginated**: Show 10-20 rows per page by default

**Example:**
```r
DT::datatable(data,
  caption = "Table Caption",
  rownames = FALSE,
  extensions = 'Buttons',
  options = list(
    dom = 'Bfrtip',
    buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
    scrollX = TRUE,
    pageLength = 15
  ))
```

---

## 3. Multi-Collection GSEA Support

### Config Changes ([config.yml](config.yml))

**New `enrichment.collections` section:**

```yaml
enrichment:
  run_enrichment: yes
  methods: gsea
  ora_pvalue: 0.1
  min_set_size: 10
  max_set_size: 500
  combine_method: fisher
  use_kegg: yes
  # Multiple gene set collections for enrichment
  collections:
    - name: "GO_BP"
      type: "GO"
      ont: "BP"
    - name: "KEGG"
      type: "KEGG"
    - name: "Hallmark"
      type: "msigdbr"
      category: "H"
    - name: "Curated_Pathways"
      type: "msigdbr"
      category: "C2"
      subcategory: "CP:KEGG"
    - name: "Custom"
      type: "gmt"
      path: "path/to/custom.gmt"
  # Top N terms to display in plots
  top_n_terms: [10, 20, 50]
```

**Supported Collection Types:**

1. **GO (Gene Ontology)**
   ```yaml
   - name: "GO_BP"
     type: "GO"
     ont: "BP"  # BP, CC, or MF
   ```

2. **KEGG Pathways**
   ```yaml
   - name: "KEGG"
     type: "KEGG"
   ```

3. **MSigDB Collections (via msigdbr)**
   ```yaml
   - name: "Hallmark"
     type: "msigdbr"
     category: "H"

   - name: "Curated_Pathways"
     type: "msigdbr"
     category: "C2"
     subcategory: "CP:KEGG"  # Optional
   ```

   **Available MSigDB Categories:**
   - `H`: Hallmark gene sets
   - `C1`: Positional gene sets
   - `C2`: Curated gene sets (includes CP:KEGG, CP:REACTOME, CP:BIOCARTA, etc.)
   - `C3`: Regulatory target gene sets
   - `C4`: Computational gene sets
   - `C5`: Ontology gene sets
   - `C6`: Oncogenic signature gene sets
   - `C7`: Immunologic signature gene sets
   - `C8`: Cell type signature gene sets

4. **Custom GMT Files**
   ```yaml
   - name: "Custom"
     type: "gmt"
     path: "path/to/custom.gmt"
   ```

### Backend Implementation ([R/10_enrichment.R](R/10_enrichment.R))

**For each collection:**
1. Maps gene IDs to ENTREZID (required for enrichment)
2. Runs enrichment using appropriate method:
   - `clusterProfiler::enrichGO()` for GO
   - `clusterProfiler::enrichKEGG()` for KEGG
   - `clusterProfiler::enricher()` with `msigdbr` term2gene for MSigDB
   - `run_simple_ora()` for custom GMT files
3. Saves results: `{omic}_{collection}_enrichment.csv`
4. Generates plots: `{omic}_{collection}_dotplot_top{N}.png`

### Frontend Implementation ([reports/analysis_report.qmd](reports/analysis_report.qmd))

**Nested Tabs Structure:**
```
Per-Omics Enrichment
├── Transcriptomics
│   ├── GO_BP
│   │   ├── [Interactive Plot Switcher]
│   │   └── [Downloadable Table]
│   ├── KEGG
│   │   ├── [Interactive Plot Switcher]
│   │   └── [Downloadable Table]
│   ├── Hallmark
│   │   ├── [Interactive Plot Switcher]
│   │   └── [Downloadable Table]
│   └── Curated_Pathways
│       ├── [Interactive Plot Switcher]
│       └── [Downloadable Table]
├── Proteomics
│   ├── GO_BP
│   └── ...
└── Metabolomics
    └── ...
```

**Fallback Behavior:**
- If `enrichment.collections` is not defined, falls back to legacy behavior (GO + KEGG)
- Ensures backward compatibility with existing configs

---

## Installation Requirements

### R Packages

**Core Packages (already in pipeline):**
- `DT`: Interactive tables
- `htmltools`: HTML rendering
- `jsonlite`: JSON handling for JavaScript integration

**New Package for MSigDB Support:**
```r
install.packages("BiocManager")
BiocManager::install("msigdbr")
```

**Organism-Specific Databases:**
```r
# For human
BiocManager::install("org.Hs.eg.db")

# For C. elegans
BiocManager::install("org.Ce.eg.db")
```

---

## Usage Examples

### Example 1: Basic Configuration with Multiple Collections

```yaml
enrichment:
  run_enrichment: yes
  collections:
    - name: "GO_BP"
      type: "GO"
      ont: "BP"
    - name: "Hallmark"
      type: "msigdbr"
      category: "H"
  top_n_terms: [10, 20, 50]
```

**Output:**
- `transcriptomics_GO_BP_enrichment.csv`
- `transcriptomics_GO_BP_dotplot_top10.png`
- `transcriptomics_GO_BP_dotplot_top20.png`
- `transcriptomics_GO_BP_dotplot_top50.png`
- `transcriptomics_Hallmark_enrichment.csv`
- `transcriptomics_Hallmark_dotplot_top10.png`
- `transcriptomics_Hallmark_dotplot_top20.png`
- `transcriptomics_Hallmark_dotplot_top50.png`

### Example 2: Using Custom GMT File

```yaml
enrichment:
  collections:
    - name: "MyPathways"
      type: "gmt"
      path: "data/my_pathways.gmt"
  top_n_terms: [10, 20]
```

**GMT File Format:**
```
PATHWAY_1	Description	GENE1	GENE2	GENE3	...
PATHWAY_2	Description	GENE4	GENE5	GENE6	...
```

### Example 3: C2 KEGG Pathways Only

```yaml
enrichment:
  collections:
    - name: "KEGG_C2"
      type: "msigdbr"
      category: "C2"
      subcategory: "CP:KEGG"
  top_n_terms: [20]
```

---

## Key Code Locations

### Modified Files:
1. **[config.yml](config.yml)** (lines 117-144)
   - Added `collections` section
   - Added `top_n_terms` parameter

2. **[R/10_enrichment.R](R/10_enrichment.R)**
   - Lines 118-245: `run_gene_enrichment()` and `run_enrichment_for_collection()`
   - Lines 929-947: `read_gmt_file()` helper function

3. **[reports/analysis_report.qmd](reports/analysis_report.qmd)**
   - Lines 136-150: DT setup in setup chunk
   - Lines 275-297: Interactive tables (sample data)
   - Lines 512-680: **Main enrichment section with interactive switcher and tabs**

---

## Inspiration Sources

All features adapted from `…\ozsol\BulkRNAseq_ATGC\render.Rmd`:

1. **Interactive PCA Switcher** (lines 220-302 in render.Rmd)
   - Adapted for enrichment plots
   - Uses base64 image embedding
   - JavaScript dropdown logic

2. **DT Tables with Download Buttons** (lines 145-164, 448-462 in render.Rmd)
   - Consistent `dom = 'Bfrtip'` configuration
   - Button set: copy, csv, excel, pdf, print
   - ScrollX for wide tables

3. **GSEA Multi-Collection Tabs** (lines 698-751 in render.Rmd)
   - Nested tabset structure
   - Dynamic tab generation
   - Collection-specific styling

---

## Testing Checklist

- [ ] Run pipeline with new config
- [ ] Verify plots generated for all top_n values
- [ ] Check interactive switcher in HTML report
- [ ] Test download buttons on all tables
- [ ] Verify all collection types work (GO, KEGG, MSigDB, GMT)
- [ ] Test with C. elegans and human organisms
- [ ] Check fallback behavior (no collections defined)
- [ ] Verify combined enrichment table displays correctly

---

## Troubleshooting

### Issue: "Package 'msigdbr' not available"
**Solution:** Install msigdbr
```r
BiocManager::install("msigdbr")
```

### Issue: "No genes mapped for [collection]"
**Solution:**
- Check organism setting in config matches your data
- Verify gene IDs are SYMBOL format (or WORMBASE for C. elegans)
- Ensure OrgDb package is installed for your organism

### Issue: "Interactive plot switcher not working"
**Solution:**
- Ensure `jsonlite` package is installed
- Check browser console for JavaScript errors
- Verify plot files exist in `output/plots/`

### Issue: "DT tables not displaying download buttons"
**Solution:**
- Ensure `DT` package version >= 0.18
- Check `extensions = 'Buttons'` is set
- Verify `dom = 'Bfrtip'` is in options

---

## Future Enhancements

Potential additions based on user feedback:
- [ ] Interactive heatmaps for enrichment results
- [ ] Network visualization for pathway relationships
- [ ] Gene set variation analysis (GSVA) scores
- [ ] Leading edge gene tables
- [ ] Enrichment map integration

---

**Date:** 2026-02-05
**Pipeline Version:** Multi-Omics Integration Pipeline v2.0
**Author:** Enhanced by Claude Sonnet 4.5
