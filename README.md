# atacreportR

> Standalone R Package for ATAC-seq Differential Accessibility Analysis

A comprehensive, modular R package for ATAC-seq data processing, statistical analysis, and visualization. Designed to work both as part of the UFHCC BCB-SR reporting pipeline and as a standalone toolkit for differential accessibility analysis.

## Features

### Data Processing & Analysis
- **Flexible Data Input**: Load from nf-core/atacseq output, existing DESeq2 objects, or raw peak/BAM files
- **Consensus Peak Generation**: Merge peaks across samples with reproducibility filtering
- **Multiple Normalization Methods**: TMM normalization or spike-in normalization
- **Differential Analysis**: edgeR quasi-likelihood framework with configurable filtering
- **Peak Annotation**: ChIPseeker integration with organism-specific databases

### Visualization
- **Interactive Plots**: PCA, MA plots, and volcano plots with plotly tooltips
- **Annotation Barplots**: Genomic distribution of differential peaks
- **Quality Control**: Comprehensive QC tables with FRiP scores and alignment statistics

### Pathway Analysis
- **GO Term Enrichment**: Biological Process, Molecular Function, Cellular Component
- **KEGG Pathway Analysis**: With gene symbol mapping and interactive visualization
- **Parallel Processing**: Configurable cores for enrichment analysis

### Browser Integration
- **WashU Epigenome Browser**: Automatic track generation with BigBed and BigWig files

## Installation

```r
# Install from GitHub
remotes::install_github("UFHCC-BCBSR/pkg-atacreportR")
```

## Quick Start

```r
library(atacreportR)
library(edgeR)

# 1. Load your DESeq2 object (from nf-core or custom)
dds <- prepare_analysis_data(
  sample_sheet = "samples.csv",
  organism = "hsa",
  output_dir = "results/",
  dds_file = "my_dds.RData",  # Or provide peak_files + bam_files
  peak_annotation = "peaks.anno.txt"
)

# 2. Convert to DGEList with TMM normalization
dge <- dds_to_dgelist(dds, norm_method = "TMM")

# 3. Run differential analysis
results_list <- run_differential_analysis(
  dge = dge,
  contrasts = c("Treatment_vs_Control", "KO_vs_WT"),
  min_count = 10,
  min_prop = 0.7
)

# 4. Visualize results
plot_volcano(results_list, sig_cutoff = 0.05, logfc_cutoff = 1)
plot_ma(results_list, sig_cutoff = 0.05, logfc_cutoff = 1)

# 5. Run pathway enrichment
enrichment_data <- prepare_gene_lists_for_enrichment(
  results_list = results_list,
  organism = "hsa",
  sig_cutoff = 0.05,
  logfc_cutoff = 1
)

GO_results <- generate_enrichment_plot_atac(
  gene_lists = enrichment_data$gene_lists,
  de_results_df = enrichment_data$results_df,
  universe_entrez = enrichment_data$universe,
  org_db = enrichment_data$org_db,
  ont_category = "BP",
  n_cores = 4
)
```

## Main Functions

### Data Preparation

#### Core Functions
- `prepare_analysis_data()` - Load or create DESeq2 object from various sources
- `load_sample_metadata()` - Parse and validate sample sheets (nf-core or simple format)
- `load_dds_file()` - Load existing DESeqDataSet from RData
- `load_peak_annotation()` - Load HOMER-format peak annotations
- `load_qc_data()` - Load FRiP scores and flagstat files

#### Building from Scratch
- `make_consensus_peaks()` - Create consensus peaks with reproducibility filtering
- `count_bam_reads_in_peaks()` - Quantify reads using Rsubread::featureCounts
- `create_dds_from_peaks()` - Build DESeqDataSet from peak/BAM files
- `generate_peak_annotation()` - Annotate peaks with ChIPseeker

### Data Conversion & Normalization
- `dds_to_dgelist()` - Convert DESeq2 to edgeR DGEList with TMM/manual/no normalization

### Differential Accessibility Analysis
- `run_differential_analysis()` - edgeR quasi-likelihood testing across contrasts

### Interactive Visualization
- `plot_pca()` - PCA with group ellipses and hover tooltips
- `plot_volcano()` - Interactive volcano plots (one per contrast)
- `plot_ma()` - Interactive MA plots (one per contrast)
- `plot_da_by_annotation()` - Genomic annotation barplots (differential vs non-differential)

### Quality Control & Reporting
- `summarize_atac_sample_qc()` - Comprehensive QC table with FRiP, alignment, peak counts

### Pathway Enrichment
- `prepare_gene_lists_for_enrichment()` - Extract gene lists and universe from results
- `generate_enrichment_plot_atac()` - GO enrichment with interactive plots
- `generate_kegg_enrichment_plot_atac()` - KEGG pathway enrichment

### Genome Browser Integration
- `create_da_bigbed()` - Generate BigBed tracks for WashU Browser

### Validation
- `validate_data_consistency()` - Check data integrity across DDS, sample_info, annotations

## Usage Examples

### Example 1: Starting from nf-core/atacseq Output

```r
# nf-core provides a DDS object and HOMER annotations
dds <- prepare_analysis_data(
  sample_sheet = "samplesheet.csv",
  organism = "mmu",
  output_dir = "analysis/",
  dds_file = "nfcore/consensus_peaks.mLb.clN.dds.RData",
  peak_annotation = "nfcore/consensus_peaks.annotatePeaks.txt",
  bigwig_files = c("s1" = "s1.bigWig", "s2" = "s2.bigWig"),
  qc_flagstat_dir = "nfcore/flagstat/",
  qc_frip_file = "nfcore/frip.txt"
)

# Convert to DGEList - all metadata is already integrated!
dge <- dds_to_dgelist(dds, norm_method = "TMM")

# Everything you need is in the DGEList:
dge$genes        # Peak annotations
dge$samples      # Sample metadata + QC metrics
dge$counts       # Read counts
```

### Example 2: Building from Peak and BAM Files

```r
# Start from scratch with your own peak calls and alignments
dds <- prepare_analysis_data(
  sample_sheet = "samples.csv",
  organism = "hsa",
  output_dir = "results/",
  peak_files = c("s1" = "s1_peaks.bed", "s2" = "s2_peaks.bed"),
  bam_files = c("s1" = "s1.bam", "s2" = "s2.bam"),
  is_paired_end = TRUE,
  n_threads = 8
)

# Annotations are generated automatically with ChIPseeker
# Consensus peaks are created with reproducibility filtering
```

### Example 3: Complete Analysis Workflow

```r
# Load data
dds <- prepare_analysis_data(
  sample_sheet = "samples.csv",
  organism = "hsa",
  output_dir = "results/",
  dds_file = "my_dds.RData"
)

# Normalize
dge <- dds_to_dgelist(dds, norm_method = "TMM")

# Quality control
qc_table <- summarize_atac_sample_qc(
  dds = dds,
  flagstat_dir = "qc/flagstat/"
)

pca_plot <- plot_pca(dge, title = "Sample PCA", plot_var = "condition")

# Differential analysis (2+ contrasts supported)
results <- run_differential_analysis(
  dge = dge,
  contrasts = c("Treated_vs_Control", "KO_vs_WT"),
  min_count = 10,
  min_prop = 0.7
)

# Visualize
plot_volcano(results, sig_cutoff = 0.05, logfc_cutoff = 1)
plot_da_by_annotation(results, sig_cutoff = 0.05, logfc_cutoff = 1)

# Enrichment (handles multiple contrasts automatically)
enrichment <- prepare_gene_lists_for_enrichment(
  results_list = results,
  organism = "hsa",
  sig_cutoff = 0.05,
  logfc_cutoff = 1
)

# GO and KEGG enrichment
GO_BP <- generate_enrichment_plot_atac(
  gene_lists = enrichment$gene_lists,
  de_results_df = enrichment$results_df,
  universe_entrez = enrichment$universe,
  org_db = enrichment$org_db,
  ont_category = "BP",
  n_cores = 4
)

kegg <- generate_kegg_enrichment_plot_atac(
  gene_lists = enrichment$gene_lists,
  de_results_df = enrichment$results_df,
  universe_entrez = enrichment$universe,
  org_db = enrichment$org_db,
  organism = "hsa",
  n_cores = 4
)
```

## Key Design Principles

### 1. Standalone Functionality
All functions work independently without requiring the full reporting pipeline. Parameters are explicit rather than bundled in `report_params`.

### 2. Single Source of Truth
- **DDS/DGEList Integration**: Sample metadata, peak annotations, and QC metrics are stored directly in the data objects
- **No Separate DataFrames**: Everything travels together, reducing data sync issues

### 3. Multi-Contrast Support
Most analysis and visualization functions accept `results_list` with multiple contrasts and process them automatically.

### 4. Interactive Visualization
All plotting functions return plotly objects with hover tooltips for data exploration.

## Organism Support

### Mouse (mmu)
- Genome: mm10 (mm39 supported in annotation)
- Gene DB: org.Mm.eg.db
- Transcript DB: TxDb.Mmusculus.UCSC.mm10.knownGene
- KEGG organism: "mmu"

### Human (hsa)
- Genome: hg38 (hg19 supported in annotation)
- Gene DB: org.Hs.eg.db
- Transcript DB: TxDb.Hsapiens.UCSC.hg38.knownGene
- KEGG organism: "hsa"

## Dependencies

### Statistical Analysis
- edgeR - Differential accessibility
- DESeq2 - Data structures
- limma - Supporting functions

### Genomics
- GenomicRanges, IRanges - Interval operations
- Rsubread - Read counting
- ChIPseeker - Peak annotation
- GenomicFeatures, GenomeInfoDb - Genome databases

### Enrichment
- clusterProfiler - GO/KEGG analysis
- org.Mm.eg.db, org.Hs.eg.db - Gene annotations
- TxDb.* packages - Transcript annotations
- BiocParallel - Parallel processing

### Visualization
- ggplot2 - Static plots
- plotly - Interactive plots
- pheatmap - Heatmaps
- DT - Interactive tables
- car - PCA ellipses

### Data Manipulation
- dplyr, tidyr, purrr - Data wrangling
- stringr - String operations

See `DESCRIPTION` for the complete dependency list.

## Pipeline Integration

While atacreportR works standalone, it integrates seamlessly with:

- **app-atacreportR** (Shiny): Configure analyses interactively and generate params files
- **atacseq-Rmd** (Rmarkdown): Automated HTML reports with full QC and analysis
- **nf-core/atacseq** (Nextflow): Upstream processing from FASTQ to consensus peaks

The Rmd pipeline uses helper functions in `scripts/rmd-helpers.R` to parse params files into function arguments, but the core atacreportR functions remain independent.

## Citation

If you use atacreportR in your research, please cite:

```
University of Florida Health Cancer Center
Bioinformatics & Computational Research (BCB-SR)
atacreportR: ATAC-seq Differential Accessibility Analysis and Reporting
https://github.com/UFHCC-BCBSR/pkg-atacreportR
```

## Development

```r
# Install development version
devtools::install_github("UFHCC-BCBSR/pkg-atacreportR", ref = "dev")

# Run tests
devtools::test()

# Check package
devtools::check()

# Build documentation
devtools::document()

# Reload after changes
devtools::load_all()
```

## Troubleshooting

### Plotly plots not rendering in console
```r
# Use print() or access specific plot
print(plot_volcano(results))
# Or extract from tagList
plot_volcano(results)[[1]]
```

### Missing packages
```r
# Install organism databases as needed
BiocManager::install("org.Hs.eg.db")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
```

### BAM index files
featureCounts works without .bai files but runs faster with them. Generate with:
```bash
samtools index file.bam
```

## Related Projects

- [app-atacreportR](https://github.com/UFHCC-BCBSR/app-atacreportR) - Shiny configuration interface
- [atacseq-Rmd](https://github.com/UFHCC-BCBSR/atacseq-Rmd) - Automated HTML reports
- [nf-core/atacseq](https://github.com/nf-core/atacseq) - Processing pipeline

## Support

For questions, bug reports, or feature requests:
- Open an issue on [GitHub](https://github.com/UFHCC-BCBSR/pkg-atacreportR/issues)
- Contact UF Health Cancer Center BCB-SR

## License

[Specify your license here]
```

The updated README reflects all the refactoring work - standalone functions, explicit parameters, integrated data structures, and the separation between core package functions and Rmd helpers!
