# atacreportR

> ATAC-seq Differential Peak Analysis and Reporting Functions

A comprehensive R package providing all data processing, statistical analysis, and visualization functions for ATAC-seq differential accessibility analysis. Includes peak processing, count normalization, edgeR-based statistical testing, ChIPseeker annotation, and interactive reporting capabilities. Required for the UFHCC BCB-SR ATAC-seq analysis pipeline.

Currently, some functions in this package can **only be run on HiPerGator** due to interaction with the /blue and /orange filesystem for use with WashU genome browser

## Features

- **Flexible Data Loading**: Works with nf-core/atacseq output, existing DDS objects, or custom peak/BAM files
- **Consensus Peak Creation**: Merges peaks across samples with reproducibility filtering
- **Multiple Normalization Methods**: TMM normalization or spike-in normalization
- **Differential Analysis**: edgeR-based quasi-likelihood framework for robust statistical testing
- **Peak Annotation**: Annotates peaks with nearby genes using ChIPseeker
- **Interactive Visualization**: PCA plots, MA plots, volcano plots with hover tooltips
- **Enrichment Analysis**: GO term and KEGG pathway analysis with gene symbol mapping
- **Genome Browser Integration**: Creates WashU Browser tracks for visualization
- **Quality Control**: Comprehensive QC metrics including FRiP scores and alignment statistics

## Installation

```r
# Install from GitHub
remotes::install_github("UFHCC-BCBSR/pkg-atacreportR")

# Or with devtools
devtools::install_github("UFHCC-BCBSR/pkg-atacreportR")
```

## Quick Start

```r
library(atacreportR)

# Parse analysis parameters
params <- parse_params("path/to/params.txt")

# Load and prepare data
analysis_data <- prepare_analysis_data(params)

# Run differential analysis
contrast_list <- parse_contrasts(params$contrasts)
contrast_strings <- sapply(seq_along(contrast_list), function(i) {
  paste0(i, ". ", contrast_list[[i]][1], " vs ", contrast_list[[i]][2])
})

results <- run_differential_analysis(
  dge = dge,
  peaks_anno = analysis_data$peaks_anno,
  contrast_strings = contrast_strings,
  report_params = params
)

# Generate visualizations
pca_plot <- plot_pca(dge, title = "Sample PCA", plot_var = "Condition")
qc_table <- summarize_atac_sample_qc(analysis_data)
```

## Main Functions

### Data Preparation

- `prepare_analysis_data()` - Main data loading and preparation function
- `parse_params()` - Parse parameter files
- `load_sample_metadata()` - Load and validate sample sheets
- `load_dds_file()` - Load existing DESeqDataSet objects
- `load_qc_data()` - Load FRiP scores and flagstat files

### Peak Processing

- `make_consensus_peaks()` - Create consensus peaks from multiple samples
- `count_bam_reads_in_peaks()` - Quantify reads in peaks using featureCounts
- `create_dds_from_peaks()` - Build DESeqDataSet from peak and BAM files

### Annotation

- `generate_peak_annotation()` - Annotate peaks with ChIPseeker
- `load_peak_annotation()` - Load existing HOMER-format annotations

### Differential Analysis

- `run_differential_analysis()` - Perform edgeR differential accessibility analysis
- `parse_contrasts()` - Parse contrast specifications

### Visualization

- `plot_pca()` - Interactive PCA plot with group ellipses
- `summarize_atac_sample_qc()` - Comprehensive QC summary table
- `preprocess_result()` - Prepare results for MA/volcano plots
- `download_button_png()` - Create download buttons for plots

### Enrichment Analysis

- `generate_enrichment_plot_atac()` - GO term enrichment with custom tooltips
- `generate_kegg_enrichment_plot_atac()` - KEGG pathway enrichment analysis

### Utilities

- `clean_group_name()` - Remove leading X from group names
- `copy_file_to_web()` - Copy files for web hosting
- `create_da_bigbed()` - Create BigBed tracks for genome browsers
- `validate_data_consistency()` - Validate data integrity

## Supported Input Formats

### nf-core/atacseq Output

```r
params$dds_file <- "path/to/consensus_peaks.mLb.clN.dds.RData"
params$peak_annotation <- "path/to/consensus_peaks.mLb.clN.annotatePeaks.txt"
analysis_data <- prepare_analysis_data(params)
```

### Custom Peak/BAM Files

```r
params$peak_files <- "sample1:/path/to/peaks1.narrowPeak,sample2:/path/to/peaks2.narrowPeak"
params$bam_files <- "sample1:/path/to/sample1.bam,sample2:/path/to/sample2.bam"
analysis_data <- prepare_analysis_data(params)
```

### Existing DDS Object

```r
params$dds_file <- "path/to/my_dds.RData"
analysis_data <- prepare_analysis_data(params)
```

## Analysis Workflow Example

Here's a complete workflow from start to finish:

```r
library(atacreportR)
library(edgeR)
library(dplyr)

# 1. Parse parameters
params <- parse_params("analysis_params.txt")

# 2. Load and prepare data
analysis_data <- prepare_analysis_data(params)

# 3. Create DGEList with normalization
targetGenome_counts <- SummarizedExperiment::assay(analysis_data$dds, "counts")

# TMM normalization
dge <- DGEList(counts = targetGenome_counts)
dge$samples <- cbind(dge$samples, as.data.frame(SummarizedExperiment::colData(analysis_data$dds)))
dge <- calcNormFactors(dge, method = "TMM")

# 4. Run differential analysis
contrast_list <- parse_contrasts(params$contrasts)
contrast_strings <- sapply(seq_along(contrast_list), function(i) {
  paste0(i, ". ", contrast_list[[i]][1], " vs ", contrast_list[[i]][2])
})

results_list <- run_differential_analysis(
  dge = dge,
  peaks_anno = analysis_data$peaks_anno,
  contrast_strings = contrast_strings,
  report_params = params
)

# 5. Extract significant peaks
sig_peaks <- results_list[[1]]$table %>%
  filter(FDR < 0.05, abs(logFC) > 1)

# 6. Run enrichment analysis
# Prepare gene lists
gene_lists <- list()
for (contrast in names(results_list)) {
  df <- results_list[[contrast]]$table
  up_genes <- df %>%
    filter(FDR < 0.05, logFC > 1) %>%
    pull(ENTREZID) %>%
    na.omit() %>% unique()
  
  gene_lists[[paste0(contrast, ".up")]] <- up_genes
}

# Get gene universe
de_results_df <- bind_rows(lapply(names(results_list), function(contrast) {
  df <- results_list[[contrast]]$table
  df$Contrast <- contrast
  return(df)
}))

universe_entrez <- na.omit(unique(de_results_df$ENTREZID))

# Run GO enrichment
GO_results <- generate_enrichment_plot_atac(
  gene_lists = gene_lists,
  de_results_df = de_results_df,
  universe_entrez = universe_entrez,
  ont_category = "BP",
  annotation_db = params$annotation_db
)

# 7. Create visualizations
qc_table <- summarize_atac_sample_qc(analysis_data)
pca_plot <- plot_pca(dge, title = "PCA", plot_var = params$variable)
```

## Organism Support

- **Mouse (mmu)**: Uses mm10 genome, org.Mm.eg.db, and TxDb.Mmusculus.UCSC.mm10.knownGene
- **Human (hsa)**: Uses hg38 genome, org.Hs.eg.db, and TxDb.Hsapiens.UCSC.hg38.knownGene

Genome chromosome size files are bundled with the package for offline use.

## Dependencies

### Core Dependencies
- edgeR, DESeq2 - Differential analysis
- GenomicRanges, IRanges - Genomic interval operations
- ChIPseeker - Peak annotation
- clusterProfiler - Enrichment analysis

### Visualization
- ggplot2, plotly - Interactive plots
- DT - Interactive tables
- pheatmap - Heatmaps

### Data Manipulation
- dplyr, tidyr, purrr - Data wrangling
- stringr - String operations

See `DESCRIPTION` for complete list.

## Citation

If you use atacreportR in your research, please cite:

```
University of Florida Health Cancer Center
Bioinformatics & Computational Research (BCB-SR)
atacreportR: ATAC-seq Differential Peak Analysis and Reporting Functions
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
```

## Related Projects

- [app-atacreportR](https://github.com/UFHCC-BCBSR/app-atacreportR) - Shiny application for interactive analysis configuration
- [nf-core/atacseq](https://github.com/nf-core/atacseq) - Upstream ATAC-seq processing pipeline

## Support

For questions, bug reports, or feature requests:
- Open an issue on [GitHub](https://github.com/UFHCC-BCBSR/pkg-atacreportR/issues)
- Contact UF Health Cancer Center BCB-SR
