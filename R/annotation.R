#' Generate peak annotation using ChIPseeker
#'
#' Annotates consensus ATAC-seq peaks with genomic features and gene information
#' using ChIPseeker and organism-specific annotation databases.
#'
#' @param dds A DESeqDataSet object containing consensus peaks in rowRanges
#' @param organism Character string specifying organism: "mmu" (mouse) or "hsa" (human)
#'
#' @return Data frame containing peak annotations including genomic location,
#'   nearest gene, distance to TSS, and gene symbols
#'
#' @importFrom ChIPseeker annotatePeak
#' @importFrom GenomicFeatures genes
#' @importFrom GenomicRanges seqnames start end width strand
#' @importFrom GenomeInfoDb seqlevels seqlevels<-
#' @importFrom AnnotationDbi mapIds
#' @importFrom SummarizedExperiment rowRanges
#' @export
generate_peak_annotation <- function(dds, organism) {
  cat("Generating peak annotation with ChIPseeker for organism:", organism, "\n")

  consensus_peaks <- rowRanges(dds)

  if (organism == "mmu") {
    requireNamespace("TxDb.Mmusculus.UCSC.mm10.knownGene", quietly = TRUE)
    requireNamespace("org.Mm.eg.db", quietly = TRUE)
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
    orgdb <- org.Mm.eg.db::org.Mm.eg.db
  } else if (organism == "hsa") {
    requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)
    requireNamespace("org.Hs.eg.db", quietly = TRUE)
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    orgdb <- org.Hs.eg.db::org.Hs.eg.db
  } else {
    stop("Unsupported organism: ", organism, ". Use 'mmu' or 'hsa'")
  }

  cat("Standardizing chromosome names...\n")
  current_chroms <- seqlevels(consensus_peaks)
  cat("Current chromosome names:", paste(head(current_chroms), collapse=", "), "\n")

  if (!any(grepl("^chr", current_chroms))) {
    seqlevels(consensus_peaks) <- paste0("chr", seqlevels(consensus_peaks))
    cat("Added 'chr' prefix to chromosome names\n")
  }

  cat("Annotating", length(consensus_peaks), "consensus peaks...\n")
  peak_anno <- annotatePeak(
    peak = consensus_peaks,
    tssRegion = c(-1000, 1000),
    TxDb = txdb,
    verbose = TRUE
  )

  anno_df <- as.data.frame(peak_anno)

  # Restore missing peaks with minimal annotation
  annotated_indices <- as.numeric(rownames(anno_df))
  all_indices <- 1:length(consensus_peaks)
  missing_indices <- setdiff(all_indices, annotated_indices)

  if(length(missing_indices) > 0) {
    missing_ranges <- consensus_peaks[missing_indices]
    missing_df <- data.frame(
      seqnames = as.character(seqnames(missing_ranges)),
      start = start(missing_ranges),
      end = end(missing_ranges),
      width = width(missing_ranges),
      strand = as.character(strand(missing_ranges)),
      annotation = "No_annotation_available",
      geneChr = as.character(seqnames(missing_ranges)),
      geneStart = start(missing_ranges),
      geneEnd = end(missing_ranges),
      geneLength = width(missing_ranges),
      geneStrand = as.character(strand(missing_ranges)),
      geneId = NA,
      distanceToTSS = NA,
      row.names = as.character(missing_indices)
    )

    missing_cols <- setdiff(colnames(anno_df), colnames(missing_df))
    missing_df[missing_cols] <- NA
    anno_df <- rbind(anno_df, missing_df)
    anno_df <- anno_df[order(as.numeric(rownames(anno_df))), ]
  }

  cat("Adding gene symbol mapping...\n")
  anno_df$Gene.Name <- mapIds(
    orgdb,
    keys = anno_df$geneId,
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"
  )

  anno_df$Entrez.ID <- anno_df$geneId
  anno_df$interval <- names(consensus_peaks)
  anno_df <- anno_df[, c("interval", setdiff(names(anno_df), "interval"))]
  anno_df$Annotation_short <- gsub(" .*", "", anno_df$annotation)

  cat("Annotation complete:", nrow(anno_df), "peaks annotated\n")
  cat("Gene symbols added for", sum(!is.na(anno_df$Gene.Name)), "peaks\n")

  return(anno_df)
}

#' Load existing peak annotation from HOMER format file
#'
#' Reads and parses a peak annotation file in HOMER format, handling
#' multiple header lines and complex formatting.
#'
#' @param annotation_file_path Character string specifying path to annotation file
#'
#' @return Data frame containing peak annotations from the file
#'
#' @export
load_peak_annotation <- function(annotation_file_path) {
  cat("Loading existing peak annotation from:", basename(annotation_file_path), "\n")

  if (!file.exists(annotation_file_path)) {
    stop("Peak annotation file not found: ", annotation_file_path)
  }

  lines <- readLines(annotation_file_path)

  header_line <- 1
  for (i in 1:min(5, length(lines))) {
    if (grepl("Chr.*Start.*End", lines[i])) {
      header_line <- i
      break
    }
  }

  cat("Found header at line", header_line, "\n")

  peaks_anno <- read.delim(
    annotation_file_path,
    header = TRUE,
    skip = header_line - 1,
    stringsAsFactors = FALSE,
    sep = "\t"
  )

  colnames(peaks_anno) <- gsub("\\.\\..*$", "", colnames(peaks_anno))
  colnames(peaks_anno) <- gsub("\\.$", "", colnames(peaks_anno))
  colnames(peaks_anno)[1] <- "interval"

  if (!"Annotation_short" %in% colnames(peaks_anno) && "Annotation" %in% colnames(peaks_anno)) {
    peaks_anno$Annotation_short <- gsub(" .*", "", peaks_anno$Annotation)
  }

  cat("Loaded annotation for", nrow(peaks_anno), "peaks\n")
  cat("First few column names:", paste(head(colnames(peaks_anno), 10), collapse = ", "), "\n")

  return(peaks_anno)
}
