#' Run differential accessibility analysis across multiple contrasts
#'
#' This function performs edgeR-based differential accessibility analysis
#' for ATAC-seq consensus peaks across multiple experimental contrasts.
#'
#' @param dge A DGEList object containing normalized count data
#' @param peaks_anno Data frame containing peak annotations
#' @param contrast_strings Character vector of contrast specifications (e.g., "1. Heat vs Control")
#' @param report_params List of analysis parameters including min_count_for_filtering and min_prop_for_filtering
#'
#' @return A named list where each element contains the differential accessibility
#'   results (topTags output) for one contrast
#'
#' @importFrom edgeR filterByExpr estimateDisp glmQLFit glmQLFTest
#' @importFrom stats model.matrix
#' @export
run_differential_analysis <- function(dge, peaks_anno, contrast_strings, report_params) {
  # Avoid NSE warnings
  FDR <- logFC <- NULL

  # Annotate peaks
  row.names(peaks_anno) <- peaks_anno$interval
  dge$genes <- peaks_anno[rownames(dge), ]

  # Convert NULL values in Gene.Name to NA_character_
  if (!is.null(dge$genes$Gene.Name)) {
    if (is.list(dge$genes$Gene.Name)) {
      dge$genes$Gene.Name <- sapply(dge$genes$Gene.Name, function(x) {
        if (is.null(x) || length(x) == 0) NA_character_ else as.character(x[1])
      })
    } else {
      dge$genes$Gene.Name <- as.character(dge$genes$Gene.Name)
    }
  } else {
    dge$genes$Gene.Name <- rep(NA_character_, nrow(dge$genes))
  }

  # Run differential analysis per contrast
  contrast_list <- lapply(contrast_strings, function(x) {
    x_clean <- sub("^\\d+\\.\\s*", "", x)
    parts <- trimws(unlist(strsplit(x_clean, "\\s+vs\\s+", perl = TRUE)))
    make.names(parts)
  })

  rownames(dge$samples) <- colnames(dge) <- make.names(colnames(dge))

  results_list <- list()
  contrast_list <- lapply(contrast_list, function(x) gsub("-", "_", x))

  for (i in seq_along(contrast_list)) {
    group1 <- make.names(contrast_list[[i]][1])
    group2 <- make.names(contrast_list[[i]][2])
    contrast_name <- paste0(group1, "_vs_", group2)

    var_match <- NULL
    for (var in colnames(dge$samples)) {
      vals <- make.names(unique(as.character(dge$samples[[var]])))
      if (all(c(group1, group2) %in% vals)) {
        var_match <- var
        break
      }
    }

    if (is.null(var_match)) stop(paste("Could not match contrast:", contrast_name))

    dge$samples$Condition <- make.names(as.character(dge$samples[[var_match]]))
    relevant_samples <- which(dge$samples$Condition %in% c(group1, group2))
    dge_contrast <- dge[, relevant_samples]

    dge_contrast$samples$Condition <- droplevels(factor(dge_contrast$samples$Condition))
    design_contrast <- model.matrix(~0 + Condition, data = dge_contrast$samples)
    colnames(design_contrast) <- levels(dge_contrast$samples$Condition)

    keep <- filterByExpr(dge_contrast, design_contrast,
                         min.count = report_params$min_count_for_filtering,
                         min.prop = report_params$min_prop_for_filtering)
    dge_contrast <- dge_contrast[keep, , keep.lib.sizes = TRUE]

    dge_contrast <- estimateDisp(dge_contrast, design_contrast)
    fit_contrast <- glmQLFit(dge_contrast, design_contrast)

    contrast_vec <- numeric(ncol(design_contrast))
    contrast_vec[which(colnames(design_contrast) == group1)] <- 1
    contrast_vec[which(colnames(design_contrast) == group2)] <- -1

    result <- glmQLFTest(fit_contrast, contrast = contrast_vec)
    results_list[[contrast_name]] <- topTags(result, n = Inf)
  }

  return(results_list)
}
