#' Run differential accessibility analysis across multiple contrasts
#'
#' This function performs edgeR-based differential accessibility analysis
#' for ATAC-seq consensus peaks across multiple experimental contrasts.
#'
#' @param dge A DGEList object containing normalized count data. Must have
#'   dge$genes populated with peak annotations including an 'interval' column.
#' @param contrast_strings Character vector of contrast specifications (e.g., "1. Heat vs Control")
#' @param min_count Minimum count threshold for filterByExpr (default: 10)
#' @param min_prop Minimum proportion threshold for filterByExpr (default: 0.7)
#'
#' @return A named list where each element contains the differential accessibility
#'   results (topTags output) for one contrast
#'
#' @importFrom edgeR filterByExpr estimateDisp glmQLFit glmQLFTest topTags
#' @importFrom stats model.matrix
#' @export
run_differential_analysis <- function(dge,
                                      contrast_strings,
                                      min_count = 10,
                                      min_prop = 0.7) {
  # Avoid NSE warnings
  FDR <- logFC <- NULL

  # Validate DGEList has annotations
  if (is.null(dge$genes)) {
    stop("dge$genes must contain peak annotations. Add annotations before calling this function.")
  }

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

  # Parse and clean contrast strings
  contrast_list <- .parse_contrast_strings(contrast_strings)

  # Ensure sample names are syntactically valid
  rownames(dge$samples) <- colnames(dge) <- make.names(colnames(dge))

  # Run differential analysis for each contrast
  results_list <- lapply(seq_along(contrast_list), function(i) {
    .run_single_contrast(dge, contrast_list[[i]], min_count, min_prop)
  })

  names(results_list) <- sapply(contrast_list, function(x) {
    paste0(x[1], "_vs_", x[2])
  })

  return(results_list)
}

#' Parse contrast string specifications
#' @keywords internal
.parse_contrast_strings <- function(contrast_strings) {
  lapply(contrast_strings, function(x) {
    x_clean <- sub("^\\d+\\.\\s*", "", x)
    parts <- trimws(unlist(strsplit(x_clean, "\\s+vs\\s+", perl = TRUE)))
    gsub("-", "_", make.names(parts))
  })
}

#' Run differential analysis for a single contrast
#' @keywords internal
.run_single_contrast <- function(dge, contrast_pair, min_count, min_prop) {
  group1 <- contrast_pair[1]
  group2 <- contrast_pair[2]

  # Find which sample metadata column contains both groups
  var_match <- .find_contrast_variable(dge$samples, group1, group2)

  if (is.null(var_match)) {
    stop(paste("Could not match contrast:", group1, "vs", group2))
  }

  # Subset to relevant samples
  dge$samples$Condition <- make.names(as.character(dge$samples[[var_match]]))
  relevant_samples <- which(dge$samples$Condition %in% c(group1, group2))
  dge_contrast <- dge[, relevant_samples]
  dge_contrast$samples$Condition <- droplevels(factor(dge_contrast$samples$Condition))

  # Create design matrix
  design_contrast <- model.matrix(~0 + Condition, data = dge_contrast$samples)
  colnames(design_contrast) <- levels(dge_contrast$samples$Condition)

  # Filter lowly expressed peaks
  keep <- filterByExpr(dge_contrast, design_contrast,
                       min.count = min_count,
                       min.prop = min_prop)
  dge_contrast <- dge_contrast[keep, , keep.lib.sizes = TRUE]

  # Estimate dispersion and fit model
  dge_contrast <- estimateDisp(dge_contrast, design_contrast)
  fit_contrast <- glmQLFit(dge_contrast, design_contrast)

  # Set up contrast vector
  contrast_vec <- numeric(ncol(design_contrast))
  contrast_vec[which(colnames(design_contrast) == group1)] <- 1
  contrast_vec[which(colnames(design_contrast) == group2)] <- -1

  # Run test
  result <- glmQLFTest(fit_contrast, contrast = contrast_vec)
  topTags(result, n = Inf)
}

#' Find which metadata column contains the contrast groups
#' @keywords internal
.find_contrast_variable <- function(sample_metadata, group1, group2) {
  for (var in colnames(sample_metadata)) {
    vals <- make.names(unique(as.character(sample_metadata[[var]])))
    if (all(c(group1, group2) %in% vals)) {
      return(var)
    }
  }
  return(NULL)
}
