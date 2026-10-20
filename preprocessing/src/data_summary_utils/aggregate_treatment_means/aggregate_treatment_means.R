aggregate_treatment_means <- function(
    cpm_mean_list,
    sample_map,
    sample_col   = "bc1_well",
    treatment_col = "treatment",
    progress = TRUE
) {
  # Check mapping validity
  if (!is.list(cpm_mean_list)) stop("cpm_mean_list must be a list of matrices")
  if (!all(c(sample_col, treatment_col) %in% colnames(sample_map))) {
    stop("sample_map must contain sample_col and treatment_col columns")
  }
  
  treatments <- unique(sample_map[[treatment_col]])
  
  out_list <- vector("list", length(cpm_mean_list))
  names(out_list) <- names(cpm_mean_list)
  
  if (progress) pb <- utils::txtProgressBar(min = 0, max = length(cpm_mean_list), style = 3)
  
  for (i in seq_along(cpm_mean_list)) {
    mat <- cpm_mean_list[[i]]    # genes × samples matrix
    
    # Ensure sample names match mapping
    valid_samples <- intersect(colnames(mat), sample_map[[sample_col]])
    mat <- mat[, valid_samples, drop = FALSE]
    
    # reorder rows according to matrix sample order
    map <- sample_map[match(colnames(mat), sample_map[[sample_col]]), ]
    
    # For each treatment, calculate gene-wise mean
    mean_treat <- sapply(treatments, function(t) {
      samp <- map[[sample_col]][map[[treatment_col]] == t]
      if (length(samp) == 0) return(rep(NA, nrow(mat)))
      rowMeans(mat[, samp, drop = FALSE], na.rm = TRUE)
    })
    
    colnames(mean_treat) <- paste0("mean_", treatments)
    
    out_list[[i]] <- mean_treat
    
    if (progress) utils::setTxtProgressBar(pb, i)
  }
  
  if (progress) close(pb)
  
  return(out_list)
}