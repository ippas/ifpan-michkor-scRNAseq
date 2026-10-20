quantile_normalize_perCluster <- function(
    count_list,         # list of matrices (e.g. from log2_transform_perCluster)
    progress = TRUE
) {
  if (!is.list(count_list)) {
    stop("Input must be a list of matrices, e.g., from log2_transform_perCluster()")
  }
  
  # Helper: manual quantile normalization of a matrix
  quantile_normalize_manual <- function(mat) {
    mat <- as.matrix(mat)
    
    # 1) sort values per sample
    mat_sorted <- apply(mat, 2, sort, na.last = TRUE)
    
    # 2) mean rank profile
    rank_means <- rowMeans(mat_sorted, na.rm = TRUE)
    
    # 3) map rank means back to original ranks
    mat_qn <- apply(mat, 2, function(col) {
      ranks <- rank(col, ties.method = "min")
      rank_means[ranks]
    })
    
    # preserve row and column names
    rownames(mat_qn) <- rownames(mat)
    colnames(mat_qn) <- colnames(mat)
    
    return(mat_qn)
  }
  
  # Output container
  out_list <- vector("list", length(count_list))
  names(out_list) <- names(count_list)
  
  # progress bar
  if (progress) pb <- utils::txtProgressBar(min = 0, max = length(count_list), style = 3)
  
  # Apply manual quantile normalization per cluster
  for (i in seq_along(count_list)) {
    m <- count_list[[i]]
    m_qn <- quantile_normalize_manual(m)
    out_list[[i]] <- m_qn
    
    if (progress) utils::setTxtProgressBar(pb, i)
  }
  
  if (progress) close(pb)
  
  return(out_list)
}

