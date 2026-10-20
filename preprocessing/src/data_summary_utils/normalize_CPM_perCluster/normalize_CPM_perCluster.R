normalize_CPM_perCluster <- function(
    count_list,
    scale_factor = 1e6,
    progress = TRUE
) {
  if (!is.list(count_list)) {
    stop("Input must be a list of matrices (e.g. from aggregate_counts_by_cluster_sample())")
  }
  
  out_list <- vector("list", length(count_list))
  names(out_list) <- names(count_list)
  
  if (progress) pb <- utils::txtProgressBar(min = 0, max = length(count_list), style = 3)
  
  for (i in seq_along(count_list)) {
    m <- count_list[[i]]
    if (!is.matrix(m) && !inherits(m, "Matrix")) {
      m <- as.matrix(m)
    }
    
    # library sizes (sum per sample)
    lib_sizes <- colSums(m)
    
    # CPM = counts / library_size * scale
    m_cpm <- sweep(m, 2, lib_sizes, "/") * scale_factor
    
    out_list[[i]] <- m_cpm
    
    if (progress) utils::setTxtProgressBar(pb, i)
  }
  
  if (progress) close(pb)
  
  return(out_list)
}
