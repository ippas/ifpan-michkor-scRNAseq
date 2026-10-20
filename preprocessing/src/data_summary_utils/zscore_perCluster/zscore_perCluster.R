zscore_perCluster <- function(
    count_list,          # list of matrices (genes x samples)
    center = TRUE,       # subtract mean
    scale  = TRUE,       # divide by sd
    clip_range = NULL,   # e.g. 2 -> clip Z-scores to [-2,2]; NULL = no clipping
    progress = TRUE      # show progress bar
) {
  
  # --- Check input ---
  if (!is.list(count_list)) {
    stop("Input must be a list of matrices, e.g., from quantile_normalize_perCluster()")
  }
  
  out_list <- vector("list", length(count_list))
  names(out_list) <- names(count_list)
  
  if (progress) pb <- utils::txtProgressBar(min = 0, max = length(count_list), style = 3)
  
  # --- Loop through clusters ---
  for (i in seq_along(count_list)) {
    m <- as.matrix(count_list[[i]])
    
    # --- Compute row-wise mean and sd ---
    means <- rowMeans(m, na.rm = TRUE)
    sds   <- apply(m, 1, sd, na.rm = TRUE)
    
    # Avoid division by zero (genes with constant expression)
    sds[sds == 0] <- 1
    
    # --- Apply Z-score normalization manually ---
    m_z <- m
    if (center) m_z <- sweep(m_z, 1, means, "-")
    if (scale)  m_z <- sweep(m_z, 1, sds, "/")
    
    # --- Optional clipping ---
    if (!is.null(clip_range)) {
      m_z[m_z >  clip_range] <-  clip_range
      m_z[m_z < -clip_range] <- -clip_range
    }
    
    # --- Restore dimnames ---
    rownames(m_z) <- rownames(m)
    colnames(m_z) <- colnames(m)
    
    out_list[[i]] <- m_z
    
    if (progress) utils::setTxtProgressBar(pb, i)
  }
  
  if (progress) close(pb)
  
  return(out_list)
}