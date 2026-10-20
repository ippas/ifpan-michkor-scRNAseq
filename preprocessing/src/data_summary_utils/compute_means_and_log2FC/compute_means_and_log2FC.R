compute_means_and_log2FC <- function(
    matrix_list,
    sample_map
) {
  if (!is.list(matrix_list)) {
    stop("matrix_list must be a list of matrices (e.g. aggregated counts per cluster)")
  }
  
  out_list <- vector("list", length(matrix_list))
  names(out_list) <- names(matrix_list)
  
  for (i in seq_along(matrix_list)) {
    m <- matrix_list[[i]]
    
    if (!is.matrix(m) && !inherits(m, "Matrix")) {
      m <- as.matrix(m)
    }
    
    sm <- sample_map[match(colnames(m), sample_map$bc1_well), ]
    
    ctrl_idx <- which(sm$treatment == "CTRL")
    cort_idx <- which(sm$treatment == "CORT")
    dex_idx  <- which(sm$treatment == "DEX")
    
    mean_ctrl <- rowMeans(m[, ctrl_idx, drop = FALSE])
    mean_cort <- rowMeans(m[, cort_idx, drop = FALSE])
    mean_dex  <- rowMeans(m[, dex_idx,  drop = FALSE])
    
    log2fc_cort <- log2(mean_cort + 1) - log2(mean_ctrl + 1)
    log2fc_dex  <- log2(mean_dex  + 1) - log2(mean_ctrl + 1)
    
    df <- data.frame(
      gene_symbol = rownames(m),      # <-- zmiana nazwy kolumny
      mean_CTRL   = mean_ctrl,
      mean_CORT   = mean_cort,
      mean_DEX    = mean_dex,
      log2FC_CORT = log2fc_cort,
      log2FC_DEX  = log2fc_dex,
      row.names   = NULL,            # <-- usuwa nazwy wierszy
      stringsAsFactors = FALSE
    )
    
    out_list[[i]] <- df
  }
  
  return(out_list)
}
