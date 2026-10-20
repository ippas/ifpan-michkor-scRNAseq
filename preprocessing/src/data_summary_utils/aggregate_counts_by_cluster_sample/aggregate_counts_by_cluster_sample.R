aggregate_counts_by_cluster_sample <- function(
    seurat_obj,
    cluster_col = "seurat_clusters", # column name in meta.data with cluster IDs
    sample_col  = "sample_id",       # column name in meta.data with sample IDs
    assay       = "RNA",             # assay to use
    layer       = "counts",          # Seurat v5 layer: "counts", "data", etc.
    agg         = c("sum", "mean"),  # aggregation function type
    return      = c("list", "long"), # output format
    progress    = TRUE               # show progress bar
) {
  # Match argument choices
  agg    <- match.arg(agg)
  return <- match.arg(return)
  
  # --- 1) Extract expression data ---
  expr <- Seurat::GetAssayData(seurat_obj, assay = assay, layer = layer)
  
  # Ensure expr is a matrix or sparse matrix
  if (!is.matrix(expr) && !inherits(expr, "Matrix")) {
    expr <- as.matrix(expr)
  }
  
  # --- 2) Extract metadata ---
  meta <- seurat_obj@meta.data
  if (!(cluster_col %in% colnames(meta))) {
    stop(paste("Column", cluster_col, "not found in meta.data"))
  }
  if (!(sample_col %in% colnames(meta))) {
    stop(paste("Column", sample_col, "not found in meta.data"))
  }
  
  # Keep meta in the same order as columns in expr
  cell_ids <- colnames(expr)
  meta     <- meta[cell_ids, , drop = FALSE]
  
  clusters <- as.factor(meta[[cluster_col]])
  samples  <- as.factor(meta[[sample_col]])
  
  # --- 3) Create group IDs for (cluster, sample) combinations ---
  grp <- interaction(clusters, samples, drop = TRUE) # e.g., "0.sampleA"
  grp_levels <- levels(grp)
  
  # One-hot encode groups (rows = cells, cols = groups)
  M <- stats::model.matrix(~ grp - 1)
  colnames(M) <- sub("^grp", "", colnames(M), perl = TRUE)
  
  # --- 4) Aggregate using matrix multiplication ---
  summed <- expr %*% M # genes x (cluster|sample) matrix
  
  # --- 5) If mean aggregation: divide by number of cells per group ---
  if (agg == "mean") {
    n_cells_per_grp <- Matrix::colSums(M)
    summed <- sweep(summed, 2, n_cells_per_grp, "/")
  }
  
  # --- 6) Split results by cluster ---
  split_parts <- do.call(rbind, strsplit(colnames(summed), "\\.", fixed = FALSE))
  colnames(split_parts) <- c("cluster", "sample")
  split_parts <- as.data.frame(split_parts, stringsAsFactors = FALSE)
  
  unique_clusters <- unique(split_parts$cluster)
  if (progress) {
    pb <- utils::txtProgressBar(min = 0, max = length(unique_clusters), style = 3)
  }
  
  out_list <- vector("list", length(unique_clusters))
  names(out_list) <- paste0("cluster_", unique_clusters)
  
  for (i in seq_along(unique_clusters)) {
    cl <- unique_clusters[i]
    idx <- which(split_parts$cluster == cl)
    
    mat_cl <- summed[, idx, drop = FALSE]
    samples_in_cl <- split_parts$sample[idx]
    colnames(mat_cl) <- samples_in_cl
    
    # Optional: sort samples alphabetically
    ord <- order(colnames(mat_cl))
    mat_cl <- mat_cl[, ord, drop = FALSE]
    
    out_list[[i]] <- mat_cl
    
    if (progress) utils::setTxtProgressBar(pb, i)
  }
  if (progress) close(pb)
  
  # --- 7) Return in requested format ---
  if (return == "list") {
    return(out_list)
  }
  
  # Long format: gene, cluster, sample, value
  long_list <- lapply(names(out_list), function(nm) {
    cl <- sub("^cluster_", "", nm)
    m  <- out_list[[nm]]
    data.frame(
      gene    = rownames(m)[row(m)],
      cluster = cl,
      sample  = colnames(m)[col(m)],
      value   = as.numeric(m),
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, long_list)
}



aggregate_counts_by_cluster_sample <- function(
    seurat_obj,
    cluster_col   = "seurat_clusters", # column name in meta.data with cluster IDs
    sample_col    = "sample_id",       # column name in meta.data with sample IDs
    assay         = "RNA",             # assay to use
    layer         = "counts",          # Seurat v5 layer: "counts", "data", etc.
    agg           = c("sum", "mean"),  # aggregation function type
    return        = c("list", "long"), # output format
    select_samples = NULL,             # vector of samples to keep (NULL = all)
    progress      = TRUE               # show progress bar
) {
  # Match argument choices
  agg    <- match.arg(agg)
  return <- match.arg(return)
  
  # --- 1) Extract expression data ---
  expr <- Seurat::GetAssayData(seurat_obj, assay = assay, layer = layer)
  if (!is.matrix(expr) && !inherits(expr, "Matrix")) {
    expr <- as.matrix(expr)
  }
  
  # --- 2) Extract metadata ---
  meta <- seurat_obj@meta.data
  if (!(cluster_col %in% colnames(meta))) {
    stop(paste("Column", cluster_col, "not found in meta.data"))
  }
  if (!(sample_col %in% colnames(meta))) {
    stop(paste("Column", sample_col, "not found in meta.data"))
  }
  
  # Keep meta in same order as expression
  cell_ids <- colnames(expr)
  meta     <- meta[cell_ids, , drop = FALSE]
  
  # --- NEW: filter samples if requested ---
  if (!is.null(select_samples)) {
    keep_idx <- meta[[sample_col]] %in% select_samples
    expr <- expr[, keep_idx, drop = FALSE]
    meta <- meta[keep_idx, , drop = FALSE]
  }
  
  clusters <- as.factor(meta[[cluster_col]])
  samples  <- as.factor(meta[[sample_col]])
  
  # --- 3) Create group IDs for (cluster, sample) combinations ---
  grp <- interaction(clusters, samples, drop = TRUE)
  M   <- stats::model.matrix(~ grp - 1)
  colnames(M) <- sub("^grp", "", colnames(M), perl = TRUE)
  
  # --- 4) Aggregate ---
  summed <- expr %*% M
  if (agg == "mean") {
    n_cells_per_grp <- Matrix::colSums(M)
    summed <- sweep(summed, 2, n_cells_per_grp, "/")
  }
  
  # --- 5) Split results by cluster ---
  split_parts <- do.call(rbind, strsplit(colnames(summed), "\\."))
  colnames(split_parts) <- c("cluster", "sample")
  split_parts <- as.data.frame(split_parts, stringsAsFactors = FALSE)
  
  unique_clusters <- unique(split_parts$cluster)
  if (progress) pb <- utils::txtProgressBar(min = 0, max = length(unique_clusters), style = 3)
  
  out_list <- vector("list", length(unique_clusters))
  names(out_list) <- paste0("cluster_", unique_clusters)
  
  for (i in seq_along(unique_clusters)) {
    cl <- unique_clusters[i]
    idx <- which(split_parts$cluster == cl)
    
    mat_cl <- summed[, idx, drop = FALSE]
    samples_in_cl <- split_parts$sample[idx]
    colnames(mat_cl) <- samples_in_cl
    
    ord <- order(colnames(mat_cl))
    mat_cl <- mat_cl[, ord, drop = FALSE]
    
    out_list[[i]] <- mat_cl
    if (progress) utils::setTxtProgressBar(pb, i)
  }
  if (progress) close(pb)
  
  # --- 6) Return ---
  if (return == "list") return(out_list)
  
  long_list <- lapply(names(out_list), function(nm) {
    cl <- sub("^cluster_", "", nm)
    m  <- out_list[[nm]]
    data.frame(
      gene    = rownames(m)[row(m)],
      cluster = cl,
      sample  = colnames(m)[col(m)],
      value   = as.numeric(m),
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, long_list)
}
