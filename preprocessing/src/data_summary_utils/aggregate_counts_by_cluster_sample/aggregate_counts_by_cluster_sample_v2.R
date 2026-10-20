aggregate_counts_by_cluster_sample <- function(
    seurat_obj,
    cluster_col   = "seurat_clusters",
    sample_col    = "sample_id",
    assay         = "RNA",
    layer         = "counts",
    agg           = c("sum", "mean"),
    return        = c("list", "long"),
    select_samples = NULL,
    genes_keep     = NULL,     # NEW: vector of gene symbols to keep
    progress       = TRUE
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
  if (!(cluster_col %in% colnames(meta))) stop(paste("Column", cluster_col, "not found in meta.data"))
  if (!(sample_col  %in% colnames(meta))) stop(paste("Column", sample_col,  "not found in meta.data"))
  
  # Keep meta order
  cell_ids <- colnames(expr)
  meta     <- meta[cell_ids, , drop = FALSE]
  
  # --- 3) Filter samples if needed ---
  if (!is.null(select_samples)) {
    keep_idx <- meta[[sample_col]] %in% select_samples
    expr <- expr[, keep_idx, drop = FALSE]
    meta <- meta[keep_idx, , drop = FALSE]
  }
  
  # --- 4) Extract gene symbols from rownames (ENSEMBL-GENE format) ---
  df_genes <- data.frame(
    original = rownames(expr),
    stringsAsFactors = FALSE
  )
  
  # split ENSG-GENE
  parts <- strsplit(df_genes$original, "-")
  df_genes$ensembl_id <- sapply(parts, `[`, 1)
  df_genes$gene_symbol <- sapply(parts, function(x) paste(x[-1], collapse = "-"))
  
  # --- 5) Filter genes if list provided (based on gene symbols) ---
  if (!is.null(genes_keep)) {
    keep_gene_idx <- df_genes$gene_symbol %in% genes_keep
    expr <- expr[keep_gene_idx, , drop = FALSE]
    df_genes <- df_genes[keep_gene_idx, , drop = FALSE]
  }
  
  # --- 6) Remove duplicated gene symbols (keep first) ---
  dup_idx <- duplicated(df_genes$gene_symbol)
  if (any(dup_idx)) {
    expr <- expr[!dup_idx, , drop = FALSE]
    df_genes <- df_genes[!dup_idx, , drop = FALSE]
  }
  
  # Rename matrix rows to gene symbols
  rownames(expr) <- df_genes$gene_symbol
  
  # --- 7) Create (cluster × sample) groups ---
  clusters <- as.factor(meta[[cluster_col]])
  samples  <- as.factor(meta[[sample_col]])
  grp      <- interaction(clusters, samples, drop = TRUE)
  
  M <- stats::model.matrix(~ grp - 1)
  colnames(M) <- sub("^grp", "", colnames(M), perl = TRUE)
  
  # --- 8) Aggregate counts using matrix multiplication ---
  summed <- expr %*% M
  
  if (agg == "mean") {
    n_cells_per_grp <- Matrix::colSums(M)
    summed <- sweep(summed, 2, n_cells_per_grp, "/")
  }
  
  # --- 9) Split result by cluster into list of matrices ---
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
    mat_cl <- mat_cl[, order(colnames(mat_cl)), drop = FALSE]
    
    out_list[[i]] <- mat_cl
    if (progress) utils::setTxtProgressBar(pb, i)
  }
  
  if (progress) close(pb)
  
  # --- 10) Return output ---
  if (return == "list") return(out_list)
  
  long_list <- lapply(names(out_list), function(nm) {
    cl <- sub("^cluster_", "", nm)
    m  <- out_list[[nm]]
    expand.grid(
      gene    = rownames(m),
      cluster = cl,
      sample  = colnames(m),
      stringsAsFactors = FALSE
    ) %>%
      dplyr::mutate(value = as.numeric(m))
  })
  
  do.call(rbind, long_list)
}
