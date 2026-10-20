aggregate_detection_by_cluster_treatment_sample <- function(
    seurat_obj,
    cluster_col   = "seurat_clusters",
    treatment_col = "treatment",
    sample_col    = "sample_id",
    assay         = "RNA",
    slot          = "counts",
    detection_threshold = 1,
    genes_keep    = NULL,
    verbose       = TRUE
) {
  
  expr <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = slot)
  if (!is.matrix(expr) && !inherits(expr, "Matrix")) expr <- as.matrix(expr)
  
  meta <- seurat_obj@meta.data
  cell_ids <- colnames(expr)
  meta <- meta[cell_ids, , drop = FALSE]
  
  # parse ENSEMBL-GENE format
  df_genes <- data.frame(original = rownames(expr), stringsAsFactors = FALSE)
  parts <- strsplit(df_genes$original, "-")
  df_genes$gene_symbol <- sapply(parts, function(x) paste(x[-1], collapse = "-"))
  
  if (!is.null(genes_keep)) {
    keep_idx <- df_genes$gene_symbol %in% genes_keep
    expr <- expr[keep_idx, , drop = FALSE]
    df_genes <- df_genes[keep_idx, , drop = FALSE]
  }
  
  rownames(expr) <- df_genes$gene_symbol
  
  # binary detection
  expr_bin <- (expr >= detection_threshold) * 1L
  
  # grouping
  grp <- interaction(
    meta[[cluster_col]],
    meta[[treatment_col]],
    meta[[sample_col]],
    drop = TRUE
  )
  
  M <- stats::model.matrix(~ grp - 1)
  colnames(M) <- sub("^grp", "", colnames(M), perl = TRUE)
  
  detected_cells <- expr_bin %*% M
  total_cells    <- Matrix::colSums(M)
  
  grouping_info <- do.call(rbind, strsplit(colnames(M), "\\."))
  colnames(grouping_info) <- c("cluster", "treatment", "sample")
  grouping_info <- as.data.frame(grouping_info, stringsAsFactors = FALSE)
  
  # RETURN LIST OF DATA.FRAMES PER CLUSTER
  clusters_unique <- unique(grouping_info$cluster)
  out_list <- setNames(vector("list", length(clusters_unique)), clusters_unique)
  
  for (cl in clusters_unique) {
    if (verbose) message("processing cluster: ", cl)
    
    idx <- grouping_info$cluster == cl
    col_ids <- which(idx)
    
    sub_detected <- detected_cells[, col_ids, drop = FALSE]
    sub_total    <- total_cells[col_ids]
    sub_treat    <- grouping_info$treatment[idx]
    sub_sample   <- grouping_info$sample[idx]
    
    # build output table for this cluster
    df_out <- data.frame(
      gene = rownames(sub_detected),
      stringsAsFactors = FALSE
    )
    
    for (j in seq_along(col_ids)) {
      name_detected <- paste0(sub_treat[j], "_", sub_sample[j], "_detected")
      name_total    <- paste0(sub_treat[j], "_", sub_sample[j], "_total")
      
      df_out[[name_detected]] <- sub_detected[, j]
      df_out[[name_total]]    <- sub_total[j]
    }
    
    out_list[[cl]] <- df_out
  }
  
  return(out_list)
}
