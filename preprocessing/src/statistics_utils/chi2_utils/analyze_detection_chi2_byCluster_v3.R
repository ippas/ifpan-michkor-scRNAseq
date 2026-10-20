seurat_binaryDetection_summaryByCluster_fast <- function(
    seurat_obj,
    clusters_to_check = NULL,
    genes_vector = NULL,
    assay = "RNA",
    slot = "counts",
    split.by = "treatment",
    split_levels = c("CTRL", "CORT", "DEX"),
    cluster_col = "res0.3_RNA_noScaled_brain",
    detection_threshold = 1,
    verbose = TRUE
) {
  requireNamespace("dplyr")
  requireNamespace("Matrix")
  
  meta_df <- seurat_obj@meta.data
  
  if (!cluster_col %in% colnames(meta_df))
    stop(sprintf("❌ cluster_col '%s' not found in metadata.", cluster_col))
  if (!split.by %in% colnames(meta_df))
    stop(sprintf("❌ split.by '%s' not found in metadata.", split.by))
  
  if (is.null(clusters_to_check)) {
    clusters_to_check <- unique(meta_df[[cluster_col]])
    if (verbose) message(sprintf("ℹ️ Using all clusters (%d): %s",
                                 length(clusters_to_check),
                                 paste0(clusters_to_check, collapse = ", ")))
  }
  
  Seurat::DefaultAssay(seurat_obj) <- assay
  
  expr <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = slot)
  if (!is.matrix(expr) && !inherits(expr, "Matrix")) expr <- as.matrix(expr)
  
  # --- Extract gene symbols correctly (ENSEMBL-GENE format)
  df_genes <- data.frame(original = rownames(expr), stringsAsFactors = FALSE)
  parts <- strsplit(df_genes$original, "-")
  df_genes$ensembl_id  <- sapply(parts, `[`, 1)
  df_genes$gene_symbol <- sapply(parts, function(x) paste(x[-1], collapse = "-"))
  
  # Optional filtering
  if (!is.null(genes_vector)) {
    keep_idx <- df_genes$gene_symbol %in% genes_vector
    expr     <- expr[keep_idx, , drop = FALSE]
    df_genes <- df_genes[keep_idx, , drop = FALSE]
  }
  
  # Remove duplicated gene symbols
  dup_idx <- duplicated(df_genes$gene_symbol)
  if (any(dup_idx)) {
    expr     <- expr[!dup_idx, , drop = FALSE]
    df_genes <- df_genes[!dup_idx, , drop = FALSE]
  }
  
  rownames(expr) <- df_genes$gene_symbol
  all_genes      <- df_genes$gene_symbol
  
  results_list <- list()
  
  for (cluster_name in clusters_to_check) {
    if (verbose) message(sprintf("\n🔹 Cluster: %s", cluster_name))
    
    cluster_cells <- rownames(meta_df[meta_df[[cluster_col]] == cluster_name, ])
    if (length(cluster_cells) == 0) {
      warning(sprintf("⚠️ No cells found for cluster: %s", cluster_name))
      next
    }
    
    expr_sub <- expr[, cluster_cells, drop = FALSE]
    
    det_mat <- (expr_sub >= detection_threshold) * 1L
    treatments <- meta_df[cluster_cells, split.by, drop = TRUE]
    
    df_summary <- data.frame(gene_symbol = all_genes, stringsAsFactors = FALSE)
    
    for (grp in split_levels) {
      grp_cells <- cluster_cells[treatments == grp]
      
      if (length(grp_cells) > 0) {
        det_sub <- det_mat[, colnames(det_mat) %in% grp_cells, drop = FALSE]
        detected_count_vec <- Matrix::rowSums(det_sub)
        total_cells_vec <- length(grp_cells)
        percent_vec <- round(100 * detected_count_vec / total_cells_vec, 1)
      } else {
        detected_count_vec <- integer(nrow(det_mat))
        total_cells_vec <- 0L
        percent_vec <- rep(0, nrow(det_mat))
      }
      
      df_summary[[paste0(grp, "_detected")]] <- detected_count_vec
      df_summary[[paste0(grp, "_total")]] <- total_cells_vec
      df_summary[[paste0(grp, "_percent")]] <- percent_vec
    }
    
    results_list[[cluster_name]] <- df_summary
  }
  
  if (verbose) message("\n✅ FAST binary detection summary generated.")
  return(results_list)
}
