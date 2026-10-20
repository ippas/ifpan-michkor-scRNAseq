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
  rn <- rownames(expr)
  
  # --- Extract gene symbols from rownames: ENSGxxxxx-SYMBOL
  gene_symbol_map <- gsub("^.*-", "", rn)
  
  # --- If genes_vector provided: matching by symbol
  if (!is.null(genes_vector)) {
    matched_idx <- which(gene_symbol_map %in% genes_vector)
    
    if (length(matched_idx) == 0)
      stop("❌ None of the requested gene symbols match rownames(expr). Check naming (expected format ENSG-GeneSymbol).")
    
    # --- Aggregate multiple rows with same symbol by sum
    symbols_to_keep <- gene_symbol_map[matched_idx]
    expr_keep <- expr[matched_idx, , drop = FALSE]
    
    # Summarize rows with identical gene symbols
    expr_grouped <- rowsum(as.matrix(expr_keep), group = symbols_to_keep)
    
    expr <- Matrix::Matrix(expr_grouped, sparse = TRUE)
    all_genes <- rownames(expr_grouped)
  } else {
    # no gene filtering — keep all genes with symbol extract
    # collapse duplicated symbols by sum anyway (rare but safe)
    expr_grouped <- rowsum(as.matrix(expr), group = gene_symbol_map)
    expr <- Matrix::Matrix(expr_grouped, sparse = TRUE)
    all_genes <- rownames(expr_grouped)
  }
  
  results_list <- list()
  
  for (cluster_name in clusters_to_check) {
    if (verbose) message(sprintf("\n🔹 Cluster: %s", cluster_name))
    
    cluster_cells <- rownames(meta_df[meta_df[[cluster_col]] == cluster_name, ])
    
    if (length(cluster_cells) == 0) {
      warning(sprintf("⚠️ No cells found for cluster: %s", cluster_name))
      next
    }
    
    expr_sub <- expr[, cluster_cells, drop = FALSE]
    
    # Binary detection: >= threshold
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
  
  if (verbose) message("\n✅ FAST list of detection summary tables generated successfully.")
  return(results_list)
}
