seurat_binaryDetection_summaryByCluster_splitList <- function(
    seurat_obj,
    genes_vector,
    clusters_to_check,
    assay = "RNA",
    slot = "counts",
    split.by = "treatment",
    split_levels = c("CTRL", "CORT", "DEX"),
    cluster_col = "res0.3_RNA_noScaled_brain",
    detection_threshold = 0,
    verbose = TRUE
) {
  requireNamespace("dplyr")
  
  meta_df <- seurat_obj@meta.data
  if (!cluster_col %in% colnames(meta_df))
    stop(sprintf("❌ Cluster column '%s' not found in metadata.", cluster_col))
  if (!split.by %in% colnames(meta_df))
    stop(sprintf("❌ Split column '%s' not found in metadata.", split.by))
  
  # Upewnij się, że assay i slot są ustawione spójnie
  if (!is.null(assay)) Seurat::DefaultAssay(seurat_obj) <- assay
  
  results_by_cluster <- list()
  n_genes <- length(genes_vector)
  
  # Prelista wszystkich feature'ów, żeby szybko dopasowywać nazwę
  all_features <- rownames(Seurat::GetAssayData(seurat_obj, assay = assay, slot = slot))
  
  for (cluster_name in clusters_to_check) {
    if (verbose) message(sprintf("\n🔹 Cluster: %s", cluster_name))
    
    cluster_cells <- rownames(meta_df[meta_df[[cluster_col]] == cluster_name, ])
    if (length(cluster_cells) == 0) {
      warning(sprintf("⚠️ No cells found for cluster: %s", cluster_name))
      next
    }
    
    # treatment dla komórek w danym klastrze (w stałej kolejności)
    treatments_vec <- meta_df[cluster_cells, split.by, drop = TRUE]
    
    cluster_rows <- list()
    
    for (i in seq_along(genes_vector)) {
      gene_pattern <- genes_vector[i]
      percent_done <- round(100 * i / n_genes, 1)
      if (verbose) message(sprintf("   ▶️ Gene: %s (%d / %d, %s%%)", gene_pattern, i, n_genes, percent_done))
      
      # dopasowanie nazwy feature na końcu (twoja konwencja "-GENE$")
      matched <- grep(paste0("-", gene_pattern, "$"), all_features, value = TRUE)
      if (length(matched) == 0) {
        warning(sprintf("❌ Gene not found: %s", gene_pattern))
        next
      }
      if (length(matched) > 1) {
        if (verbose) message("   ⚠️ Multiple matches, using first: ", matched[1])
        matched <- matched[1]
      }
      
      # KLUCZOWA ZMIANA: FetchData gwarantuje gęsty wektor w kolejności komórek
      df_gene <- Seurat::FetchData(
        object = seurat_obj,
        vars = matched,
        cells = cluster_cells,
        slot = slot
      )
      # wektor ekspresji zgodny z kolejnością 'cluster_cells'
      expr_vec <- df_gene[[matched]]
      detected_binary <- as.integer(expr_vec > detection_threshold)
      
      df <- data.frame(
        gene_symbol = gene_pattern,
        treatment = treatments_vec,
        detected = detected_binary,
        stringsAsFactors = FALSE
      )
      
      summary_df <- df |>
        dplyr::group_by(treatment) |>
        dplyr::summarise(
          detected_count = sum(detected, na.rm = TRUE),
          total_cells    = dplyr::n(),
          percent        = ifelse(total_cells > 0,
                                  round(100 * detected_count / total_cells, 1),
                                  0),
          .groups = "drop"
        ) |>
        dplyr::mutate(treatment = factor(treatment, levels = split_levels))
      
      # szeroki format: CTRL/CORT/DEX x (detected/total/percent)
      wide <- data.frame(gene_symbol = gene_pattern, stringsAsFactors = FALSE)
      for (grp in split_levels) {
        if (grp %in% summary_df$treatment) {
          row_vals <- summary_df[summary_df$treatment == grp, , drop = FALSE]
          wide[[paste0(grp, "_detected")]] <- row_vals$detected_count
          wide[[paste0(grp, "_total")]]    <- row_vals$total_cells
          wide[[paste0(grp, "_percent")]]  <- row_vals$percent
        } else {
          wide[[paste0(grp, "_detected")]] <- 0L
          wide[[paste0(grp, "_total")]]    <- 0L
          wide[[paste0(grp, "_percent")]]  <- 0
        }
      }
      
      cluster_rows[[gene_pattern]] <- wide
    }
    
    results_by_cluster[[cluster_name]] <- dplyr::bind_rows(cluster_rows)
  }
  
  if (verbose) message("\n✅ List of data frames (per cluster) generated (fixed detection).")
  return(results_by_cluster)
}