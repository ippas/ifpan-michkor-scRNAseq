seurat_binaryDetection_summaryTable <- function(
    seurat_obj,
    genes_vector,
    assay = "RNA",
    slot = "counts",
    split.by = "treatment",
    split_levels = c("CTRL", "CORT", "DEX"),
    detection_threshold = 0,
    verbose = TRUE
) {
  requireNamespace("dplyr")
  
  # 🔹 Przygotowanie wynikowej listy
  results_list <- list()
  
  # 🔹 Iteracja po genach
  for (gene_pattern in genes_vector) {
    if (verbose) message(sprintf("▶️ Processing gene: %s", gene_pattern))
    
    all_features <- rownames(GetAssayData(seurat_obj, assay = assay, slot = slot))
    matched <- grep(paste0("-", gene_pattern, "$"), all_features, value = TRUE)
    if (length(matched) == 0) {
      warning(sprintf("❌ Gene not found: %s", gene_pattern))
      next
    }
    if (length(matched) > 1) {
      message("⚠️ Multiple matches found for ", gene_pattern, ", using first: ", matched[1])
      matched <- matched[1]
    }
    
    # 🔹 Ekspresja i detekcja
    expr_values <- GetAssayData(seurat_obj, assay = assay, slot = slot)[matched, ]
    detected_binary <- as.numeric(expr_values > detection_threshold)
    
    # 🔹 Dodaj metadane z grupami
    meta_df <- seurat_obj@meta.data
    if (!split.by %in% colnames(meta_df))
      stop(sprintf("❌ Column '%s' not found in metadata", split.by))
    
    df <- data.frame(
      gene_symbol = gene_pattern,
      treatment = meta_df[[split.by]],
      detected = detected_binary
    )
    
    # 🔹 Agregacja po treatment
    summary_df <- df |>
      dplyr::group_by(treatment) |>
      dplyr::summarise(
        detected_count = sum(detected, na.rm = TRUE),
        total_cells = dplyr::n(),
        percent = round(100 * detected_count / total_cells, 1),
        .groups = "drop"
      )
    
    # 🔹 Uporządkuj kolejność grup
    summary_df$treatment <- factor(summary_df$treatment, levels = split_levels)
    summary_df <- summary_df |> dplyr::arrange(treatment)
    
    # 🔹 Rozszerz w jedną linię
    row_out <- data.frame(gene_symbol = gene_pattern)
    for (grp in split_levels) {
      if (grp %in% summary_df$treatment) {
        row_vals <- summary_df |> dplyr::filter(treatment == grp)
        row_out[[paste0(grp, "_detected")]] <- row_vals$detected_count
        row_out[[paste0(grp, "_total")]] <- row_vals$total_cells
        row_out[[paste0(grp, "_percent")]] <- row_vals$percent
      } else {
        row_out[[paste0(grp, "_detected")]] <- 0
        row_out[[paste0(grp, "_total")]] <- 0
        row_out[[paste0(grp, "_percent")]] <- 0
      }
    }
    
    results_list[[gene_pattern]] <- row_out
  }
  
  # 🔹 Połącz w jedną tabelę
  result_df <- dplyr::bind_rows(results_list)
  
  if (verbose) message("✅ Summary table generated successfully.")
  
  return(result_df)
}