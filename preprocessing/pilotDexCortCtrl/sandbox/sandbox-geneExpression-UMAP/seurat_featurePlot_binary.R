seurat_featurePlot_binary <- function(
    seurat_obj,
    gene_pattern,
    assay = "RNA",
    slot = "counts",
    split.by = NULL,
    split_levels = NULL,
    cluster_col = NULL,
    reduction = "umap",
    cols = c("gray90", "firebrick"),
    point_size = 0.8,
    legend_position = "none",
    detection_threshold = 0,
    return_data = FALSE,
    # 🔹 nowe argumenty
    title_size = 16,
    facet_title_size = 14
) {
  library(ggplot2)
  library(dplyr)
  
  # create local copy
  obj_tmp <- seurat_obj
  if (!is.null(assay)) DefaultAssay(obj_tmp) <- assay
  
  # find matching genes
  all_features <- rownames(GetAssayData(obj_tmp, assay = assay, slot = slot))
  matched <- grep(paste0("-", gene_pattern, "$"), all_features, value = TRUE)
  
  if (length(matched) == 0)
    stop(paste0("❌ No gene found matching pattern: ", gene_pattern))
  if (length(matched) > 1) {
    message("⚠️ Multiple matches found, using first: ", paste(matched, collapse = ", "))
    matched <- matched[1]
  } else {
    message("✅ Found feature: ", matched)
  }
  
  # embeddings
  emb <- as.data.frame(Embeddings(obj_tmp, reduction = reduction))
  colnames(emb) <- c("dim1", "dim2")[seq_len(ncol(emb))]
  
  # expression values
  expr_values <- GetAssayData(obj_tmp, assay = assay, slot = slot)[matched, ]
  
  # binary detection
  detected_binary <- as.numeric(expr_values > detection_threshold)
  
  # combine data
  data_df <- emb %>%
    mutate(
      expression = expr_values,
      detected = detected_binary
    )
  
  # metadata
  if (!is.null(split.by)) data_df$treatment <- obj_tmp@meta.data[[split.by]]
  if (!is.null(cluster_col) && cluster_col %in% colnames(obj_tmp@meta.data))
    data_df$cluster <- obj_tmp@meta.data[[cluster_col]]
  
  # reorder split.by
  if (!is.null(split.by)) {
    current_values <- unique(data_df$treatment)
    if (is.null(split_levels)) {
      preferred_order <- c("CTRL", "CORT", "DEX")
      split_levels <- intersect(preferred_order, current_values)
      split_levels <- c(split_levels, setdiff(current_values, split_levels))
    }
    data_df$treatment <- factor(data_df$treatment, levels = split_levels)
  }
  
  # compute detection stats
  if (is.null(split.by)) {
    n_detected <- sum(data_df$detected)
    n_total <- nrow(data_df)
    percent_detected <- round(100 * n_detected / n_total, 1)
    title_text <- paste0(
      matched, " — detected in ", n_detected, "/", n_total,
      " cells (", percent_detected, "%)"
    )
  } else {
    summary_df <- as.data.frame(table(data_df$treatment))
    colnames(summary_df) <- c("treatment", "total_cells")
    
    det_summary <- as.data.frame(table(data_df$treatment[data_df$detected == 1]))
    colnames(det_summary) <- c("treatment", "detected_count")
    
    summary_df <- merge(summary_df, det_summary, by = "treatment", all.x = TRUE)
    summary_df$detected_count[is.na(summary_df$detected_count)] <- 0
    summary_df$percent <- round(100 * summary_df$detected_count / summary_df$total_cells, 1)
    
    data_df <- merge(data_df, summary_df, by = "treatment", all.x = TRUE)
    data_df$facet_label <- paste0(
      data_df$treatment, "\n",
      data_df$detected_count, "/", data_df$total_cells,
      " (", data_df$percent, "%)"
    )
    
    # zachowaj kolejność facetów zgodną z split_levels
    if (!is.null(split_levels)) {
      data_df$treatment <- factor(data_df$treatment, levels = split_levels)
      facet_labels_ordered <- data_df %>%
        distinct(treatment, facet_label) %>%
        arrange(match(treatment, split_levels)) %>%
        pull(facet_label)
      data_df$facet_label <- factor(data_df$facet_label, levels = facet_labels_ordered)
    }
    
    title_text <- paste0(matched, " (binary detection)")
  }
  
  # return data if requested
  if (return_data) return(data_df)
  
  # plotting
  p <- ggplot(data_df, aes(x = dim1, y = dim2, color = factor(detected))) +
    geom_point(size = point_size) +
    scale_color_manual(
      values = c("0" = cols[1], "1" = cols[2]),
      labels = c("Not detected", "Detected")
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = legend_position,
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(size = title_size, color = "black"),
      strip.text = element_text(size = facet_title_size, color = "black"),  # 🔹 facet labels
      text = element_text(color = "black")  # 🔹 cały tekst czarny
    ) +
    labs(
      title = title_text,
      color = "Detection"
    )
  
  if (!is.null(split.by)) {
    p <- p + facet_wrap(~facet_label) + labs(title = title_text)
  }
  
  return(p)
}

# seurat_featurePlot_binary(
#   seurat_obj = integrated_pilotDexCortCtrl_filtered,
#   gene_pattern = "SGK1",
#   split.by = "treatment",
#   cluster_col = "res0.3_RNA_noScaled_brain",
#   cols = c("gray90", "firebrick"),
#   point_size = 0.2,
#   title_size = 14,
#   facet_title_size = 14
# )
