seurat_featurePlot_multiBinary <- function(
    seurat_obj,
    gene_patterns,
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
    title_size = 16,
    facet_title_size = 14
) {
  library(ggplot2)
  library(dplyr)
  
  obj_tmp <- seurat_obj
  if (!is.null(assay)) DefaultAssay(obj_tmp) <- assay
  
  # --- znajdź wszystkie geny pasujące do patternów ---
  all_features <- rownames(GetAssayData(obj_tmp, assay = assay, slot = slot))
  matched <- unlist(lapply(gene_patterns, function(pat)
    grep(paste0("-", pat, "$"), all_features, value = TRUE)
  ))
  matched <- unique(matched)
  
  if (length(matched) == 0)
    stop("❌ No genes found matching provided patterns.")
  message("✅ Found ", length(matched), " matching genes: ", paste(matched, collapse = ", "))
  
  # --- embeddings ---
  emb <- as.data.frame(Embeddings(obj_tmp, reduction = reduction))
  colnames(emb) <- c("dim1", "dim2")[seq_len(ncol(emb))]
  
  # --- expression matrix ---
  expr_mat <- GetAssayData(obj_tmp, assay = assay, slot = slot)[matched, , drop = FALSE]
  
  # --- binary matrix (1 = detected, 0 = not detected) ---
  detected_mat <- (expr_mat > detection_threshold) * 1
  
  # --- sum of detections per cell ---
  detection_sum <- colSums(detected_mat)
  
  # --- combine ---
  data_df <- emb %>%
    mutate(
      detected_sum = detection_sum,
      detected_any = as.numeric(detection_sum > 0)
    )
  
  # --- metadata ---
  if (!is.null(split.by)) data_df$treatment <- obj_tmp@meta.data[[split.by]]
  if (!is.null(cluster_col) && cluster_col %in% colnames(obj_tmp@meta.data))
    data_df$cluster <- obj_tmp@meta.data[[cluster_col]]
  
  # --- reorder split.by ---
  if (!is.null(split.by)) {
    current_values <- unique(data_df$treatment)
    if (is.null(split_levels)) {
      preferred_order <- c("CTRL", "CORT", "DEX")
      split_levels <- intersect(preferred_order, current_values)
      split_levels <- c(split_levels, setdiff(current_values, split_levels))
    }
    data_df$treatment <- factor(data_df$treatment, levels = split_levels)
  }
  
  # --- compute stats ---
  if (is.null(split.by)) {
    n_any <- sum(data_df$detected_any)
    n_total <- nrow(data_df)
    percent_any <- round(100 * n_any / n_total, 1)
    title_text <- paste0(
      length(matched), " genes detected in ", n_any, "/", n_total,
      " cells (", percent_any, "%)"
    )
  } else {
    summary_df <- as.data.frame(table(data_df$treatment))
    colnames(summary_df) <- c("treatment", "total_cells")
    
    det_summary <- as.data.frame(table(data_df$treatment[data_df$detected_any == 1]))
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
    
    if (!is.null(split_levels)) {
      data_df$treatment <- factor(data_df$treatment, levels = split_levels)
      facet_labels_ordered <- data_df %>%
        distinct(treatment, facet_label) %>%
        arrange(match(treatment, split_levels)) %>%
        pull(facet_label)
      data_df$facet_label <- factor(data_df$facet_label, levels = facet_labels_ordered)
    }
    
    title_text <- paste0("Combined detection of ", length(matched), " genes")
  }
  
  # --- return data if needed ---
  if (return_data) return(data_df)
  
  # --- plotting ---
  p <- ggplot(data_df, aes(x = dim1, y = dim2, color = detected_sum)) +
    geom_point(size = point_size) +
    scale_color_gradientn(
      colors = cols,
      name = paste0("Detected genes (", length(matched), ")")
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = legend_position,
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(size = title_size, color = "black"),
      strip.text = element_text(size = facet_title_size, color = "black"),
      text = element_text(color = "black")
    ) +
    labs(title = title_text)
  
  if (!is.null(split.by)) {
    p <- p + facet_wrap(~facet_label)
  }
  
  return(p)
}


# genes <- c(
#   "ADAMTS1",
#   "ADARB1",
#   "AFF1",
#   "ARHGEF3",
#   "ARRDC2",
#   "BCL2L1",
#   "BCL6",
#   "CDKN1A",
#   "CDO1",
#   "CHKA",
#   "CRISPLD2",
#   "DCUN1D3",
#   "DDIT4",
#   "EGFR",
#   "ERRFI1",
#   "FAM107A",
#   "FKBP5",
#   "FOXO3",
#   "HSD17B11",
#   "IL6R",
#   "KLF13",
#   "KLF15",
#   "KLF9",
#   "LIMS2",
#   "LPIN1",
#   "MAP3K6",
#   "MBNL2",
#   "MERTK",
#   "NEDD9",
#   "NFKBIA",
#   "PDK4",
#   "PER1",
#   "PIK3R1",
#   "RASSF4",
#   "RHOU",
#   "RIN2",
#   "SESN1",
#   "SGK1",
#   "SPSB1",
#   "STIM1",
#   "TSC22D3",
#   "USP53",
#   "ZBTB16",
#   "ZFP36",
#   "ZHX3"
# )
# 
# seurat_featurePlot_multiBinary(
#   seurat_obj = integrated_pilotDexCortCtrl_filtered,
#   gene_patterns = genes,
#   split.by = "treatment",
#   split_levels = c("CTRL", "CORT", "DEX"),
#   cluster_col = "res0.3_RNA_noScaled_brain",
#   cols = c("gray90", "firebrick"),
#   detection_threshold = 0,
#   legend_position = "bottom"
# )

