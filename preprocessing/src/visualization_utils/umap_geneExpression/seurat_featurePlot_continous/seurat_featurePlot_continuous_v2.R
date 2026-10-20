seurat_featurePlot_continuous_v2 <- function(
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
    legend_position = "bottom",
    log_transform = FALSE,
    return_data = FALSE,
    winsorize_top = 0,        # 🔥 NEW argument
    title_size = 16,
    facet_title_size = 14
) {
  library(ggplot2)
  library(dplyr)
  library(scales)
  
  obj_tmp <- seurat_obj
  if (!is.null(assay)) DefaultAssay(obj_tmp) <- assay
  
  all_features <- rownames(GetAssayData(obj_tmp, assay = assay, slot = slot))
  matched <- grep(paste0("-", gene_pattern, "$"), all_features, value = TRUE)
  
  if (length(matched) == 0) stop(paste0("❌ No gene found matching: ", gene_pattern))
  if (length(matched) > 1) {
    message("⚠ Multiple matches found, using first: ", paste(matched, collapse = ", "))
    matched <- matched[1]
  } else {
    message("✅ Found feature: ", matched)
  }
  
  # UMAP embedding
  emb <- as.data.frame(Embeddings(obj_tmp, reduction = reduction))
  colnames(emb) <- c("dim1", "dim2")[seq_len(ncol(emb))]
  
  expr_values <- GetAssayData(obj_tmp, assay = assay, slot = slot)[matched, ]
  if (log_transform) expr_values <- log1p(expr_values)
  
  # --- 🔥 NEW WINSORIZATION ---
  if (winsorize_top > 0 && winsorize_top < 1) {
    cutoff <- quantile(expr_values, probs = 1 - winsorize_top)
    expr_values <- pmin(expr_values, cutoff)
    message("✂ Winsorized top ", winsorize_top * 100, "% at cutoff=", round(cutoff, 3))
  }
  
  data_df <- emb %>% mutate(expression = expr_values)
  
  if (!is.null(split.by)) data_df$treatment <- obj_tmp@meta.data[[split.by]]
  if (!is.null(cluster_col)) data_df$cluster <- obj_tmp@meta.data[[cluster_col]]
  
  if (!is.null(split.by)) {
    current <- unique(data_df$treatment)
    if (is.null(split_levels)) {
      preferred <- c("CTRL", "CORT", "DEX")
      split_levels <- c(intersect(preferred, current), setdiff(current, preferred))
    }
    data_df$treatment <- factor(data_df$treatment, levels = split_levels)
  }
  
  if (return_data) return(data_df)
  
  p <- ggplot(data_df, aes(x = dim1, y = dim2, color = expression)) +
    geom_point(size = point_size) +
    scale_color_gradient(low = cols[1], high = cols[2]) +
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
    labs(title = paste0(matched, " (", assay, ")"), color = "Expression")
  
  if (!is.null(split.by)) p <- p + facet_wrap(~ treatment)
  
  return(p)
}
