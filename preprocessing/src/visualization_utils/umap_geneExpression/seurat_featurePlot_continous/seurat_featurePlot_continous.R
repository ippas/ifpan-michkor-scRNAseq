seurat_featurePlot_continuous <- function(
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
    # 🔹 nowe argumenty
    title_size = 16,
    facet_title_size = 14
) {
  library(ggplot2)
  library(dplyr)
  
  # create local copy
  obj_tmp <- seurat_obj
  if (!is.null(assay)) DefaultAssay(obj_tmp) <- assay
  
  # find matching genes by pattern (e.g. "-ZBTB16$")
  all_features <- rownames(GetAssayData(obj_tmp, assay = assay, slot = slot))
  matched <- grep(paste0("-", gene_pattern, "$"), all_features, value = TRUE)
  
  if (length(matched) == 0) {
    stop(paste0("❌ No gene found matching pattern: ", gene_pattern))
  }
  if (length(matched) > 1) {
    message("⚠️ Multiple matches found, using first: ", paste(matched, collapse = ", "))
    matched <- matched[1]
  } else {
    message("✅ Found feature: ", matched)
  }
  
  # extract UMAP embeddings
  emb <- as.data.frame(Embeddings(obj_tmp, reduction = reduction))
  colnames(emb) <- c("dim1", "dim2")[seq_len(ncol(emb))]
  
  # extract expression values
  expr_values <- GetAssayData(obj_tmp, assay = assay, slot = slot)[matched, ]
  if (log_transform) expr_values <- log1p(expr_values)
  
  # combine into dataframe
  data_df <- emb %>% mutate(expression = expr_values)
  
  # add metadata
  if (!is.null(split.by)) data_df$treatment <- obj_tmp@meta.data[[split.by]]
  if (!is.null(cluster_col) && cluster_col %in% colnames(obj_tmp@meta.data))
    data_df$cluster <- obj_tmp@meta.data[[cluster_col]]
  
  # reorder split levels if applicable
  if (!is.null(split.by)) {
    current_values <- unique(data_df$treatment)
    if (is.null(split_levels)) {
      preferred_order <- c("CTRL", "CORT", "DEX")
      split_levels <- intersect(preferred_order, current_values)
      split_levels <- c(split_levels, setdiff(current_values, split_levels))
    }
    data_df$treatment <- factor(data_df$treatment, levels = split_levels)
  }
  
  # return data if requested
  if (return_data) return(data_df)
  
  # ggplot visualization
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
      text = element_text(color = "black")  # 🔹 cały tekst czarny
    ) +
    labs(
      title = paste0(matched, " (", assay, " assay)"),
      color = "Expression"
    )
  
  # facet split
  if (!is.null(split.by)) {
    p <- p + facet_wrap(~ treatment)
  }
  
  return(p)
}


# seurat_featurePlot_continuous(
#   seurat_obj = integrated_pilotDexCortCtrl_filtered,
#   gene_pattern = "ZBTB16",
#   split.by = "treatment",
#   split_levels = c("CTRL", "CORT", "DEX"),
#   cols = c("gray90", "firebrick"),
#   log_transform = TRUE,
#   title_size = 14,
#   facet_title_size = 14
# )

