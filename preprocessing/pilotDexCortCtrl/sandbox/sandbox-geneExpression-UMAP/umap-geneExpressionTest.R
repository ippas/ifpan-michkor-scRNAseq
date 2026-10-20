library(ggplot2)
library(dplyr)

#' Custom ggplot-based FeaturePlot for Seurat
#'
#' Extracts UMAP coordinates, expression values, metadata, and cluster information
#' to generate a publication-quality UMAP colored by gene expression.
#'
#' @param seurat_obj A Seurat object.
#' @param gene_id Full row name of gene (e.g. "ENSG00000109906-ZBTB16").
#' @param assay Which assay to use (default = "RNA").
#' @param slot Which slot to use from assay (default = "counts").
#' @param split.by Metadata column name for faceting (e.g. "treatment").
#' @param split_levels Optional vector to set order of panels (default = CTRL → CORT → DEX).
#' @param cluster_col Metadata column name used as cluster ID (optional).
#' @param reduction Dimensionality reduction to use (default = "umap").
#' @param cols Two-color gradient (default gray90 → firebrick).
#' @param point_size Size of UMAP points.
#' @param legend_position Where to place legend ("bottom", "right", etc.).
#' @param log_transform Logical; whether to apply log1p() to expression.
#' @return A ggplot object.
#' @examples
#' seurat_featurePlot_custom(seurat_obj, "ENSG00000109906-ZBTB16", split.by="treatment")
seurat_featurePlot_custom <- function(
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
    return_data = FALSE
) {
  # create local copy
  obj_tmp <- seurat_obj
  
  # set assay if needed
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
  
  # extract UMAP
  emb <- as.data.frame(Embeddings(obj_tmp, reduction = reduction))
  colnames(emb) <- c("dim1", "dim2")[seq_len(ncol(emb))]
  
  # extract expression values
  expr_values <- GetAssayData(obj_tmp, assay = assay, slot = slot)[matched, ]
  if (log_transform) expr_values <- log1p(expr_values)
  
  # combine info
  data_df <- emb %>%
    mutate(expression = expr_values)
  
  # add metadata
  if (!is.null(split.by)) {
    data_df$treatment <- obj_tmp@meta.data[[split.by]]
  }
  if (!is.null(cluster_col) && cluster_col %in% colnames(obj_tmp@meta.data)) {
    data_df$cluster <- obj_tmp@meta.data[[cluster_col]]
  }
  
  # reorder split.by if applicable
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
    theme_minimal(base_size = 12) +
    theme(
      legend.position = legend_position,
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank()
    ) +
    theme_classic() +
    labs(
      title = paste0(matched, " (", assay, " assay)"),
      color = "Expression"
    )
  
  if (!is.null(split.by)) {
    p <- p + facet_wrap(~ treatment)
  }
  
  return(p)
}

seurat_featurePlot_custom(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  gene_pattern = "FKBP5",
  assay = "RNA",
  slot = "counts",
  split.by = "treatment",
  split_levels = c("CTRL", "CORT", "DEX"),
  cluster_col = "res0.3_RNA_noScaled_brain",
  cols = c("gray90", "firebrick"),
  point_size = 0.4,
  log_transform = FALSE
)

