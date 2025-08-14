plot_umap_facet_by_sample <- function(
    seurat_data,
    plot_title = "UMAP split by sample",
    sample_order = NULL,
    facet_nrow = NULL,
    facet_ncol = NULL,
    label_clusters = TRUE,
    text_size_axis_title = 10,
    text_size_axis_text = 9,
    text_size_legend_title = 10,
    text_size_legend_text = 9,
    text_size_facet = 10,
    text_size_plot_title = 12
) {
  stopifnot("umap" %in% names(seurat_data@reductions))
  
  # 1. UMAP i metadane
  umap_df <- as.data.frame(Embeddings(seurat_data, reduction = "umap")) %>%
    rownames_to_column("cell_id")
  
  metadata_df <- seurat_data@meta.data %>%
    rownames_to_column("cell_id") %>%
    select(cell_id, seurat_clusters, bc1_well, condition)
  
  umap_df <- umap_df %>%
    left_join(metadata_df, by = "cell_id") %>%
    mutate(
      cluster = factor(seurat_clusters),
      sample_label = paste(condition, bc1_well, sep = " - ")
    )
  
  # 2. Kolejność próbek
  if (!is.null(sample_order)) {
    umap_df$sample_label <- factor(umap_df$sample_label, levels = sample_order)
  } else {
    umap_df$sample_label <- factor(umap_df$sample_label)
  }
  
  # 3. Pozycje etykiet klastrów (jeśli potrzebne)
  if (label_clusters) {
    label_df <- umap_df %>%
      group_by(cluster, sample_label) %>%
      summarise(
        x = median(umap_1),
        y = median(umap_2),
        .groups = "drop"
      )
  }
  
  # 4. Wykres
  p <- ggplot(umap_df, aes(x = umap_1, y = umap_2, color = cluster)) +
    geom_point(size = 0.5, alpha = 0.7)
  
  if (label_clusters) {
    p <- p + geom_text(
      data = label_df,
      aes(x = x, y = y, label = cluster),
      inherit.aes = FALSE,
      size = 3,
      color = "black"
    )
  }
  
  p <- p +
    labs(
      title = plot_title,
      color = "Cluster"
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(size = text_size_legend_title, face = "bold"),
      legend.text = element_text(size = text_size_legend_text),
      axis.title = element_text(size = text_size_axis_title),
      axis.text = element_text(size = text_size_axis_text),
      strip.text = element_text(size = text_size_facet),
      plot.title = element_text(size = text_size_plot_title, hjust = 0.5, face = "bold")
    ) +
    guides(
      color = guide_legend(nrow = 1, override.aes = list(size = 3))
    )
  
  # 5. Facetowanie
  facet_formula <- as.formula("~ sample_label")
  if (!is.null(facet_nrow) || !is.null(facet_ncol)) {
    p <- p + facet_wrap(facet_formula, nrow = facet_nrow, ncol = facet_ncol)
  } else {
    p <- p + facet_wrap(facet_formula)
  }
  
  return(p)
}
