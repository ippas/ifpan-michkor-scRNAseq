plot_umap_cluster_faceted <- function(
    seurat_data,
    facet_var = "sample_label",
    cluster_var = "seurat_clusters",
    plot_title = "UMAP faceted plot",
    facet_order = NULL,
    facet_rename = NULL,
    facet_nrow = NULL,
    facet_ncol = NULL,
    label_clusters = TRUE,
    text_size_label = 3,
    text_size_axis_title = 10,
    text_size_axis_text = 9,
    text_size_legend_title = 10,
    text_size_legend_text = 9,
    text_size_facet = 10,
    text_size_plot_title = 12,
    n_row_legend = 1,
    mapper_color = NULL   # ★ NOWY ARGUMENT ★
) {
  stopifnot("umap" %in% names(seurat_data@reductions))
  
  # UMAP coords
  umap_df <- as.data.frame(Seurat::Embeddings(seurat_data, reduction = "umap")) |>
    tibble::rownames_to_column("cell_id")
  
  metadata_df <- seurat_data@meta.data |>
    tibble::rownames_to_column("cell_id")
  
  facet_group_present <- TRUE
  if (is.null(facet_var)) {
    facet_group_present <- FALSE
  } else if (!(facet_var %in% colnames(metadata_df))) {
    stop(paste0("Column '", facet_var, "' not found in metadata"))
  }
  
  if (!(cluster_var %in% colnames(metadata_df))) {
    stop(paste0("Column '", cluster_var, "' not found in metadata"))
  }
  
  # Merge data
  umap_df <- umap_df |>
    dplyr::left_join(metadata_df, by = "cell_id") |>
    dplyr::mutate(cluster = factor(!!rlang::sym(cluster_var)))
  
  # Facet group
  if (facet_group_present) {
    umap_df <- umap_df |>
      dplyr::mutate(facet_group = !!rlang::sym(facet_var))
    
    # Order facet groups
    if (!is.null(facet_order)) {
      umap_df$facet_group <- factor(umap_df$facet_group, levels = facet_order)
    } else {
      umap_df$facet_group <- factor(umap_df$facet_group)
    }
    
    # Rename facet groups
    if (!is.null(facet_rename)) {
      levels(umap_df$facet_group) <- sapply(
        levels(umap_df$facet_group),
        function(lv) if (lv %in% names(facet_rename)) facet_rename[[lv]] else lv
      )
    }
  }
  
  # Cluster label positions
  if (label_clusters) {
    if (facet_group_present) {
      label_df <- umap_df |>
        dplyr::group_by(cluster, facet_group) |>
        dplyr::summarise(
          x = median(umap_1),
          y = median(umap_2),
          .groups = "drop"
        )
    } else {
      label_df <- umap_df |>
        dplyr::group_by(cluster) |>
        dplyr::summarise(
          x = median(umap_1),
          y = median(umap_2),
          .groups = "drop"
        )
    }
  }
  
  # Base plot
  p <- ggplot2::ggplot(umap_df, ggplot2::aes(x = umap_1, y = umap_2, color = cluster)) +
    ggplot2::geom_point(size = 0.5, alpha = 0.7)
  
  # ★★ APPLY MAPPER_COLOR IF GIVEN ★★
  if (!is.null(mapper_color)) {
    
    # Ensure mapper_color is named
    if (is.null(names(mapper_color))) {
      stop("mapper_color must be a named vector: names = clusters, values = colors")
    }
    
    clusters_present <- levels(umap_df$cluster)
    clusters_in_mapper <- names(mapper_color)
    
    # Warn if missing
    missing_clusters <- setdiff(clusters_present, clusters_in_mapper)
    if (length(missing_clusters) > 0) {
      warning("mapper_color is missing colors for clusters: ",
              paste(missing_clusters, collapse = ", "))
    }
    
    p <- p + ggplot2::scale_color_manual(
      values = mapper_color,
      drop = FALSE
    )
  }
  
  # Cluster labels
  if (label_clusters) {
    p <- p + ggplot2::geom_text(
      data = label_df,
      ggplot2::aes(x = x, y = y, label = cluster),
      inherit.aes = FALSE,
      size = text_size_label,
      color = "black"
    )
  }
  
  # Themes & labels
  p <- p +
    ggplot2::labs(
      title = plot_title,
      color = "Cluster"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = ggplot2::element_text(size = text_size_legend_title, face = "bold", color = "black"),
      legend.text = ggplot2::element_text(size = text_size_legend_text, color = "black"),
      axis.title = ggplot2::element_text(size = text_size_axis_title, color = "black"),
      axis.text = ggplot2::element_text(size = text_size_axis_text, color = "black"),
      strip.text = ggplot2::element_text(size = text_size_facet, color = "black"),
      plot.title = ggplot2::element_text(size = text_size_plot_title, hjust = 0.5, face = "bold", color = "black")
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(nrow = n_row_legend, override.aes = list(size = 3))
    )
  
  # Add faceting
  if (facet_group_present) {
    facet_formula <- stats::as.formula("~ facet_group")
    if (!is.null(facet_nrow) || !is.null(facet_ncol)) {
      p <- p + ggplot2::facet_wrap(facet_formula, nrow = facet_nrow, ncol = facet_ncol)
    } else {
      p <- p + ggplot2::facet_wrap(facet_formula)
    }
  }
  
  return(p)
}
