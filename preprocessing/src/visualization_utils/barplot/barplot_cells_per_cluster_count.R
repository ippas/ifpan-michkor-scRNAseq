barplot_cells_per_cluster_count <- function(
    seurat_obj,
    resolution = 0.5,
    facet_var = "treatment",          # set to NULL for a single combined plot
    facet_order = c("CTRL", "DEX", "CORT"),
    cluster_colors = NULL,
    # y-axis control
    fix_y_axis = FALSE,
    y_min = 0,
    y_max = NULL,
    # text sizes (points)
    axis_text_pt    = 12,
    axis_title_pt   = 14,
    plot_title_pt   = 16,
    legend_title_pt = 14,
    legend_text_pt  = 12,
    bar_label_pt    = 12,
    sigma_label_pt  = 18,
    # labels
    plot_title = "Cell counts per cluster",
    x_label    = NULL,
    y_label    = "n cells"
) {
  # --- Basic checks ---
  stopifnot(inherits(seurat_obj, "Seurat"))
  meta <- seurat_obj@meta.data
  
  # Cluster column from resolution
  cluster_col <- paste0("seurat_clusters_res", resolution)
  if (!cluster_col %in% colnames(meta)) {
    stop("Missing column '", cluster_col, "' in Seurat meta.data.")
  }
  
  cluster_sym <- rlang::sym(cluster_col)
  pt_to_geom  <- function(pt) pt / 2.845
  if (is.null(x_label)) x_label <- paste0("Cluster (res", resolution, ")")
  
  # Color builder helper
  build_colors <- function(levels, cluster_colors) {
    if (is.null(cluster_colors)) {
      cols <- scales::hue_pal()(length(levels))
      names(cols) <- levels
      return(cols)
    }
    if (is.null(names(cluster_colors))) {
      if (length(cluster_colors) < length(levels)) {
        stop("Provided 'cluster_colors' has fewer entries than clusters.")
      }
      names(cluster_colors) <- levels
      return(cluster_colors)
    }
    missing_levels <- setdiff(levels, names(cluster_colors))
    if (length(missing_levels) > 0) {
      stop("Missing colors for cluster levels: ", paste(missing_levels, collapse = ", "))
    }
    cluster_colors[levels]
  }
  
  # ----------------------------------
  # Case 1: Single plot (facet_var = NULL)
  # ----------------------------------
  if (is.null(facet_var)) {
    counts_df <- meta %>%
      dplyr::count(!!cluster_sym, name = "n_cells")
    stats_df <- counts_df %>%
      dplyr::summarise(
        total_cells = sum(n_cells),
        max_cells   = max(n_cells),
        n_clusters  = dplyr::n_distinct(!!cluster_sym),
        .groups     = "drop"
      )
    cluster_levels <- sort(unique(dplyr::pull(counts_df, !!cluster_sym)))
    cluster_colors <- build_colors(cluster_levels, cluster_colors)
    
    p <- ggplot2::ggplot(
      counts_df,
      ggplot2::aes(x = !!cluster_sym, y = n_cells, fill = !!cluster_sym)
    ) +
      ggplot2::geom_col() +
      ggplot2::geom_text(
        ggplot2::aes(label = n_cells),
        angle = 45, hjust = -0.1, vjust = -0.5,
        size = pt_to_geom(bar_label_pt)
      ) +
      ggplot2::annotate(
        "text",
        x = (stats_df$n_clusters + 1) / 2,
        y = stats_df$max_cells * 1.06,
        label = paste0("Sigma[cells] == ", stats_df$total_cells),
        parse = TRUE,
        fontface = "bold",
        size = pt_to_geom(sigma_label_pt)
      ) +
      ggplot2::scale_fill_manual(values = cluster_colors) +
      ggplot2::coord_cartesian(
        clip = "off",
        ylim = if (fix_y_axis) c(y_min, y_max %||% (max(counts_df$n_cells) * 1.12)) else NULL
      ) +
      ggplot2::labs(
        x = x_label, y = y_label, fill = "Cluster", title = plot_title
      ) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        axis.text.x  = ggplot2::element_text(size = axis_text_pt),
        axis.text.y  = ggplot2::element_text(size = axis_text_pt),
        axis.title.x = ggplot2::element_text(size = axis_title_pt),
        axis.title.y = ggplot2::element_text(size = axis_title_pt),
        plot.title   = ggplot2::element_text(size = plot_title_pt),
        legend.title = ggplot2::element_text(size = legend_title_pt),
        legend.text  = ggplot2::element_text(size = legend_text_pt)
      )
    return(p)
  }
  
  # ----------------------------------
  # Case 2: Faceted plot
  # ----------------------------------
  if (!facet_var %in% colnames(meta)) {
    stop("Missing facet variable '", facet_var, "' in Seurat meta.data.")
  }
  facet_sym <- rlang::sym(facet_var)
  
  counts_df <- meta %>%
    dplyr::count(!!cluster_sym, !!facet_sym, name = "n_cells") %>%
    dplyr::mutate(!!facet_sym := factor(as.character(!!facet_sym), levels = facet_order))
  
  stats_df <- counts_df %>%
    dplyr::group_by(!!facet_sym) %>%
    dplyr::summarise(
      total_cells = sum(n_cells),
      max_cells   = max(n_cells),
      n_clusters  = dplyr::n_distinct(!!cluster_sym),
      .groups     = "drop"
    ) %>%
    dplyr::mutate(!!facet_sym := factor(as.character(!!facet_sym), levels = facet_order))
  
  cluster_levels <- sort(unique(dplyr::pull(counts_df, !!cluster_sym)))
  cluster_colors <- build_colors(cluster_levels, cluster_colors)
  
  p <- ggplot2::ggplot(
    counts_df,
    ggplot2::aes(x = !!cluster_sym, y = n_cells, fill = !!cluster_sym)
  ) +
    ggplot2::geom_col() +
    ggplot2::geom_text(
      ggplot2::aes(label = n_cells),
      angle = 45, hjust = -0.1, vjust = -0.5,
      size = pt_to_geom(bar_label_pt)
    ) +
    ggplot2::geom_text(
      data = stats_df,
      ggplot2::aes(
        x = (n_clusters + 1) / 2,
        y = max_cells * 1.06,
        label = paste0("Sigma[cells] == ", total_cells)
      ),
      inherit.aes = FALSE,
      parse = TRUE,
      fontface = "bold",
      size = pt_to_geom(sigma_label_pt)
    ) +
    ggplot2::scale_fill_manual(values = cluster_colors) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(!!facet_sym),
      scales = if (fix_y_axis) "fixed" else "free_y"
    ) +
    ggplot2::coord_cartesian(
      clip = "off",
      ylim = if (fix_y_axis) c(y_min, y_max %||% (max(counts_df$n_cells) * 1.12)) else NULL
    ) +
    ggplot2::labs(
      x = x_label, y = y_label, fill = "Cluster", title = plot_title
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(size = axis_text_pt),
      axis.text.y  = ggplot2::element_text(size = axis_text_pt),
      axis.title.x = ggplot2::element_text(size = axis_title_pt),
      axis.title.y = ggplot2::element_text(size = axis_title_pt),
      plot.title   = ggplot2::element_text(size = plot_title_pt),
      legend.title = ggplot2::element_text(size = legend_title_pt),
      legend.text  = ggplot2::element_text(size = legend_text_pt)
    )
  
  return(p)
}
