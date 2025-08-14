
# ##############################################################################
# ---- UMAP visualization ----
# ##############################################################################

plot_umap_cluster_faceted(
  merged_pilotDexCortCtrl,
  cluster_var = "seurat_clusters_res0.5",
  plot_title = "UMAP Projection Across All Samplest",
  facet_var = NULL,
  # facet_order = c("CTRL", "DEX", "CORT"),
  label_clusters = TRUE,               # czy pokazywać etykiety klastrów
  text_size_axis_title = 18,           # rozmiar tytułów osi
  text_size_axis_text = 14,             # rozmiar tekstu na osiach (tick labels)
  text_size_legend_title = 16,         # rozmiar tytułu legendy
  text_size_legend_text = 16,          # rozmiar tekstu legendy
  text_size_facet = 14,                # rozmiar tekstu podtytułów paneli (facetów)
  text_size_plot_title = 18 ,           # rozmiar tytułu wykresu
  text_size_label = 6,
  n_row_legend = 2
)

plot_umap_cluster_faceted(
  merged_pilotDexCortCtrl,
  facet_var = "sublibrary",
  cluster_var = "seurat_clusters_res0.5",
  plot_title = "UMAP Projection by Sublibrary",
  # facet_order = c("CTRL", "DEX", "CORT"),
  label_clusters = TRUE,               # czy pokazywać etykiety klastrów
  text_size_axis_title = 14,           # rozmiar tytułów osi
  text_size_axis_text = 10,             # rozmiar tekstu na osiach (tick labels)
  text_size_legend_title = 12,         # rozmiar tytułu legendy
  text_size_legend_text = 12,          # rozmiar tekstu legendy
  text_size_facet = 12,                # rozmiar tekstu podtytułów paneli (facetów)
  text_size_plot_title = 14 ,           # rozmiar tytułu wykresu
  text_size_label = 4,
  n_row_legend = 2
)

plot_umap_cluster_faceted(
  merged_pilotDexCortCtrl,
  facet_var = "treatment",
  cluster_var = "seurat_clusters_res0.022",
  plot_title = "UMAP Projection by Treatment, res=0.022",
  facet_order = c("CTRL", "DEX", "CORT"),
  label_clusters = TRUE,               # czy pokazywać etykiety klastrów
  text_size_axis_title = 14,           # rozmiar tytułów osi
  text_size_axis_text = 10,             # rozmiar tekstu na osiach (tick labels)
  text_size_legend_title = 12,         # rozmiar tytułu legendy
  text_size_legend_text = 12,          # rozmiar tekstu legendy
  text_size_facet = 12,                # rozmiar tekstu podtytułów paneli (facetów)
  text_size_plot_title = 14 ,           # rozmiar tytułu wykresu
  text_size_label = 4
)

plot_umap_cluster_faceted(
  merged_pilotDexCortCtrl,
  facet_var = "bc1_well",
  cluster_var = "seurat_clusters_res0.5",
  facet_order = paste0("A", 1:12),
  facet_rename = c(
    A1 = "A1 - Ctrl",
    A2 = "A2 - Ctrl",
    A3 = "A3 - Ctrl",
    A4 = "A4 - Ctrl",
    A5 = "A5 - Dex",
    A6 = "A6 - Dex",
    A7 = "A7 - Dex",
    A8 = "A8 - Dex",
    A9 = "A9 - Cort",
    A10 = "A10 - Cort",
    A11 = "A11 - Cort",
    A12 = "A12 - Cort"
  ),
  plot_title = "UMAP, cluster per sample",
  facet_nrow = 3,
  label_clusters = TRUE,               # czy pokazywać etykiety klastrów
  text_size_axis_title = 14,           # rozmiar tytułów osi
  text_size_axis_text = 10,             # rozmiar tekstu na osiach (tick labels)
  text_size_legend_title = 12,         # rozmiar tytułu legendy
  text_size_legend_text = 12,          # rozmiar tekstu legendy
  text_size_facet = 12,                # rozmiar tekstu podtytułów paneli (facetów)
  text_size_plot_title = 14 ,           # rozmiar tytułu wykresu
  text_size_label = 4
)


plot_umap_cluster_faceted <- 
  function(
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
    n_row_legend = 1
  ) {
    stopifnot("umap" %in% names(seurat_data@reductions))
    
    umap_df <- base::as.data.frame(Seurat::Embeddings(seurat_data, reduction = "umap")) |>
      tibble::rownames_to_column("cell_id")
    
    metadata_df <- seurat_data@meta.data |>
      tibble::rownames_to_column("cell_id")
    
    facet_group_present <- TRUE
    if (base::is.null(facet_var)) {
      facet_group_present <- FALSE
    } else if (!(facet_var %in% base::colnames(metadata_df))) {
      base::stop(base::paste0("Column '", facet_var, "' not found in metadata"))
    }
    
    if (!(cluster_var %in% base::colnames(metadata_df))) {
      base::stop(base::paste0("Column '", cluster_var, "' not found in metadata"))
    }
    
    if (facet_group_present) {
      umap_df <- umap_df |>
        dplyr::left_join(metadata_df, by = "cell_id") |>
        dplyr::mutate(
          cluster = base::factor(rlang::.data[[cluster_var]]),
          facet_group = rlang::.data[[facet_var]]
        )
      
      if (!base::is.null(facet_order)) {
        umap_df[["facet_group"]] <- base::factor(umap_df[["facet_group"]], levels = facet_order)
      } else {
        umap_df[["facet_group"]] <- base::factor(umap_df[["facet_group"]])
      }
      
      if (!base::is.null(facet_rename)) {
        base::levels(umap_df[["facet_group"]]) <- base::sapply(
          base::levels(umap_df[["facet_group"]]),
          function(lv) if (lv %in% base::names(facet_rename)) facet_rename[[lv]] else lv
        )
      }
      
      if (label_clusters) {
        label_df <- umap_df |>
          dplyr::group_by(cluster, facet_group) |>
          dplyr::summarise(
            x = stats::median(umap_1),
            y = stats::median(umap_2),
            .groups = "drop"
          )
      }
      
    } else {
      umap_df <- umap_df |>
        dplyr::left_join(metadata_df, by = "cell_id") |>
        dplyr::mutate(cluster = base::factor(rlang::.data[[cluster_var]]))
      
      if (label_clusters) {
        label_df <- umap_df |>
          dplyr::group_by(cluster) |>
          dplyr::summarise(
            x = stats::median(umap_1),
            y = stats::median(umap_2),
            .groups = "drop"
          )
      }
    }
    
    p <- ggplot2::ggplot(umap_df, ggplot2::aes(x = umap_1, y = umap_2, color = cluster)) +
      ggplot2::geom_point(size = 0.5, alpha = 0.7)
    
    if (label_clusters) {
      p <- p + ggplot2::geom_text(
        data = label_df,
        ggplot2::aes(x = x, y = y, label = cluster),
        inherit.aes = FALSE,
        size = text_size_label,
        color = "black"
      )
    }
    
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
        axis.title  = ggplot2::element_text(size = text_size_axis_title, color = "black"),
        axis.text   = ggplot2::element_text(size = text_size_axis_text, color = "black"),
        strip.text  = ggplot2::element_text(size = text_size_facet, color = "black"),
        plot.title  = ggplot2::element_text(size = text_size_plot_title, hjust = 0.5, face = "bold", color = "black")
      ) +
      ggplot2::guides(
        color = ggplot2::guide_legend(nrow = n_row_legend, override.aes = list(size = 3))
      )
    
    if (facet_group_present) {
      facet_formula <- stats::as.formula("~ facet_group")
      if (!base::is.null(facet_nrow) || !base::is.null(facet_ncol)) {
        p <- p + ggplot2::facet_wrap(facet_formula, nrow = facet_nrow, ncol = facet_ncol)
      } else {
        p <- p + ggplot2::facet_wrap(facet_formula)
      }
    }
    
    return(p)
  }

