svg(
  filename = file.path("results/figures/umapCellType_07.12.2025/", "umap_res0.3_cellType11_allColors_07.12.2025.svg"),
  width = 8,  # cale
  height = 12,  # cale
  pointsize = 12,
  bg = "white"
)

plot_umap_cluster_faceted(
  integrated_pilotDexCortCtrl_filtered,
  cluster_var = "res0.3_RNA_noScaled_brain",
  plot_title = "UMAP Projection Across All Samples",
  facet_var = NULL,
  label_clusters = TRUE,
  text_size_axis_title = 18,
  text_size_axis_text = 18,
  text_size_legend_title = 18,
  text_size_legend_text = 18,
  text_size_facet = 18,
  text_size_plot_title = 18,
  text_size_label = 7,
  n_row_legend = 6
)

dev.off()

mapper_color <- c(
  "Radial glial cells"              = "#b385ffff",  # blue
  "Neural stem cells"               = "gray90",  # light blue
  "Myelinating Schwann cells"       = "gray90",  # dark green
  "GABAergic neurons"               = "gray90",  # orange
  "Oligodendrocyte precursor cells" = "#00badeff",  # light purple
  "Immature neurons"                = "#aea200ff",  # light orange
  "Astrocytes"                      = "#f6766dff",  # purple
  "Mature neurons"                  = "#64b100ff",  # red
  "Oligodendrocytes"                = "#01a6ffff",  # deep purple
  "Tanycytes"                       = "gray90",  # green
  "Schwann precursor cells"         = "gray90"   # light green
)


svg(
  filename = file.path("results/figures/umapCellType_07.12.2025/", "umap_res0.3_cellType11_6cellTypesColors_07.12.2025.svg"),
  width = 8,  # cale
  height = 12,  # cale
  pointsize = 12,
  bg = "white"
)

plot_umap_cluster_faceted(
  integrated_pilotDexCortCtrl_filtered,
  cluster_var = "res0.3_RNA_noScaled_brain",
  plot_title = "UMAP Projection Across All Samples",
  facet_var = NULL,
  label_clusters = TRUE,
  text_size_axis_title = 18,
  text_size_axis_text = 18,
  text_size_legend_title = 18,
  text_size_legend_text = 18,
  text_size_facet = 18,
  text_size_plot_title = 18,
  text_size_label = 7,
  n_row_legend = 6,
  mapper_color = mapper_color
)

dev.off()
