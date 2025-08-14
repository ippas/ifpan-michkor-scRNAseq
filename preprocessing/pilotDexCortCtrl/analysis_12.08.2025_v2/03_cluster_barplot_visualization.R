barplot_cells_per_cluster_count(
  seurat_obj = integrated_pilotDexCortCtrl,
  resolution = 0.04
)


barplot_cells_per_cluster_fraction(
  seurat_obj = integrated_pilotDexCortCtrl,
  resolution = 0.04,
  facet_var  = "treatment",
  fix_y_axis = TRUE,
  y_min = 0
)


barplot_cells_per_cluster_count(
  seurat_obj = integrated_pilotDexCortCtrl,
  resolution = 0.13
)


barplot_cells_per_cluster_fraction(
  seurat_obj = integrated_pilotDexCortCtrl,
  resolution = 0.13,
  facet_var  = "treatment",
  fix_y_axis = TRUE,
  y_min = 0
)


barplot_cells_per_cluster_count(
  seurat_obj = integrated_pilotDexCortCtrl,
  resolution = 0.2
)


barplot_cells_per_cluster_fraction(
  seurat_obj = integrated_pilotDexCortCtrl,
  resolution = 0.2,
  facet_var  = "treatment",
  fix_y_axis = TRUE,
  y_min = 0
)
