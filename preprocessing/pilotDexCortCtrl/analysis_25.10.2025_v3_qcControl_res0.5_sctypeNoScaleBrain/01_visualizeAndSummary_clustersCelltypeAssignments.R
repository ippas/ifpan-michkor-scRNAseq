# ##############################################################################
# ---- uses data ----
# ##############################################################################
# gdzie zapisać:
outdir <- "results/results_analysis_25.10.2025_v3_qcControl_res0.5_sctypeNoScaleBrain/UMAP/"

dev.off()

# ##############################################################################
svg(
  filename = file.path(outdir, "umap_res0.5_clusters20.svg"),
  width = 8,  # cale
  height = 11,  # cale
  pointsize = 12,
  bg = "white"
)

# rysunek
plot_umap_cluster_faceted(
  integrated_pilotDexCortCtrl_filtered,
  cluster_var = "seurat_clusters_res0.3",
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
  n_row_legend = 3
)

# zamknij urządzenie
dev.off()

# ##############################################################################
svg(
  filename = file.path(outdir, "umap_res0.5_clusters20_wrapTreatment.svg"),
  width = 8*3,  # cale
  height = 11,  # cale
  pointsize = 12,
  bg = "white"
)

# rysunek
plot_umap_cluster_faceted(
  integrated_pilotDexCortCtrl_filtered,
  cluster_var = "seurat_clusters_res0.5",
  plot_title = "UMAP Projection Across All Samples",
  facet_var = "treatment",
  label_clusters = TRUE,
  text_size_axis_title = 18,
  text_size_axis_text = 18,
  text_size_legend_title = 18,
  text_size_legend_text = 18,
  text_size_facet = 18,
  text_size_plot_title = 18,
  text_size_label = 7,
  n_row_legend = 3
)

# zamknij urządzenie
dev.off()

# ##############################################################################
svg(
  filename = file.path(outdir, "umap_res0.5_cellType15.svg"),
  width = 8,  # cale
  height = 12,  # cale
  pointsize = 12,
  bg = "white"
)

# rysunek
plot_umap_cluster_faceted(
  integrated_pilotDexCortCtrl_filtered,
  cluster_var = "res0.5_integrated_noScaled_brain",
  plot_title = "UMAP Projection Across All Samples",
  facet_var = NULL,
  label_clusters = TRUE,
  text_size_axis_title = 18,
  text_size_axis_text = 18,
  text_size_legend_title = 18,
  text_size_legend_text = 18,
  text_size_facet = 18,
  text_size_plot_title = 18,
  text_size_label = 6,
  n_row_legend = 7
)

# zamknij urządzenie
dev.off()


# ##############################################################################
svg(
  filename = file.path(outdir, "umap_res0.5_cellType15_wrapTreatment.svg"),
  width = 8*3,  # cale
  height = 11,  # cale
  pointsize = 12,
  bg = "white"
)

# rysunek
plot_umap_cluster_faceted(
  integrated_pilotDexCortCtrl_filtered,
  cluster_var = "res0.3_RNA_noScaled_brain",
  plot_title = "UMAP Projection Across All Samples",
  facet_var = "treatment",
  label_clusters = TRUE,
  text_size_axis_title = 18,
  text_size_axis_text = 18,
  text_size_legend_title = 18,
  text_size_legend_text = 18,
  text_size_facet = 18,
  text_size_plot_title = 18,
  text_size_label = 7,
  n_row_legend = 3
)

# zamknij urządzenie
dev.off()



# ##############################################################################
# ---- summary number of cellstypes ----
# ##############################################################################
integrated_pilotDexCortCtrl_filtered@meta.data$res0.5_integrated_noScaled_brain %>% 
  table %>% 
  as.data.frame() %>% 
  set_colnames(c("cett type", "n cells"))


integrated_pilotDexCortCtrl_filtered@meta.data %>% 
  select(c(res0.5_integrated_noScaled_brain, treatment)) %>% 
  rownames_to_column(var = "cell") %>% 
  dplyr::rename(cellType = res0.5_integrated_noScaled_brain) %>% 
  group_by(cellType, treatment) %>% 
  nest() %>% 
  mutate(n_cells = map(data, ~ .x %>% nrow)) %>% 
  unnest(n_cells) %>% 
  select(-data) %>% 
  tidyr::pivot_wider(
    names_from = treatment,
    values_from = n_cells,
    values_fill = 0             # 👈 jeśli niektóre kombinacje nie występują
  ) %>% 
  mutate(sum_cells = CTRL + DEX + CORT) %>% 
  arrange(desc(sum_cells))



integrated_pilotDexCortCtrl_filtered@meta.data %>% 
  select(c(res0.5_integrated_noScaled_brain, seurat_clusters_res0.5)) %>% 
  rownames_to_column() %>% 
  select(-rowname) %>% 
  unique() %>% 
  group_by(res0.5_integrated_noScaled_brain) %>% nest %>% 
  mutate(n_cluster = map(data, ~ .x %>% nrow)) %>% 
  unnest(n_cluster) %>% arrange(n_cluster)




# ##############################################################################
plot_umap_cluster_faceted(
  integrated_pilotDexCortCtrl_filtered,
  cluster_var = "treatment",
  plot_title = "",
  facet_var = "sample",
  label_clusters = TRUE,
  text_size_axis_title = 18,
  text_size_axis_text = 18,
  text_size_legend_title = 18,
  text_size_legend_text = 18,
  text_size_facet = 18,
  text_size_plot_title = 18,
  text_size_label = 7,
  n_row_legend = 4
)
  



svg(
  filename = file.path(outdir, "umap_res0.3_cellType11.svg"),
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

svg(
  filename = file.path(outdir, "umap_res0.3_clusters16.svg"),
  width = 8,  # cale
  height = 12,  # cale
  pointsize = 12,
  bg = "white"
)

plot_umap_cluster_faceted(
  integrated_pilotDexCortCtrl_filtered,
  cluster_var = "res0.3_RNA_noScaled_allMarkers",
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
  n_row_legend = 3
)

dev.off()



plot_umap_cluster_faceted(
  integrated_pilotDexCortCtrl_filtered,
  cluster_var = "seurat_clusters_res0.3",
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
  n_row_legend = 3
)




svg(
  filename = file.path(outdir, "umap_res0.3_RNA_clusters16_allMarkers.svg"),
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
  n_row_legend = 3
)


dev.off()
