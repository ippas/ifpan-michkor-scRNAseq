# ##############################################################################
# ---- preapare data ----
# ##############################################################################
integrated_pilotDexCortCtrl

read.delim(file = "data/pilotDexCortCtrl/r_outputs/celltypes_trailmaker/celltypes_trailmaker_ScType.csv", 
           sep = ",") %>% 
  mutate(cell_type = str_replace_all(cell_type, " ", "_")) %>% 
  set_colnames(c("id", "celltype_trailmaker")) -> celltypes_trailmaker


new_meta <- integrated_pilotDexCortCtrl@meta.data %>% 
  left_join(celltypes_trailmaker, by = c("bc_wells" = "id")) %>% 
  mutate(celltype_trailmaker = ifelse(is.na(celltype_trailmaker), "missing", celltype_trailmaker))

integrated_pilotDexCortCtrl <- AddMetaData(
  object   = integrated_pilotDexCortCtrl,
  metadata = new_meta$celltype_trailmaker,
  col.name = "celltype_trailmaker"
)

save(
  integrated_pilotDexCortCtrl,
  file = "data/pilotDexCortCtrl/r_outputs/integrated_pilotDexCortCtrl_07.09.2025_v2.RData"
)


plot_umap_cluster_faceted(
  integrated_pilotDexCortCtrl,
  cluster_var = "celltype_trailmaker",
  # plot_title = "UMAP Projection Across All Samplest",
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
  n_row_legend = 2,
  filter_exclude = c("missing","unknown"),
  filter_var = "celltype_trailmaker"
)

plot_umap_cluster_faceted(
  integrated_pilotDexCortCtrl,
  cluster_var = "celltype_trailmaker",
  # plot_title = "UMAP Projection Across All Samplest",
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
  n_row_legend = 2,
  filter_include = c("neural_stem_cells"),
  filter_var = "celltype_trailmaker"
)


plot_umap_cluster_faceted(
  integrated_pilotDexCortCtrl,
  cluster_var = "celltype_trailmaker",
  # plot_title = "UMAP Projection Across All Samplest",
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
  n_row_legend = 2,
  filter_include = c("neural_stem_cells"),
  filter_var = "celltype_trailmaker"
)

integrated_pilotDexCortCtrl@meta.data$celltype_trailmaker
