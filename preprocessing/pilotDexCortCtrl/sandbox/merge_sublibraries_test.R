
# ##############################################################################
# ---- read data ----
# ##############################################################################
dge_dir <- "data/pilotDexCortCtrl/combined_BG60BLUE/all-well/DGE_filtered"
combined_BG60BLUE <- read_parse_to_seurat(dge_dir = dge_dir, project = "ParsePilot", use_gene_symbol = FALSE)

dge_dir <- "data/pilotDexCortCtrl/combined_BG60GREEN/all-well/DGE_filtered"
combined_BG60GREEN <- read_parse_to_seurat(dge_dir = dge_dir, project = "ParsePilot", use_gene_symbol = FALSE)

sample_metadata <- data.frame(
  sample_id = sprintf("%02d", 1:12),  # czyli "01", "02", ..., "12"
  well_label = paste0("A", 1:12),
  condition = c(
    rep("CTRL", 4),   # A1–A4
    rep("DEX", 4),    # A5–A8
    rep("CORT", 4)    # A9–A12
  )
)
# ##############################################################################
# ---- merged sublibraries ----
# ##############################################################################
prefixes <- c("BG60BLUE", "BG60GREEN")

merged_pilotDexCortCtrl <- merge_multiple_seurat_objects(
  seurat_list = list(combined_BG60BLUE, combined_BG60GREEN),
  project_name = "mergedParsePilot",
  prefixes = prefixes, 
  metadata_column_name = "sublibrary"
)


merged_pilotDexCortCtrl$treatment <- case_when(
  merged_pilotDexCortCtrl$bc1_well %in% c("A1", "A2", "A3", "A4")  ~ "CTRL",
  merged_pilotDexCortCtrl$bc1_well %in% c("A5", "A6", "A7", "A8")  ~ "DEX",
  merged_pilotDexCortCtrl$bc1_well %in% c("A9", "A10", "A11", "A12") ~ "CORT",
  TRUE ~ NA_character_
)


# ##############################################################################
# ---- run seruat analysis ----
# ##############################################################################
# Preprocessing step
merged_pilotDexCortCtrl <- seurat_precluster_preprocessing(
  seurat_obj = merged_pilotDexCortCtrl,
  n_variable_features = 2000,
  selection_method = "vst",
  # vars_to_regress = "treatment",
  verbose = TRUE
)

# Clustering + UMAP step with multiple resolutions
merged_pilotDexCortCtrl <- seurat_multiRes_clustering_umap(
  seurat_obj = merged_pilotDexCortCtrl,
  dims = 1:20,
  resolution = c(0.05, 0.1, 0.2, 0.3, 0.5, 0.4, 0.8, 1),
  k.param = 20,
  nn.method = "rann",
  algorithm = 1,
  reduction = "pca",
  do_umap = TRUE,
  umap_n_neighbors = 30,
  umap_min_dist = 0.3,
  umap_spread = 1.0,
  umap_metric = "cosine",
  verbose = TRUE,
  seed = 777
)


save(
  merged_pilotDexCortCtrl,
  file = "data/pilotDexCortCtrl/r_outputs/merged_pilotDexCortCtrl.RData"
)

# save(
#   merged_pilotDexCortCtrl,
#   file = "data/pilotDexCortCtrl/r_outputs/merged_pilotDexCortCtrl_covTreatment.RData",
#   compress = "xz",
#   version = 3
# )

merged_pilotDexCortCtrl@meta.data %>% head 



# ##############################################################################
# ---- UMAP visualization ----
# ##############################################################################

plot_umap_cluster_faceted(
  merged_pilotDexCortCtrl,
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
  cluster_var = "seurat_clusters_res0.1",
  plot_title = "UMAP Projection by Treatment",
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
