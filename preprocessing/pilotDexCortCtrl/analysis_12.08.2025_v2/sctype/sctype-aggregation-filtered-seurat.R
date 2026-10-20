combine_sctype_results <- function(
    seurat_obj,
    res,
    all_markers_list,
    brain_markers_list,
    assay = c("RNA", "integrated"),
    layer = "data",
    biomart_version = "v112"
) {
  assay <- match.arg(assay)
  
  # --- przygotowanie macierzy ekspresji ---
  exprMat_cluster_meanData <- prepare_proteinCodingGenesMatrix_fromSeurat_v2(
    seurat_obj,
    assay           = assay,
    layer           = layer,
    biomart_version = biomart_version,
    method          = "mean",
    aggregate_cluster = TRUE,
    group.by = paste0("seurat_clusters_res", res)
  )
  
  # --- liczenie sctype ---
  sctype_scaled_all      <- sctype_score(exprMat_cluster_meanData, TRUE,  gs = all_markers_list,   gene_names_to_uppercase = TRUE)
  sctype_noscaled_all    <- sctype_score(exprMat_cluster_meanData, FALSE, gs = all_markers_list,   gene_names_to_uppercase = TRUE)
  sctype_scaled_brain    <- sctype_score(exprMat_cluster_meanData, TRUE,  gs = brain_markers_list, gene_names_to_uppercase = TRUE)
  sctype_noscaled_brain  <- sctype_score(exprMat_cluster_meanData, FALSE, gs = brain_markers_list, gene_names_to_uppercase = TRUE)
  
  # --- helper do top celltype ---
  get_top_celltype <- function(sctype_mat) {
    if (is.null(sctype_mat)) return(rep(NA_character_, length = ncol(sctype_mat)))
    apply(sctype_mat, 2, function(x) rownames(sctype_mat)[which.max(x)])
  }
  
  # --- zwracamy wynikowy df ---
  tibble(
    cluster                          = names(get_top_celltype(sctype_scaled_all)),
    cluster_id                       = gsub("cluster_", "", names(get_top_celltype(sctype_scaled_all))),
    
    !!paste0("res", res, "_", assay, "_scaled_allMarkers")   := unname(get_top_celltype(sctype_scaled_all)),
    !!paste0("res", res, "_", assay, "_noScaled_allMarkers") := unname(get_top_celltype(sctype_noscaled_all)),
    !!paste0("res", res, "_", assay, "_scaled_brain")        := unname(get_top_celltype(sctype_scaled_brain)),
    !!paste0("res", res, "_", assay, "_noScaled_brain")      := unname(get_top_celltype(sctype_noscaled_brain))
  ) %>%
    arrange(as.numeric(cluster_id))
}

df_top_sctype_res0.5_RNA <- combine_sctype_results(
  integrated_pilotDexCortCtrl_filtered,
  res = "0.5",
  all_markers_list = all_markers_list,
  brain_markers_list = brain_markers_list,
  assay = "RNA",
  biomart_version = "v114"
)

df_top_sctype_res0.13_RNA <- combine_sctype_results(
  integrated_pilotDexCortCtrl_filtered,
  res = "0.13",
  all_markers_list = all_markers_list,
  brain_markers_list = brain_markers_list,
  assay = "RNA",
  biomart_version = "v114"
)

df_top_sctype_res0.5_Integrated <- combine_sctype_results(
  integrated_pilotDexCortCtrl_filtered,
  res = "0.5",
  all_markers_list = all_markers_list,
  brain_markers_list = brain_markers_list,
  assay = "integrated",
  biomart_version = "v114"
)

df_top_sctype_res0.13_Integrated <- combine_sctype_results(
  integrated_pilotDexCortCtrl_filtered,
  res = "0.13",
  all_markers_list = all_markers_list,
  brain_markers_list = brain_markers_list,
  assay = "integrated", 
  biomart_version = "v114"
)


# --- res0.5: łączymy po cluster i cluster_id ---
df_top_sctype_res0.5_combined <- df_top_sctype_res0.5_RNA %>%
  full_join(
    df_top_sctype_res0.5_Integrated %>% select(-cluster),
    by = "cluster_id"
  ) %>%
  arrange(as.numeric(cluster_id))

# --- res0.13: łączymy po cluster i cluster_id ---
df_top_sctype_res0.13_combined <- df_top_sctype_res0.13_RNA %>%
  full_join(
    df_top_sctype_res0.13_Integrated %>% select(-cluster),
    by = "cluster_id"
  ) %>%
  arrange(as.numeric(cluster_id))



# ---- remove tmp tables ----
rm(
  df_top_sctype_res0.5_RNA,
  df_top_sctype_res0.13_RNA,
  df_top_sctype_res0.5_Integrated,
  df_top_sctype_res0.13_Integrated
)


# ##############################################################################
# ---- add to metadata ----
# ##############################################################################
# integrated_pilotDexCortCtrl_filtered <- add_df_to_seurat_metadata(
#   seurat_obj      = integrated_pilotDexCortCtrl_filtered,
#   df              = df_top_sctype_res0.5_combined,
#   join_by_meta    = "seurat_clusters_res0.5",
#   join_by_df      = "cluster_id",
#   exclude_columns = "cluster"  # pomijamy kolumnę 'cluster'
# )
# 

integrated_pilotDexCortCtrl_filtered <- add_df_to_seurat_metadata(
  seurat_obj      = integrated_pilotDexCortCtrl_filtered,
  df              = df_top_sctype_res0.13_combined,
  join_by_meta    = "seurat_clusters_res0.13",
  join_by_df      = "cluster_id",
  exclude_columns = "cluster"
)

integrated_pilotDexCortCtrl@meta.data


# ##############################################################################
# ---- visualization on umap ----
# ##############################################################################

plot_umap_cluster_faceted(
  integrated_pilotDexCortCtrl_filtered,
  cluster_var = "res0.5_integrated_noScaled_brain",
  plot_title = "UMAP, res0.5; noScalled BrainMarkers",
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
) -> p1


plot_umap_cluster_faceted(
  integrated_pilotDexCortCtrl,
  cluster_var = "res0.5_integrated_scaled_brain",
  plot_title = "UMAP, res0.5; scaled BrainMarkers",
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
)-> p2


plot_umap_cluster_faceted(
  integrated_pilotDexCortCtrl,
  cluster_var = "res0.5_integrated_noScaled_allMarkers",
  plot_title = "UMAP, res0.5; noScalled AllMarkers",
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
) -> p3


plot_umap_cluster_faceted(
  integrated_pilotDexCortCtrl,
  cluster_var = "res0.5_integrated_scaled_allMarkers",
  plot_title = "UMAP, res0.5; scalled AllMarkers",
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
) -> p4

(p3 + p4) / (p1 + p2)



plot_umap_cluster_faceted(
  integrated_pilotDexCortCtrl,
  cluster_var = "res0.5_RNA_noScaled_allMarkers",
  plot_title = "UMAP, res0.5; allProteinCoding; noScalled AllMarkers",
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
) -> p5

plot_umap_cluster_faceted(
  integrated_pilotDexCortCtrl,
  cluster_var = "res0.5_RNA_noScaled_brain",
  plot_title = "UMAP, res0.5; allProteinCoding; noScalled brainMarkers",
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
) -> p6



(p3 + p1) / (p5 + p6)



plot_umap_cluster_faceted(
  integrated_pilotDexCortCtrl,
  cluster_var = "res0.13_RNA_noScaled_brain",
  plot_title = "UMAP, res0.13; noScalled; Brain + Hippocampus markers",
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
) -> p1

plot_umap_cluster_faceted(
  integrated_pilotDexCortCtrl,
  cluster_var = "res0.13_RNA_noScaled_allMarkers",
  plot_title = "UMAP, res0.13; noScalled AllMarkers",
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
) -> p2

p1 + p2
