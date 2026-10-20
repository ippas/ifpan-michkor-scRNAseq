# ##############################################################################
# ---- prepare data ----
# ##############################################################################

panglaoDB_all_markers <- read.table(
  "data/celltype_references/PanglaoDB/PanglaoDB_markers_27_Mar_2020.tsv",
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  # check.names = FALSE,
  row.names = NULL
)

panglaoDB_all_markers %>% 
  filter(organ == "Brain") %>% .$cell.type %>% unique

panglaoDB_all_markers %>% 
  select(c(official.gene.symbol, cell.type)) %>% 
  mutate(cell.type = str_replace_all(cell.type, " ", "_")) %>% 
  group_by(cell.type) %>% 
  split(.$cell.type) %>% 
  lapply(function(x) unique(x$official.gene.symbol)) -> panglaoDB_list_allMarkers


panglaoDB_all_markers %>% 
  filter(organ %in% c("Brain", "Vasculature")) %>% 
  # filter(organ == "Vasculature") %>%
  filter(!(cell.type %in% c("Endothelial cells", "Endothelial cells (aorta)", "Endothelial cells (blood brain barrier)" ))) %>% 
  select(c(official.gene.symbol, cell.type)) %>% 
  mutate(cell.type = str_replace_all(cell.type, " ", "_")) %>% 
  group_by(cell.type) %>% 
  split(.$cell.type) %>% 
  lapply(function(x) unique(x$official.gene.symbol)) -> panglaoDB_list_brainVasculatureMarkers

# ##############################################################################
# ---- run analysis ----
# ##############################################################################

avg_mat <- prepare_proteincoding_avg_expr(
  seurat_obj       = integrated_pilotDexCortCtrl,
  cluster_col      = "seurat_clusters_res0.13",
  sample_col       = "sample_id",
  assay            = "RNA",
  layer            = "counts",
  biomart_version  = "114",
  verbose          = TRUE
)


# --- wyliczenie specyficzności ---
gene_specificity <- calculate_gene_specificity(avg_mat)

fgsea_out <- run_fgsea_on_score_matrix(
  # mat        =  tmp %>% na.omit(),
  mat = gene_specificity$combined,
  gene_lists = panglaoDB_list_allMarkers,
  # gene_lists = panglaoDB_list_allMarkers[c("Pericytes",   "Astrocytes", "Bergmann_glia", "Microglia",
  #                                          "Oligodendrocyte_progenitor_cells", "Radial_glia_cells",
  #                                          "Satellite_glial_cells", "Tanycytes",
  #                                          "Neurons", "Immature_neurons", "Neural_stem/precursor_cells")],   # or your custom marker sets
  nperm      = 100,
  do_fisher = TRUE,
  export_xlsx = TRUE,
  xlsx_path = "data/pilotDexCortCtrl/r_outputs/gsea_results/gsea_CellTypePangloDB_proteinCodingBiomartv112_res0.13.xlsx",
  write_first_by_NES = TRUE,
  write_first_by_pval = TRUE,
  write_first_by_padj = TRUE,
  seed = 777
)

fgsea_out$results$cluster_0 %>% select(-genes)



res_multi_aprioriSelectedCellType2 <- run_celltype_fgseaEnrichemnt_multClusters(
  seurat_obj   = integrated_pilotDexCortCtrl,
  biomart_version = "112",
  cluster_cols = c(
    # "seurat_clusters_res0.04", 
                   "seurat_clusters_res0.13" 
                   # "seurat_clusters_res0.2",
                   # "seurat_clusters_res0.5"
                   ),
  # gene_lists   = panglaoDB_list_allMarkers[c(
  #   "Pericytes","Astrocytes","Bergmann_glia","Microglia",
  #   "Oligodendrocyte_progenitor_cells","Radial_glia_cells",
  #   "Satellite_glial_cells","Tanycytes",
  #   # "Neurons",
  #   "Oligodendrocyte",
  #   "Immature_neurons","Neural_stem/precursor_cells"
  # )],
  gene_lists =  panglaoDB_list_brainVasculatureMarkers,
  nperm        = 10000,
  do_fisher    = TRUE,
  export_xlsx  = TRUE,
  select_samples = c(
    "BG60BLUE_A1", "BG60BLUE_A2", "BG60BLUE_A3", "BG60BLUE_A4",
    "BG60GREEN_A1", "BG60GREEN_A2", "BG60GREEN_A3", "BG60GREEN_A4"
  ),
  xlsx_path    = "data/pilotDexCortCtrl/r_outputs/gsea_results/gsea_aprioriSelectedCellTypePanglaoDB_{cluster}.xlsx"
)

res_multi_allCellType <- run_celltype_fgseaEnrichemnt_multClusters(
  seurat_obj   = integrated_pilotDexCortCtrl,
  biomart_version = "112",
  cluster_cols = c("seurat_clusters_res0.04", 
                   "seurat_clusters_res0.13", 
                   "seurat_clusters_res0.2",
                   "seurat_clusters_res0.5"),
  gene_lists   = panglaoDB_list_allMarkers,
  nperm        = 100000,
  do_fisher    = TRUE,
  export_xlsx  = TRUE,
  xlsx_path    = "data/pilotDexCortCtrl/r_outputs/gsea_results/gsea_allCellTypePanglaoDB_{cluster}.xlsx"
)

integrated_pilotDexCortCtrl@meta.data$seurat_clusters_res0.5


res_multi2$seurat_clusters_res0.13$cell_counts_vec



fgsea_heatmap_log10p_from_results(
  fgsea_results  = res_multi_aprioriSelectedCellType$seurat_clusters_res0.04$fgsea,
  path_col       = "predict_cell_type",
  pval_col       = "pval",
  axis_x         = "cluster",
  axis_y         = "celltype",
  outline_max    = TRUE,
  outline_reduce = "by_col",
  outline_mode   = "max_nes",
  outer_border_lwd = 3,
  outline_lwd = 3,
  numbers_thresh = 3,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  color_palette  = list(low = "white", mid = "#E47A66", high = "#6B2113"),
  title_prefix = "resolution = 0.04;   ",
  celltype_rename = list("Oligodendrocyte_progenitor_cells" = "OPC", 
                         "Neural_stem/precursor_cells" = "NPCs")
  )

fgsea_heatmap_log10p_from_results(
  fgsea_results  = res_multi_aprioriSelectedCellType$seurat_clusters_res0.13$fgsea,
  path_col       = "predict_cell_type",
  pval_col       = "pval",
  axis_x         = "cluster",
  axis_y         = "celltype",
  outline_max    = TRUE,
  outline_reduce = "by_col",
  outline_mode   = "max_nes",
  outer_border_lwd = 3,
  outline_lwd = 3,
  numbers_thresh = 3,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  color_palette  = list(low = "white", mid = "#E47A66", high = "#6B2113"),
  title_prefix = "resolution = 0.13;   ",
  celltype_rename = list("Oligodendrocyte_progenitor_cells" = "OPC", 
                         "Neural_stem/precursor_cells" = "NPCs")
)

fgsea_heatmap_log10p_from_results(
  fgsea_results  = res_multi_aprioriSelectedCellType$seurat_clusters_res0.2$fgsea,
  path_col       = "predict_cell_type",
  pval_col       = "pval",
  axis_x         = "cluster",
  axis_y         = "celltype",
  outline_max    = TRUE,
  outline_reduce = "by_col",
  outline_mode   = "max_nes",
  outer_border_lwd = 3,
  outline_lwd = 3,
  numbers_thresh = 3,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  color_palette  = list(low = "white", mid = "#E47A66", high = "#6B2113"),
  title_prefix = "resolution = 0.2;   ",
  celltype_rename = list("Oligodendrocyte_progenitor_cells" = "OPC", 
                         "Neural_stem/precursor_cells" = "NPCs")
)

fgsea_heatmap_log10p_from_results(
  fgsea_results  = res_multi_aprioriSelectedCellType$seurat_clusters_res0.5$fgsea,
  path_col       = "predict_cell_type",
  pval_col       = "pval",
  axis_x         = "cluster",
  axis_y         = "celltype",
  outline_max    = TRUE,
  outline_reduce = "by_col",
  outline_mode   = "max_nes",
  outer_border_lwd = 3,
  outline_lwd = 3,
  numbers_thresh = 3,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  color_palette  = list(low = "white", mid = "#E47A66", high = "#6B2113"),
  title_prefix = "resolution = 0.5;   ",
  celltype_rename = list("Oligodendrocyte_progenitor_cells" = "OPC", 
                         "Neural_stem/precursor_cells" = "NPCs")
)


# ##############################################################################
fgsea_heatmap_log10p_from_results(
  fgsea_results  = res_multi_allCellType$seurat_clusters_res0.13$fgsea,
  path_col       = "predict_cell_type",
  pval_col       = "pval",
  axis_x         = "cluster",
  axis_y         = "celltype",
  outline_max    = TRUE,
  outline_reduce = "by_col",
  outline_mode   = "max_nes",
  outer_border_lwd = 3,
  outline_lwd = 3,
  numbers_thresh = 3,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  color_palette  = list(low = "white", mid = "#E47A66", high = "#6B2113"),
  title_prefix = "resolution = 0.13;   ",
  celltype_rename = list("Oligodendrocyte_progenitor_cells" = "OPC", 
                         "Neural_stem/precursor_cells" = "NPCs")
)













# 
# fgsea_plot_NES_heatmap(
#   fgsea_results = res_multi$seurat_clusters_res0.13$fgsea,              # albo res_multi[["seurat_clusters_res0.13"]]$fgsea
#   path_col      = "predict_cell_type",    # albo "pathway"
#   padj_max      = 0.05,                   # np. pokaż tylko istotne
#   scale_rows    = TRUE                    # standaryzacja wierszy
# )
# 
# res_multi$seurat_clusters_res0.13$fgse$results$cluster_0
# 
# fgsea_heatmap_from_results(
#   fgsea_results = res_multi$seurat_clusters_res0.13$fgsea,
#   path_col      = "predict_cell_type",
#   value_col     = "NES",
#   axis_x        = "celltype",
#   axis_y        = "cluster",
#   na_color = "gray",
#   zscore_row    = F,
#   zscore_limit  = c(-4, 4)   # przycięcie
# )
# 
# res_multi$seurat_clusters_res0.13$fgsea$results %>%
#   purrr::imap_dfr(~ .x %>% dplyr::filter(is.na(NES)) %>%
#                     dplyr::mutate(cluster = .y))
# res_multi$seurat_clusters_res0.13$fgsea$rank_stats
# 
# 
# 
# fgsea_heatmap_from_results(
#   fgsea_results = res_multi$seurat_clusters_res0.13$fgsea,
#   path_col      = "predict_cell_type",
#   value_col     = "NES",
#   axis_x        = "celltype",
#   axis_y        = "cluster",
#   zscore_row    = FALSE,       # kolory = surowy NES (polecane do porównań)
#   numbers       = TRUE,
#   numbers_thresh= 1.4,
#   outline_max   = TRUE,
#   numbers_bold = TRUE,
#   gap_factor      = 0.97,  
#   celltype_rename = list("Oligodendrocyte_progenitor_cells" = "OPC", 
#                          "Neural_stem/precursor_cells" = "NPCs")
#   # heatmap_colors = list(low = "blue", mid = "white", high = "firebrick")
#   # row_dend_side_space = grid::unit(15, "mm")
# )
# 
# res_multi$seurat_clusters_res0.13$fgsea$results$cluster_7
# # 
# 
# 
# fgsea_heatmap_log10p_from_results(
#   fgsea_results  = res_multi$seurat_clusters_res0.13$fgsea,
#   path_col       = "predict_cell_type",
#   pval_col       = "pval",
#   numbers        = TRUE,
#   numbers_digits = 2,
#   nes_digits     = 2,
#   numbers_size_pval = 11,
#   numbers_size_nes  = 11,
#   numbers_bold   = TRUE,
#   outline_max    = TRUE,
#   outline_mode   = "max_nes",
#   
#   celltype_rename = list("Oligodendrocyte_progenitor_cells" = "OPC", 
#                          "Neural_stem/precursor_cells" = "NPCs")
# )
# 
# 
# fgsea_heatmap_log10p_from_results(
#   fgsea_results  = res_multi$seurat_clusters_res0.13$fgsea,
#   path_col       = "predict_cell_type",
#   pval_col       = "pval",
#   axis_x         = "celltype",
#   axis_y         = "cluster",
#   outline_max    = TRUE,
#   outline_mode   = "max_nes",
#   celltype_rename = list("Oligodendrocyte_progenitor_cells" = "OPC", 
#                          "Neural_stem/precursor_cells" = "NPCs"),
#   outline_reduce = "by_row"
# )