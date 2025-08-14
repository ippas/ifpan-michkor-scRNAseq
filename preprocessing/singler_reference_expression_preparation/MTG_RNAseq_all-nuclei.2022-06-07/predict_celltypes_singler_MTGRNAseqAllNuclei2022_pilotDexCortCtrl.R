
# ##############################################################################
# ---- load data ----
# ##############################################################################
load("data/reference_celltype_expression_processed/sce_ref_MTGRNAseqAllNuclei2022.RData")

# ##############################################################################
# ---- prepare data do sce ----
# ##############################################################################
sce_integrated_pilotDexCortCtrl_multi <- prepare_expression_query_multi_singler_v1(
  seurat_obj = integrated_pilotDexCortCtrl,
  assay = "RNA",
  data_layer_regex = "^data",
  normalize_if_missing = TRUE,
  gene_name_from_rownames = TRUE,
  gene_sep = "-",
  agg_fun_for_duplicates = "mean",
  cluster_cols = c( "seurat_clusters_res0.04", "seurat_clusters_res0.13", "seurat_clusters_res0.2"),
  alias_first_as = "cluster",     # SingleR-friendly alias
  save_rds = "data/reference_celltype_expression_processed/sce_integrated_pilotDexCortCtrl_res_0.04_0.13_0.2.rds",
  verbose = TRUE
)


colnames(sce_ref_MTGRNAseqAllNuclei2022)
# colData(sce_ref_MTGRNAseqAllNuclei2022) %>% head
# 
# sce_ref_MTGRNAseqAllNuclei2022$CA_subclass_label %>% table
# 
# sce_ref_MTGRNAseqAllNuclei2022$GA_cluster_label %>% table



predict_MTGRNAseqAllNuclei2022_pilotDexCortCtrl <- predict_celltypes_singler_multi_v1(
  sce_query = sce_integrated_pilotDexCortCtrl_multi,
  sce_ref   = sce_ref_MTGRNAseqAllNuclei2022,
  label_col = "CA_subclass_label",
  cluster_cols = c( "seurat_clusters_res0.04", "seurat_clusters_res0.13", "seurat_clusters_res0.2"),
  min_common_genes = 1000,
  save_dir = "data/processed_reference_celltype_expression/singler_MTGRNAseqAllNuclei2022_pilotDexCortCtrl",
  verbose = TRUE
)



