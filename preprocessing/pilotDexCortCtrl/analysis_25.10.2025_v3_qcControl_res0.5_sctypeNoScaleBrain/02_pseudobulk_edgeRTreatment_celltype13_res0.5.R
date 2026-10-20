integrated_pilotDexCortCtrl_filtered <- seurat_merge_all_counts_layers(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  assay = "RNA",
  remove_old = FALSE
)

agg_list <- aggregate_counts_by_cluster_sample(
  seurat_obj  = integrated_pilotDexCortCtrl_filtered,
  cluster_col = "res0.3_RNA_noScaled_brain",
  sample_col  = "bc1_well",
  assay       = "RNA",
  layer       = "counts",  # <- Seurat v5
  agg         = "sum",
  return      = "list",
  progress    = TRUE
)

sample_map <- integrated_pilotDexCortCtrl_filtered@meta.data %>% distinct(bc1_well, treatment)

# Wszystkie klastry, z komunikatami i paskiem postępu, bez gatekeep:
edgeR_anovaOneFactor_treatment_pseudobulk_res0.5_clusters13_results <- edger_pseudobulk_anova_per_cluster_v2(
  agg_list,
  sample_meta  = sample_map,
  sample_col   = "bc1_well",
  group_levels = c("CTRL","DEX","CORT"),
  verbose      = TRUE,
  progress     = TRUE,
  gatekeep     = F
)



integrated_pilotDexCortCtrl_filtered <- seurat_merge_all_counts_layers(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  assay = "RNA",
  remove_old = FALSE
)

agg_list <- aggregate_counts_by_cluster_sample(
  seurat_obj  = integrated_pilotDexCortCtrl_filtered,
  cluster_col = "res0.3_RNA_noScaled_allMarkers",
  sample_col  = "bc1_well",
  assay       = "RNA",
  layer       = "counts",  # <- Seurat v5
  agg         = "sum",
  return      = "list",
  progress    = TRUE
)

sample_map <- integrated_pilotDexCortCtrl_filtered@meta.data %>% distinct(bc1_well, treatment)

# Wszystkie klastry, z komunikatami i paskiem postępu, bez gatekeep:
edgeR_anovaOneFactor_treatment_pseudobulk_res0.3_clusters16cells15_results <- edger_pseudobulk_anova_per_cluster_v2(
  agg_list,
  sample_meta  = sample_map,
  sample_col   = "bc1_well",
  group_levels = c("CTRL","DEX","CORT"),
  verbose      = TRUE,
  progress     = TRUE,
  gatekeep     = F
)




edgeR_anovaOneFactor_treatment_pseudobulk_res0.3_clusters16cells15_results$cluster_Oligodendrocytes$f_test_global


edgeR_anovaOneFactor_treatment_pseudobulk_res0.3_clusters16cells15_results$`cluster_Oligodendrocyte precursor cells`$f_test_global %>% filter(FDR < 0.01)


edgeR_anovaOneFactor_treatment_pseudobulk_res0.3_clusters16cells15_results$cluster_Oligodendrocytes$f_test_global %>% filter(FDR < 0.01) %>% 
  .$gene_symbol


edgeR_anovaOneFactor_treatment_pseudobulk_res0.3_clusters16cells15_results$`cluster_Immature neurons`$DEX_vs_CTRL %>% filter(FDR < 0.01) %>% filter(abs(logFC) > 1)

edgeR_anovaOneFactor_treatment_pseudobulk_res0.3_clusters16cells15_results$`cluster_Immature neurons`$CORT_vs_CTRL %>% filter(FDR < 0.01) %>% filter(abs(logFC) > 1)



edgeR_anovaOneFactor_treatment_pseudobulk_res0.3_clusters16cells15_results$cluster_Oligodendrocytes$DEX_vs_CTRL %>% filter(FDR < 0.01) %>% filter(abs(logFC) > 1)

edgeR_anovaOneFactor_treatment_pseudobulk_res0.3_clusters16cells15_results$cluster_Oligodendrocytes$CORT_vs_CTRL %>% filter(FDR < 0.01) %>% filter(abs(logFC) > 1)


edgeR_anovaOneFactor_treatment_pseudobulk_res0.3_clusters16cells15_results$`cluster_Radial glial cells`$DEX_vs_CTRL %>% filter(FDR < 0.01) %>% filter(abs(logFC) > 1)

edgeR_anovaOneFactor_treatment_pseudobulk_res0.3_clusters16cells15_results$`cluster_Radial glial cells`$CORT_vs_CTRL %>% filter(FDR < 0.01) %>% filter(abs(logFC) > 1)



lapply(
  edgeR_anovaOneFactor_treatment_pseudobulk_res0.3_clusters16cells15_results,
  function(cluster) {
    if (!is.null(cluster$f_test_global)) {
      cluster$f_test_global %>%
        dplyr::filter(FDR < 0.01)
    } else {
      NULL
    }
  }
) %>% lapply(., function(x){
  x %>% .$gene_symbol
}) %>% 
  unname %>% unlist %>% table %>% as.data.frame() %>% arrange(desc(Freq))



# ##############################################################################
agg_list <- aggregate_counts_by_cluster_sample(
  seurat_obj  = integrated_pilotDexCortCtrl_filtered,
  cluster_col = "res0.3_RNA_noScaled_brain",
  sample_col  = "bc1_well",
  assay       = "RNA",
  layer       = "counts",  # <- Seurat v5
  agg         = "sum",
  return      = "list",
  progress    = TRUE
)

sample_map <- integrated_pilotDexCortCtrl_filtered@meta.data %>% distinct(bc1_well, treatment)

# Wszystkie klastry, z komunikatami i paskiem postępu, bez gatekeep:
edgeR_anovaOneFactor_treatment_pseudobulk_res0.4_clusters16_celltypes11_results <- edger_pseudobulk_anova_per_cluster_v2(
  agg_list,
  sample_meta  = sample_map,
  sample_col   = "bc1_well",
  group_levels = c("CTRL","DEX","CORT"),
  verbose      = TRUE,
  progress     = TRUE,
  gatekeep     = F
)


integrated_pilotDexCortCtrl_filtered <- seurat_merge_all_counts_layers(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  assay = "RNA",
  remove_old = FALSE
)

agg_list <- aggregate_counts_by_cluster_sample(
  seurat_obj  = integrated_pilotDexCortCtrl_filtered,
  cluster_col = "res0.3_RNA_noScaled_allMarkers",
  sample_col  = "bc1_well",
  assay       = "RNA",
  layer       = "counts",  # <- Seurat v5
  agg         = "sum",
  return      = "list",
  progress    = TRUE
)


                                                              