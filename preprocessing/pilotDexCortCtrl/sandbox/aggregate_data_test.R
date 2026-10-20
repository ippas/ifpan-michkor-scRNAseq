integrated_pilotDexCortCtrl <- seurat_merge_all_counts_layers(
  seurat_obj = integrated_pilotDexCortCtrl,
  assay = "RNA",
  remove_old = FALSE
)

agg_list <- aggregate_counts_by_cluster_sample(
  seurat_obj  = integrated_pilotDexCortCtrl,
  cluster_col = "seurat_clusters_res0.04",
  sample_col  = "sample_id",
  assay       = "RNA",
  layer       = "counts",  # <- Seurat v5
  agg         = "sum",
  return      = "list",
  progress    = TRUE
)

sample_map <- integrated_pilotDexCortCtrl@meta.data %>% distinct(sample_id, treatment)

# Wszystkie klastry, z komunikatami i paskiem postępu, bez gatekeep:
res <- edger_pseudobulk_anova_per_cluster(
  agg_list,
  sample_meta  = sample_map,
  group_levels = c("CTRL","DEX","CORT"),
  verbose      = TRUE,
  progress     = TRUE,
  gatekeep     = FALSE
)

res$cluster_0$f_test_global %>% filter(FDR < 0.01) %>% 
  filter(logFC.treatmentDEX > 1 | logFC.treatmentCORT > 1)


res %>%
  keep(~ "f_test_global" %in% names(.x)) %>%                
  imap(~ .x$f_test_global %>%
         filter(FDR < 0.05 & (logFC.treatmentDEX > 1 | logFC.treatmentCORT > 1))) %>% 
  write_cluster_list_to_xlsx(file_path = "data/pilotDexCortCtrl/r_outputs/statistics_results/edger_analysis_15.08.2025_v1/edgeR_AllTreatments_res0.04_FDR0.01log2FC1.xlsx")


res$cluster_0$CORT_vs_DEX %>% filter(gene_symbol == "FKBP5")

res %>% 
  keep(~ "CORT_vs_CTRL" %in% names(.x)) %>%                
  imap(~ .x$CORT_vs_CTRL %>%
         filter(FDR < 0.05) %>% 
         filter(logFC > 1 | logFC < -1)) %>% 
  write_cluster_list_to_xlsx(file_path = "data/pilotDexCortCtrl/r_outputs/statistics_results/edger_analysis_15.08.2025_v1/edgeR_CTRLvsCORT_res0.04_FDR0.01log2FC1.xlsx")

res %>% 
  keep(~ "DEX_vs_CTRL" %in% names(.x)) %>%                
  imap(~ .x$DEX_vs_CTRL %>%
         filter(FDR < 0.05) %>% 
         filter(logFC > 1 | logFC < -1)) %>% 
  write_cluster_list_to_xlsx(file_path = "data/pilotDexCortCtrl/r_outputs/statistics_results/edger_analysis_15.08.2025_v1/edgeR_CTRLvsDEX_res0.04_FDR0.01log2FC1.xlsx")

res %>% 
  keep(~ "CORT_vs_DEX" %in% names(.x)) %>%                
  imap(~ .x$DEX_vs_CTRL %>%
         filter(FDR < 0.05) %>% 
         filter(logFC > 1 | logFC < -1)) %>% 
  write_cluster_list_to_xlsx(file_path = "data/pilotDexCortCtrl/r_outputs/statistics_results/edger_analysis_15.08.2025_v1/edgeR_CORTvsDEX_res0.04_FDR0.01log2FC1.xlsx")

# ##############################################################################
# ---- for res 0.13 ----
# ##############################################################################
agg_list <- aggregate_counts_by_cluster_sample(
  seurat_obj  = integrated_pilotDexCortCtrl,
  cluster_col = "seurat_clusters_res0.13",
  sample_col  = "sample_id",
  assay       = "RNA",
  layer       = "counts",  # <- Seurat v5
  agg         = "sum",
  return      = "list",
  progress    = TRUE
)

sample_map <- integrated_pilotDexCortCtrl@meta.data %>% distinct(sample_id, treatment)

# Wszystkie klastry, z komunikatami i paskiem postępu, bez gatekeep:
res <- edger_pseudobulk_anova_per_cluster(
  agg_list,
  sample_meta  = sample_map,
  group_levels = c("CTRL","DEX","CORT"),
  verbose      = TRUE,
  progress     = TRUE,
  gatekeep     = FALSE
)

res$cluster_0$f_test_global %>% filter(FDR < 0.01) %>% 
  filter(logFC.treatmentDEX > 1 | logFC.treatmentCORT > 1)


res %>%
  keep(~ "f_test_global" %in% names(.x)) %>%                
  imap(~ .x$f_test_global %>%
         filter(FDR < 0.05 & (logFC.treatmentDEX > 1 | logFC.treatmentCORT > 1))) %>% 
  write_cluster_list_to_xlsx(file_path = "data/pilotDexCortCtrl/r_outputs/statistics_results/edger_analysis_15.08.2025_v1/edgeR_AllTreatments_res0.13_FDR0.01log2FC1.xlsx")


res$cluster_0$CORT_vs_DEX %>% filter(gene_symbol == "FKBP5")

res %>% 
  keep(~ "CORT_vs_CTRL" %in% names(.x)) %>%                
  imap(~ .x$CORT_vs_CTRL %>%
         filter(FDR < 0.05) %>% 
         filter(logFC > 1 | logFC < -1)) %>% 
  write_cluster_list_to_xlsx(file_path = "data/pilotDexCortCtrl/r_outputs/statistics_results/edger_analysis_15.08.2025_v1/edgeR_CTRLvsCORT_res0.13_FDR0.01log2FC1.xlsx")

res %>% 
  keep(~ "DEX_vs_CTRL" %in% names(.x)) %>%                
  imap(~ .x$DEX_vs_CTRL %>%
         filter(FDR < 0.05) %>% 
         filter(logFC > 1 | logFC < -1)) %>% 
  write_cluster_list_to_xlsx(file_path = "data/pilotDexCortCtrl/r_outputs/statistics_results/edger_analysis_15.08.2025_v1/edgeR_CTRLvsDEX_res0.13_FDR0.01log2FC1.xlsx")

res %>% 
  keep(~ "CORT_vs_DEX" %in% names(.x)) %>%                
  imap(~ .x$DEX_vs_CTRL %>%
         filter(FDR < 0.05) %>% 
         filter(logFC > 1 | logFC < -1)) %>% 
  write_cluster_list_to_xlsx(file_path = "data/pilotDexCortCtrl/r_outputs/statistics_results/edger_analysis_15.08.2025_v1/edgeR_CORTvsDEX_res0.13_FDR0.01log2FC1.xlsx")



# ##############################################################################
# ---- for res 0.2 ----
# ##############################################################################
agg_list <- aggregate_counts_by_cluster_sample(
  seurat_obj  = integrated_pilotDexCortCtrl,
  cluster_col = "seurat_clusters_res0.2",
  sample_col  = "sample_id",
  assay       = "RNA",
  layer       = "counts",  # <- Seurat v5
  agg         = "sum",
  return      = "list",
  progress    = TRUE
)

sample_map <- integrated_pilotDexCortCtrl@meta.data %>% distinct(sample_id, treatment)

# Wszystkie klastry, z komunikatami i paskiem postępu, bez gatekeep:
res <- edger_pseudobulk_anova_per_cluster(
  agg_list,
  sample_meta  = sample_map,
  group_levels = c("CTRL","DEX","CORT"),
  verbose      = TRUE,
  progress     = TRUE,
  gatekeep     = FALSE
)

res$cluster_0$f_test_global %>% filter(FDR < 0.01) %>% 
  filter(logFC.treatmentDEX > 1 | logFC.treatmentCORT > 1)


res %>%
  keep(~ "f_test_global" %in% names(.x)) %>%                
  imap(~ .x$f_test_global %>%
         filter(FDR < 0.05 & (logFC.treatmentDEX > 1 | logFC.treatmentCORT > 1))) %>% 
  write_cluster_list_to_xlsx(file_path = "data/pilotDexCortCtrl/r_outputs/statistics_results/edger_analysis_15.08.2025_v1/edgeR_AllTreatments_res0.2_FDR0.01log2FC1.xlsx")


res$cluster_0$CORT_vs_DEX %>% filter(gene_symbol == "FKBP5")

res %>% 
  keep(~ "CORT_vs_CTRL" %in% names(.x)) %>%                
  imap(~ .x$CORT_vs_CTRL %>%
         filter(FDR < 0.05) %>% 
         filter(logFC > 1 | logFC < -1)) %>% 
  write_cluster_list_to_xlsx(file_path = "data/pilotDexCortCtrl/r_outputs/statistics_results/edger_analysis_15.08.2025_v1/edgeR_CTRLvsCORT_res0.2_FDR0.01log2FC1.xlsx")

res %>% 
  keep(~ "DEX_vs_CTRL" %in% names(.x)) %>%                
  imap(~ .x$DEX_vs_CTRL %>%
         filter(FDR < 0.05) %>% 
         filter(logFC > 1 | logFC < -1)) %>% 
  write_cluster_list_to_xlsx(file_path = "data/pilotDexCortCtrl/r_outputs/statistics_results/edger_analysis_15.08.2025_v1/edgeR_CTRLvsDEX_res0.2_FDR0.01log2FC1.xlsx")

res %>% 
  keep(~ "CORT_vs_DEX" %in% names(.x)) %>%                
  imap(~ .x$DEX_vs_CTRL %>%
         filter(FDR < 0.05) %>% 
         filter(logFC > 1 | logFC < -1)) %>% 
  write_cluster_list_to_xlsx(file_path = "data/pilotDexCortCtrl/r_outputs/statistics_results/edger_analysis_15.08.2025_v1/edgeR_CORTvsDEX_res0.2_FDR0.01log2FC1.xlsx")
