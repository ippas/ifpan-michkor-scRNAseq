rawSumCounts_res0.3_RNA_noScaledBrain <- aggregate_counts_by_cluster_sample(
  seurat_obj  = integrated_pilotDexCortCtrl_filtered,
  cluster_col = "res0.3_RNA_noScaled_brain",
  sample_col  = "bc1_well",
  assay       = "RNA",
  layer       = "counts",  # <- Seurat v5
  agg         = "sum",
  return      = "list",
  genes_keep  = hgnc_symbols_vector_v110,
  progress    = TRUE
)


normCPM_res0.3_RNA_noScaledBrain <- normalize_CPM_perCluster(
  count_list = rawSumCounts_res0.3_RNA_noScaledBrain,
  scale_factor = 1e6,
  progress = TRUE
)

log2CPM_res0.3_RNA_noScaledBrain <- log2_transform_perCluster(
  count_list = normCPM_res0.3_RNA_noScaledBrain,
  pseudocount = 1,
  progress = TRUE
)

quantileNormalize_log2CPM_res0.3_RNA_noScaledBrain <- quantile_normalize_perCluster(
  count_list = log2CPM_res0.3_RNA_noScaledBrain,
  prefix2colname = "QN",
  progress = TRUE
)

zcore_QN_log2CPM_res0.3_RNA_noScaledBrain <- zscore_perCluster(
  count_list = quantileNormalize_log2CPM_res0.3_RNA_noScaledBrain,
  center = TRUE,
  scale = TRUE,
  progress = TRUE
)

mean_CPM_perClusterTreatment <- aggregate_CPMmean_list(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  cluster_col = "res0.3_RNA_noScaled_brain",
  sample_col  = "bc1_well",
  assay         = "RNA",
  layer         = "counts",
  select_samples = NULL,
  genes_keep     = hgnc_symbols_vector_v110,
  scale_factor   = 1e6,
  progress       = TRUE
)


# rawSumCounts_res0.3_RNA_noScaledBrain$Astrocytes[,1] %>% sum
# normCPM_res0.3_RNA_noScaledBrain$cluster_Astrocytes %>% head
# mean_CPM_perClusterTreatment$cluster_Astrocytes %>% head()

treatmentMeain_meanCPM_res0.3_RNA_noScaled_brain <- aggregate_treatment_means(
  cpm_mean_list = mean_CPM_perClusterTreatment,
  sample_map    = sample_map,
  sample_col    = "bc1_well",
  treatment_col = "treatment"
)


log2FC_normCPM_res0.3RNAnoScaledBrain <- compute_means_and_log2FC(
  matrix_list = normCPM_res0.3_RNA_noScaledBrain,
  sample_map = sample_map
)


edgerResults_res0.3RNAnoScaledBrain <- edger_pseudobulk_anova_per_cluster_v3(
  agg_list     = rawSumCounts_res0.3_RNA_noScaledBrain,
  sample_meta  = sample_map,
  group_levels = c("CTRL","DEX","CORT"),
  group_col    = "treatment",
  sample_col   = "bc1_well",
  verbose      = TRUE,
  progress     = TRUE,
  FDR_global_threshold = 1,
  gatekeep     = F,
  fixed_dispersion_value = 0.1
)


edgerResults_res0.3RNAnoScaledBrain$Astrocytes$DEX_vs_CTRL
edgerResults_res0.3RNAnoScaledBrain$`cluster_Radial glial cells`$DEX_vs_CTRL %>% 
  filter(logFC > 1, PValue < 0.01, logCPM > 4 )

seurat_featurePlot_binary_clusterSummary(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  gene_pattern = "GRM8",
  split.by = "treatment",
  split_levels = c("CTRL", "CORT", "DEX"),
  cluster_col = "res0.3_RNA_noScaled_brain",
  label_size = 12,
  point_size = 0.4,
  cluster_summary_figure = c("Astrocytes" = "Astrocytes", "Oligodendrocytes" = "Oligodendrocytes",
                             "Oligodendrocyte precursor cells" = "OPCs", "Radial glial cells" = "RGCs")
)


seurat_featurePlot_continuous_v3(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  gene_pattern = "LINGO2",
  split.by = "treatment",
  split_levels = c("CTRL", "CORT", "DEX"),
  cluster_col = "res0.3_RNA_noScaled_brain",
  point_size = 0.4,
  winsorize_top = 0.005,
  scale_zeroOne = TRUE
  # cluster_summary_figure = c("Astrocytes" = "Astrocytes", "Oligodendrocytes" = "Oligodendrocytes",
  #                            "Oligodendrocyte precursor cells" = "OPCs", "Radial glial cells" = "RGCs")
)

