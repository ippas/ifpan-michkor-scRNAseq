binaryGeneDetectionTreatment_res0.3RNAnoScaledBrain <- seurat_binaryDetection_summaryByCluster_fast(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  clusters_to_check = NULL,  
  genes_vector = hgnc_symbols_vector_v110,                # NEW
  assay = "RNA",
  slot = "counts",
  split.by = "treatment",
  split_levels = c("CTRL", "CORT", "DEX"),
  cluster_col = "res0.3_RNA_noScaled_brain",
  detection_threshold = 1,            # counts >= 1
  verbose = TRUE
)

binaryGeneDetectionTreatment_res0.3RNAnoScaledBrain$Astrocytes %>% 
  filter(gene_symbol %in% c("FKBP5", "TSC22D3", "STON1-GTF2A1L"))

chi2Results_res0.3RNAnoScaledBrain <- analyze_detection_chi2_byCluster(
  detection_list = binaryGeneDetectionTreatment_res0.3RNAnoScaledBrain,
  return_long = FALSE
)


df_sample <- aggregate_detection_by_cluster_treatment_sample(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  cluster_col = "res0.3_RNA_noScaled_brain",
  treatment_col = "treatment",
  sample_col = "bc1_well",
  
  # clusters_to_check = NULL,  
  genes_keep = hgnc_symbols_vector_v110,                # NEW
  assay = "RNA",
  slot = "counts",
  # split_levels = c("CTRL", "CORT", "DEX"),
  detection_threshold = 1
)
