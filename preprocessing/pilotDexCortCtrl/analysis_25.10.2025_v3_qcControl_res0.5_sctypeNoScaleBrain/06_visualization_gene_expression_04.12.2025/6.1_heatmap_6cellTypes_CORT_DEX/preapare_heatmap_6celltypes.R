# ##############################################################################
# ---- prepare heatmap and data ----
# ##############################################################################
res_heatmap <- generate_heatmap_ctrlNormalized_zscore_per_block(
  data_list        = combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
  celltypes        = c(
    "Radial glial cells",
    "Astrocytes",
    "OPCs",
    "Oligodendrocytes",
    "Immature neurons",
    "Mature neurons"
  ),
  distance_method  = "spearman",
  linkage_method   = "complete",
  logfc_mode       = "both",
  contrast         = "both",
  statistics_mode  = "pvalue",
  pval_threshold   = 1e-5,
  logfc_threshold  = 1,
  zscore_range     = c(-2,0,2),
  cluster_rows_k   = 5,
  show_rows        = TRUE,
  show_cols        = FALSE,
  row_dend_width   = unit(20,"mm"), 
  cpm_filter_mode  = "max",
  seed = 9,
  # svg_filename = "results/figures/heatmap_DGE_6celltypes_08.12.2025/heatmapDGE_268genesPvalue10e5CPM50log2FC1_09.12.2025_v1.svg",
  width = 11,
  height = 6
)



# fdr
res_heatmap2 <- generate_heatmap_ctrlNormalized_zscore_per_block(
  data_list        = combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
  celltypes        = c(
    # "Radial glial cells",
    "Astrocytes",
    # "OPCs",
    "Oligodendrocytes",
    "Immature neurons"
    # "Mature neurons"
  ),
  distance_method  = "spearman",
  linkage_method   = "complete",
  logfc_mode       = "both",
  contrast         = "both",
  statistics_mode  = "fdr",
  pval_threshold   = 1e-3,
  logfc_threshold  = 0.5,
  zscore_range     = c(-2,0,2),
  cluster_rows_k   = 4,
  show_rows        = TRUE,
  show_cols        = FALSE,
  row_dend_width   = unit(20,"mm"), 
  cpm_filter_mode  = "max",
  seed = 10,
  # svg_filename = "results/figures/heatmap_DGE_6celltypes_08.12.2025/heatmapDGE_268genesPvalue10e5CPM50log2FC1_09.12.2025_v1.svg",
  width = 11,
  height = 6
)



generate_heatmap_ctrlNormalized_zscore_per_block(
  data_list        = combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
  celltypes        = c(
    "Radial glial cells",
    "Astrocytes",
    "OPCs",
    "Oligodendrocytes",
    "Immature neurons",
    "Mature neurons"
  ),
  distance_method  = "spearman",
  linkage_method   = "complete",
  logfc_mode       = "both",
  contrast         = "both",
  statistics_mode  = "pvalue",
  pval_threshold   = 1e-5,
  logfc_threshold  = 1,
  zscore_range     = c(-2,0,2),
  cluster_rows_k   = 5,
  show_rows        = TRUE,
  show_cols        = FALSE,
  row_dend_width   = unit(20,"mm"), 
  cpm_filter_mode  = "max",
  seed = 9,
  svg_filename = "results/figures/heatmap_DGE_6celltypes_11.12.2025/heatmapDGE_268genesPvalue10e5CPM50log2FC1_11.12.2025_v2_long.svg",
  width = 11,
  height = 30,
  fontsize_row_names = 7
)



res_heatmap$genes_per_cluster %>% lapply(length)

intersect(
  res_heatmap$genes_per_cluster$cluster_1,
  genes_global
)


res_heatmap$cluster_genes$Cluster_2 %>% cat(sep = "\n")


combined_summaryEdgerChi2_res0.3RNAnoScaledBrain %>% 
  lapply(., function(x){
    x %>% filter(gene_symbol %in% res_heatmap$cluster_genes$Cluster_1)
  }) %>% 
  .[celltypes] %>% 
  bind_rows(.id = "celltype") %>% 
  filter(gene_symbol == "CARMIL1")


c(
  genes_global,
  genes_GlobalDown,
  genes_neural,
  genes_neuralDown
) %>% unique %>% length()


intersect(
  signif268genes,
  c(
    genes_global,
    genes_GlobalDown,
    genes_neural,
    genes_neuralDown
  )
) %>% length()

combined_summaryEdgerChi2_res0.3RNAnoScaledBrain[
  c(
    "Astrocytes",
    "Immature neurons", 
    "Mature neurons", 
    "OPCs", 
    "Oligodendrocytes",
    "Radial glial cells"
  )
] %>% 
  lapply(., function(x){
    x %>% filter(max_CPM > 50) %>% 
      filter(PValue_edgerDEX < 0.00001 | PValue_edgerCORT < 0.00001) %>% .$gene_symbol
      # filter(abs(log2FC_edgerDEX) > 1 | abs(log2FC_edgerCORT) > 1) %>% .$gene_symbol
  }) %>% unname %>% unlist %>% unique -> signif268genes

combined_summaryEdgerChi2_res0.3RNAnoScaledBrain$Astrocytes$gene_symbol %>% unique %>% length()
