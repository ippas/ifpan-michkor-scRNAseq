clusterAll_summary_list_neuralGlobal$Astrocytes %>% 
  filter(DEX_percent  < 5) %>% .$gene_symbol
  filter(gene_symbol == "FKBP5")


seurat_featurePlot_binary_clusterSummary(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  gene_pattern = "PER1",
  split.by = "treatment",
  split_levels = c("CTRL", "CORT", "DEX"),
  cluster_col = "res0.3_RNA_noScaled_brain",
  label_size = 12,
  point_size = 0.4,
  cluster_summary_figure = c("Astrocytes" = "Astrocytes", "Oligodendrocytes" = "Oligodendrocytes",
                             "Oligodendrocyte precursor cells" = "OPCs", "Radial glial cells" = "RGCs")
)

# PDK4, TSC22D3
results_list$Astrocytes %>% 
  filter(p_chi2 < 0.05) %>% 
  filter(treatment == "DEX") %>% .$gene_symbol %>% unique
