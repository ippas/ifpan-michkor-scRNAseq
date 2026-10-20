plot_multiGenes_clusters_multiGenesSvg(
  summary_list  = combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
  cluster_names = c("Radial glial cells","Astrocytes",
                    "Oligodendrocyte precursor cells","Oligodendrocytes",
                    "Immature neurons","Mature neurons"),
  genes_vector  = genes_global,
  output_dir    = "results/results_figures_26.11.2025_umapBarplotContinousBinary/globalUpGrSignatures_umapBarplot_CPM_percent/",
  width  = 16,
  height = 20,
  force_create_dir = TRUE,
  gene_title_size = 14,
  n_rows = 3,
  log2_transform = FALSE,
  border_size = 2
)

plot_multiGenes_clusters_multiGenesSvg(
  summary_list  = combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
  cluster_names = c("Radial glial cells","Astrocytes",
                    "Oligodendrocyte precursor cells","Oligodendrocytes",
                    "Immature neurons","Mature neurons"),
  genes_vector  = genes_neural,
  output_dir    = "results/results_figures_26.11.2025_umapBarplotContinousBinary/neuralUpGrSignatures_umapBarplot_CPM_percent/",
  width  = 16,
  height = 20,
  force_create_dir = TRUE,
  gene_title_size = 14,
  n_rows = 3,
  log2_transform = FALSE,
  border_size = 2
)

# ##############################################################################
# ---- DGE genes in RGS cells ----
# ##############################################################################
combined_summaryEdgerChi2_res0.3RNAnoScaledBrain$`Radial glial cells` %>% 
  filter(meanCPM_CORT > 50 | meanCPM_DEX > 50) %>% 
  filter(log2FC_edgerCORT > 1 | log2FC_edgerDEX > 1) %>%
  filter(PValue_edgerCORT < 0.01 | PValue_edgerDEX < 0.01) %>% .$gene_symbol -> genes_vector

plot_multiGenes_clusters_multiGenesSvg(
  summary_list  = combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
  cluster_names = c("Radial glial cells","Astrocytes",
                    "Oligodendrocyte precursor cells","Oligodendrocytes",
                    "Immature neurons","Mature neurons"),
  genes_vector  = genes_vector,
  output_dir    = "results/results_figures_26.11.2025_umapBarplotContinousBinary/RGS_CPM50p0.01log2FC1_umapBarplot_CPM_percent/",
  width  = 16,
  height = 20,
  force_create_dir = TRUE,
  gene_title_size = 14,
  n_rows = 3,
  log2_transform = FALSE,
  border_size = 2
)


# ##############################################################################
# ---- DGE genes in Astrocytes cells ----
# ##############################################################################
combined_summaryEdgerChi2_res0.3RNAnoScaledBrain$Astrocytes%>% 
  filter(meanCPM_CORT > 50 | meanCPM_DEX > 50) %>% 
  filter(log2FC_edgerCORT > 1 | log2FC_edgerDEX > 1) %>%
  filter(PValue_edgerCORT < 0.01 | PValue_edgerDEX < 0.01) %>% .$gene_symbol -> genes_vector

plot_multiGenes_clusters_multiGenesSvg(
  summary_list  = combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
  cluster_names = c("Radial glial cells","Astrocytes",
                    "Oligodendrocyte precursor cells","Oligodendrocytes",
                    "Immature neurons","Mature neurons"),
  genes_vector  = genes_vector,
  output_dir    = "results/results_figures_26.11.2025_umapBarplotContinousBinary/Astrocytes_CPM50p0.01log2FC1_umapBarplot_CPM_percent/",
  width  = 16,
  height = 20,
  force_create_dir = TRUE,
  gene_title_size = 14,
  n_rows = 3,
  log2_transform = FALSE,
  border_size = 2
)


# ##############################################################################
# ---- DGE genes in OPC cells ----
# ##############################################################################
combined_summaryEdgerChi2_res0.3RNAnoScaledBrain$`Oligodendrocyte precursor cells` %>% 
  filter(meanCPM_CORT > 50 | meanCPM_DEX > 50) %>% 
  filter(log2FC_edgerCORT > 1 | log2FC_edgerDEX > 1) %>%
  filter(PValue_edgerCORT < 0.01 | PValue_edgerDEX < 0.01) %>% .$gene_symbol -> genes_vector

plot_multiGenes_clusters_multiGenesSvg(
  summary_list  = combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
  cluster_names = c("Radial glial cells","Astrocytes",
                    "Oligodendrocyte precursor cells","Oligodendrocytes",
                    "Immature neurons","Mature neurons"),
  genes_vector  = genes_vector,
  output_dir    = "results/results_figures_26.11.2025_umapBarplotContinousBinary/OPC_CPM50p0.01log2FC1_umapBarplot_CPM_percent/",
  width  = 16,
  height = 20,
  force_create_dir = TRUE,
  gene_title_size = 14,
  n_rows = 3,
  log2_transform = FALSE,
  border_size = 2
)


# ##############################################################################
# ---- DGE genes in Oligodendrocytes cells ----
# ##############################################################################
combined_summaryEdgerChi2_res0.3RNAnoScaledBrain$Oligodendrocytes %>% 
  filter(meanCPM_CORT > 50 | meanCPM_DEX > 50) %>% 
  filter(log2FC_edgerCORT > 1 | log2FC_edgerDEX > 1) %>%
  filter(PValue_edgerCORT < 0.01 | PValue_edgerDEX < 0.01) %>% .$gene_symbol -> genes_vector

plot_multiGenes_clusters_multiGenesSvg(
  summary_list  = combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
  cluster_names = c("Radial glial cells","Astrocytes",
                    "Oligodendrocyte precursor cells","Oligodendrocytes",
                    "Immature neurons","Mature neurons"),
  genes_vector  = genes_vector,
  output_dir    = "results/results_figures_26.11.2025_umapBarplotContinousBinary/Oligodendrocytes_CPM50p0.01log2FC1_umapBarplot_CPM_percent/",
  width  = 16,
  height = 20,
  force_create_dir = TRUE,
  gene_title_size = 14,
  n_rows = 3,
  log2_transform = FALSE,
  border_size = 2
)


# ##############################################################################
# ---- DGE genes in Immature Neurons cells ----
# ##############################################################################
combined_summaryEdgerChi2_res0.3RNAnoScaledBrain$`Immature neurons` %>% 
  filter(meanCPM_CORT > 50 | meanCPM_DEX > 50) %>% 
  filter(log2FC_edgerCORT > 1 | log2FC_edgerDEX > 1) %>%
  filter(PValue_edgerCORT < 0.01 | PValue_edgerDEX < 0.01) %>% .$gene_symbol -> genes_vector

plot_multiGenes_clusters_multiGenesSvg(
  summary_list  = combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
  cluster_names = c("Radial glial cells","Astrocytes",
                    "Oligodendrocyte precursor cells","Oligodendrocytes",
                    "Immature neurons","Mature neurons"),
  genes_vector  = genes_vector,
  output_dir    = "results/results_figures_26.11.2025_umapBarplotContinousBinary/ImmatureNeurons_CPM50p0.01log2FC1_umapBarplot_CPM_percent/",
  width  = 16,
  height = 20,
  force_create_dir = TRUE,
  gene_title_size = 14,
  n_rows = 3,
  log2_transform = FALSE,
  border_size = 2
)


# ##############################################################################
# ---- DGE genes in Mature Neurons cells ----
# ##############################################################################
combined_summaryEdgerChi2_res0.3RNAnoScaledBrain$`Mature neurons` %>% 
  filter(meanCPM_CORT > 50 | meanCPM_DEX > 50) %>% 
  filter(log2FC_edgerCORT > 1 | log2FC_edgerDEX > 1) %>%
  filter(PValue_edgerCORT < 0.01 | PValue_edgerDEX < 0.01) %>% .$gene_symbol -> genes_vector

plot_multiGenes_clusters_multiGenesSvg(
  summary_list  = combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
  cluster_names = c("Radial glial cells","Astrocytes",
                    "Oligodendrocyte precursor cells","Oligodendrocytes",
                    "Immature neurons","Mature neurons"),
  genes_vector  = genes_vector,
  output_dir    = "results/results_figures_26.11.2025_umapBarplotContinousBinary/MatureNeurons_CPM50p0.01log2FC1_umapBarplot_CPM_percent/",
  width  = 16,
  height = 20,
  force_create_dir = TRUE,
  gene_title_size = 14,
  n_rows = 3,
  log2_transform = FALSE,
  border_size = 2
)


