seurat_featurePlot_binary_clusterSummary(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  gene_pattern = "TSC22D1",
  split.by = "treatment",
  split_levels = c("CTRL", "CORT", "DEX"),
  cluster_col = "res0.3_RNA_noScaled_brain",
  label_size = 12,
  point_size = 0.4
) -> p1


seurat_featurePlot_continuous_v3(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  gene_pattern = "TSC22D1",
  split.by = "treatment",
  split_levels = c("CTRL", "CORT", "DEX"),
  cluster_col = "res0.3_RNA_noScaled_brain",
  point_size = 0.4,
  winsorize_top = 0.01,
  scale_zeroOne = FALSE
) -> p2

p1 / p2


combined_summaryEdgerChi2_res0.3RNAnoScaledBrain$Astrocytes %>% 
  filter(gene_symbol == "FKBP5")

gene_symbol <- "FKBP5"

# ---- p1 ----
barplot_meanCPM_gene(
  summary_list = combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
  cluster_name = "Radial glial cells",
  gene_symbol = gene_symbol,
  bar_color = "black",
  point_size = 5,
  bar_width = 0.7,
  bar_size = 1,
  log2_transform = TRUE,
  add_points = T,
  se_width = 0.1,
  add_log2FC_CPM = TRUE,
  show_p = TRUE,
  text_size_values = 0,  # size in pt
  text_size_p = 12,       # size in pt
  text_size_title = 10,   # size in pt
  text_size_axis = 12,     # size in pt
  custom_title = "",
  stats_bracket_mult = 0.75,   # shift of bracket height (relative multiplier)
  stats_text_mult    = 0.35    # shift of text height above bracket
)

barplot_percentDetection_gene(
  summary_list = combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
  cluster_name = "Radial glial cells",
  gene_symbol = gene_symbol,
  bar_color = "black",
  bar_width = 0.7,
  show_p = TRUE,
  text_size_values = 10,
  text_size_p = 12,
  text_size_title = 10,
  text_size_axis = 12,
  custom_title = ""
)

p1_CPM + p1_percent -> p1

# ---- p2 ----

barplot_meanCPM_gene(
  summary_list = combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
  cluster_name = "Astrocytes",
  gene_symbol = gene_symbol,
  bar_color = "black",
  point_size = 5,
  bar_width = 0.7,
  bar_size = 1,
  log2_transform = TRUE,
  add_points = T,
  se_width = 0.1,
  show_p = TRUE,
  text_size_values = 0,  # size in pt
  text_size_p = 12,       # size in pt
  text_size_title = 10,   # size in pt
  text_size_axis = 12     # size in pt
) -> p2_CPM

barplot_percentDetection_gene(
  summary_list = combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
  cluster_name = "Astrocytes",
  gene_symbol = gene_symbol,
  bar_color = "black",
  bar_width = 0.7,
  show_p = TRUE,
  # y_limits = c(0, 60),
  text_size_values = 10,
  text_size_p = 12,
  text_size_title = 10,
  text_size_axis = 12
) -> p2_percent

p2_CPM + p2_percent -> p2

# ---- p3 ----

barplot_meanCPM_gene(
  summary_list = combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
  cluster_name = "Oligodendrocyte precursor cells",
  gene_symbol = gene_symbol,
  bar_color = "black",
  point_size = 5,
  bar_width = 0.7,
  bar_size = 1,
  log2_transform = TRUE,
  add_points = T,
  se_width = 0.1,
  show_p = TRUE,
  text_size_values = 0,  # size in pt
  text_size_p = 12,       # size in pt
  text_size_title = 10,   # size in pt
  text_size_axis = 12     # size in pt
) -> p3_CPM

barplot_percentDetection_gene(
  summary_list = combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
  cluster_name = "Oligodendrocyte precursor cells",
  gene_symbol = gene_symbol,
  bar_color = "black",
  bar_width = 0.7,
  show_p = TRUE,
  # y_limits = c(0, 60),
  text_size_values = 10,
  text_size_p = 12,
  text_size_title = 10,
  text_size_axis = 12
) -> p3_percent

p3_CPM + p3_percent -> p3

# ---- p4 ----

barplot_meanCPM_gene(
  summary_list = combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
  cluster_name = "Oligodendrocytes",
  gene_symbol = gene_symbol,
  bar_color = "black",
  point_size = 5,
  bar_width = 0.7,
  bar_size = 1,
  log2_transform = TRUE,
  add_points = T,
  se_width = 0.1,
  show_p = TRUE,
  text_size_values = 0,  # size in pt
  text_size_p = 12,       # size in pt
  text_size_title = 10,   # size in pt
  text_size_axis = 12     # size in pt
) -> p4_CPM

barplot_percentDetection_gene(
  summary_list = combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
  cluster_name = "Oligodendrocytes",
  gene_symbol = gene_symbol,
  bar_color = "black",
  bar_width = 0.7,
  show_p = TRUE,
  # y_limits = c(0, 60),
  text_size_values = 10,
  text_size_p = 12,
  text_size_title = 10,
  text_size_axis = 12
) -> p4_percent

p4_CPM + p4_percent -> p4

# ---- p5 ----

barplot_meanCPM_gene(
  summary_list = combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
  cluster_name = "Immature neurons",
  gene_symbol = gene_symbol,
  bar_color = "black",
  point_size = 5,
  bar_width = 0.7,
  bar_size = 1,
  log2_transform = TRUE,
  add_points = T,
  se_width = 0.1,
  show_p = TRUE,
  text_size_values = 0,  # size in pt
  text_size_p = 12,       # size in pt
  text_size_title = 10,   # size in pt
  text_size_axis = 12     # size in pt
) -> p5_CPM

barplot_percentDetection_gene(
  summary_list = combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
  cluster_name = "Immature neurons",
  gene_symbol = gene_symbol,
  bar_color = "black",
  bar_width = 0.7,
  show_p = TRUE,
  # y_limits = c(0, 60),
  text_size_values = 10,
  text_size_p = 12,
  text_size_title = 10,
  text_size_axis = 12
) -> p5_percent

p5_CPM + p5_percent -> p5


# ---- p6 ----

barplot_meanCPM_gene(
  summary_list = combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
  cluster_name = "Mature neurons",
  gene_symbol = gene_symbol,
  bar_color = "black",
  point_size = 5,
  bar_width = 0.7,
  bar_size = 1,
  log2_transform = TRUE,
  add_points = T,
  se_width = 0.1,
  show_p = TRUE,
  text_size_values = 0,  # size in pt
  text_size_p = 12,       # size in pt
  text_size_title = 10,   # size in pt
  text_size_axis = 12     # size in pt
) -> p6_CPM

barplot_percentDetection_gene(
  summary_list = combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
  cluster_name = "Mature neurons",
  gene_symbol = gene_symbol,
  bar_color = "black",
  bar_width = 0.7,
  show_p = TRUE,
  # y_limits = c(0, 60),
  text_size_values = 10,
  text_size_p = 12,
  text_size_title = 10,
  text_size_axis = 12
) -> p6_percent

p6_CPM + p6_percent -> p6



(
  p1_CPM    | p1_percent | p2_CPM    | p2_percent
) /
(
  p3_CPM    | p3_percent | p4_CPM    | p4_percent
)  /
(
  p5_CPM    | p5_percent | p6_CPM    | p6_percent
) 


