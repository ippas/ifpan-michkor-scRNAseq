seurat_featurePlot_binaryChi2(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  gene_pattern = "MERTK",
  split.by = "treatment",
  split_levels = c("CTRL", "CORT", "DEX"),  # zachowane, ale kolejność i tak wymuszona
  cluster_col = "res0.3_RNA_noScaled_brain",
  detection_threshold = 0,
  legend_position = "none",
  epsilon = 1e-6,
  reference_level = "CTRL",
  add_OR = TRUE,
  add_table_on_plot = F,
  font_size_stat_px = 12,
  font_size_table_px = 12,
  font_size_facet_title_px = 14,
  mode = "cell"
)

# GFAP, NR3C1

seurat_featurePlot_binary(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  gene_pattern = "OLIG1",
  split.by = "treatment",
  split_levels = c("CTRL", "CORT", "DEX"),  # zachowane, ale kolejność i tak wymuszona
  cluster_col = "res0.3_RNA_noScaled_brain",
  detection_threshold = 0,
  point_size = 0.4,
  return_data = T
) -> tmp


seurat_featurePlot_continuous(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  gene_pattern = "HIF3A",
  split.by = "treatment",
  split_levels = c("CTRL", "CORT", "DEX"),  # zachowane, ale kolejność i tak wymuszona
  cluster_col = "res0.3_RNA_noScaled_brain",
  # detection_threshold = 0,
  point_size = 0.4
)


edgeR_anovaOneFactor_treatment_pseudobulk_res0.4_clusters16_celltypes11_results $`cluster_Oligodendrocyte precursor cells`$f_test_global %>% filter(logFC.treatmentCORT  > 2)
