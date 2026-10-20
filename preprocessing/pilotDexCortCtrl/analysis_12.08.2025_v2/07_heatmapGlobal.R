pseudobulk_heatmap_global_complexHeatmap(
  agg = agg_list,
  sample_meta = sample_map,
  sample_col = "sample_id",
  group_col  = "treatment",
  group_levels = c("CTRL","DEX","CORT"),
  transform  = "logCPM",   # "logCPM" (edgeR), "log2", "none"
  # top_n_var  = 100,
  scale_rows = TRUE,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  legend_title = "Expression (z)",
  genes = res_global$f_test_global %>% 
    filter(FDR < 0.05) %>%
    filter(logFC.treatmentCORT > 1 | logFC.treatmentCORT < -1 | logFC.treatmentDEX > 1 | logFC.treatmentDEX < -1) %>%
    head(100) %>%
    .$name_id,
  quantile_normalize = TRUE,
  z_score_limits = c(-2, 2),
  heatmap_colors = c("#6986acff", "white", "#802226ff")
)


# col_low = "#6986acff",
# col_mid = "white",
# col_high = "#802226ff"