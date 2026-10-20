# sample_map
# 
# combined_summaryEdgerChi2_res0.3RNAnoScaledBrain$`Radial glial cells` %>% 
#   filter(meanCPM_CORT > 50 | meanCPM_DEX > 50) %>% 
#   filter(log2FC_edgerCORT > 1 | log2FC_edgerDEX > 1) %>%
#   filter(PValue_edgerCORT < 0.01 | PValue_edgerDEX < 0.01) %>% .$gene_symbol
#   # dplyr::select(gene_symbol, meanCPM_CORT, meanCPM_DEX, log2FC_edgerCORT, log2VC_edgerDEX,
#   #               PValue_edgerCORT, PValue_edgerDEX)
#   # 


  
  
# ##############################################################################
# ---- test heatmap ----
# ##############################################################################
# ---- 1) Filtrowanie ------------------------------------------
filtered_df <- combined_summaryEdgerChi2_res0.3RNAnoScaledBrain$`Mature neurons` %>%
  filter(meanCPM_CORT > 50 | meanCPM_DEX > 50) %>%
  filter(log2FC_edgerCORT > 1 | log2FC_edgerDEX > 1) %>%
  filter(PValue_edgerCORT < 0.01 | PValue_edgerDEX < 0.01)

# ---- 2) Wybór kolumn QN_A1:QN_A12 + gene_symbol ---------------
qn_cols <- paste0("QN_A", 1:12)
df_qn <- filtered_df %>%
  select(gene_symbol, all_of(qn_cols))

# ---- 3) Matryca danych ----------------------------------------
mat <- df_qn %>%
  column_to_rownames("gene_symbol") %>%
  as.matrix()

# ---- 4) Z-score per gene (wiersz) i ograniczenie do [-2, 2] ---
z_mat <- t(scale(t(mat)))     # z-score by rows
z_mat <- pmax(pmin(z_mat, 2), -2)  # clamp to [-2,2]

# ---- 5) Kolejność próbek (CTRL → CORT → DEX) ------------------
new_order <- c(paste0("QN_A", 1:4),
               paste0("QN_A", 5:8),
               paste0("QN_A", 9:12))
z_mat <- z_mat[, new_order]

condition <- data.frame(
  group = c(
    rep("CTRL", 4),
    rep("DEX", 4),
    rep("CORT", 4)
  )
)
rownames(condition) <- colnames(z_mat)

ha_col <- list(
  group = c("CTRL" = "grey40", "CORT" = "firebrick2", "DEX" = "steelblue3")
)

ha <- HeatmapAnnotation(
  df = condition,
  col = ha_col,
  annotation_name_side = "left"
)

# ---- 6) HEATMAPA ----------------------------------------------
Heatmap(
  z_mat,
  name = "z-score",
  top_annotation = ha,
  show_row_names = TRUE,
  show_column_names = TRUE,
  cluster_rows = TRUE,
  cluster_columns = FALSE,   # bez klastrowania kolumn
  column_title = paste0("QN Z-score (n=", nrow(z_mat), ")"),
  row_title = "Genes",
  heatmap_legend_param = list(title = "z", legend_height = unit(4, "cm")),
  col = colorRamp2(
    c(-2, 0, 2),
    c("navy", "white", "darkred")
  )
)

