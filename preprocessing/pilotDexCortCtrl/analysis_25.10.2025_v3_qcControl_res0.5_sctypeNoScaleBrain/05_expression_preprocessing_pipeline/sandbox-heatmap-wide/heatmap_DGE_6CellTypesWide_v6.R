library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(tibble)
library(purrr)
library(grid)

celltypes <- c(
  "Astrocytes",
  "Radial glial cells",
  "Oligodendrocyte precursor cells",
  "Oligodendrocytes",
  "Immature neurons",
  "Mature neurons"
)

# =========================================
# 1) UNION genów spełniających kryteria
# =========================================
genes2heatmap <- combined_summaryEdgerChi2_res0.3RNAnoScaledBrain[celltypes] %>% 
  lapply(function(x){
    x %>%
      filter(meanCPM_CORT > 50 | meanCPM_DEX > 50) %>%
      filter(log2FC_edgerCORT > 1 | log2FC_edgerDEX > 1) %>%
      filter(PValue_edgerCORT < 0.01 | PValue_edgerDEX < 0.01) %>%
      .$gene_symbol
  }) %>% 
  unlist() %>% unique()

cat("Liczba genów do heatmapy:", length(genes2heatmap), "\n")


# =========================================
# 2) Normalizacja na CTRL + Z-score względem CTRL
# =========================================
df_list <- lapply(celltypes, function(cl) {
  
  df_raw <- combined_summaryEdgerChi2_res0.3RNAnoScaledBrain[[cl]] %>%
    filter(gene_symbol %in% genes2heatmap) %>%
    select(gene_symbol, paste0("QN_A", 1:12))
  
  # średnia CTRL A1-A4
  ctrl_mean <- df_raw %>%
    select(QN_A1, QN_A2, QN_A3, QN_A4) %>%
    rowMeans()
  
  # normalizacja do mean CTRL (CTRL powinno być ≈1)
  df_norm <- df_raw
  df_norm[, -1] <- sweep(df_norm[, -1], 1, ctrl_mean, "/")
  
  # odchylenie std CTRL do z-score
  ctrl_sd <- df_norm %>%
    select(QN_A1, QN_A2, QN_A3, QN_A4) %>%
    apply(1, sd)
  
  # Z-score względem CTRL: (X-1)/SD_CTRL
  zscores <- as.matrix(df_norm[, -1])
  zscores <- sweep(zscores, 1, 1, "-")
  zscores <- sweep(zscores, 1, ctrl_sd, "/")
  zscores[!is.finite(zscores)] <- 0
  
  df_z <- data.frame(
    gene_symbol = df_norm$gene_symbol,
    as.data.frame(zscores)
  )
  
  colnames(df_z)[-1] <- paste0(gsub(" ", "_", cl), "_", colnames(df_z)[-1])
  
  df_z
})

names(df_list) <- celltypes


# =========================================
# 3) Full join & matrix
# =========================================
merged <- purrr::reduce(df_list, full_join, by = "gene_symbol") %>%
  column_to_rownames("gene_symbol")

mat <- as.matrix(merged)
mat[!is.finite(mat)] <- 0


# =========================================
# 4) Kolejność kolumn CTRL → CORT → DEX
# =========================================
ctrl_cols <- grep("A([1-4])$", colnames(mat), value = TRUE)
cort_cols <- grep("A(5|6|7|8)$", colnames(mat), value = TRUE)
dex_cols  <- grep("A(9|1[0-2])$", colnames(mat), value = TRUE)

mat <- mat[, c(ctrl_cols, cort_cols, dex_cols)]

groups <- ifelse(colnames(mat) %in% ctrl_cols, "CTRL",
                 ifelse(colnames(mat) %in% cort_cols, "CORT", "DEX"))

ha <- HeatmapAnnotation(
  df = data.frame(group = groups),
  col = list(group = c("CTRL"="grey40", "CORT"="firebrick2", "DEX"="steelblue3")),
  annotation_name_gp = gpar(fontsize = 16, fontface = "bold")
)


# =========================================
# 5) HEATMAPA
# =========================================
Heatmap(
  mat,
  name = "Z-score\nvs CTRL",
  top_annotation = ha,
  show_row_names = TRUE,
  show_column_names = TRUE,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_title = paste0("Genes (n=", nrow(mat), ")"),
  row_title_gp = gpar(fontsize = 18, fontface="bold"),
  column_title = "Normalized QN vs CTRL + Z-score relative to CTRL",
  column_title_gp = gpar(fontsize = 20, fontface="bold"),
  col = colorRamp2(c(-10, 0, 10), c("navy", "white", "firebrick4")),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 16, fontface="bold"),
    labels_gp = gpar(fontsize = 14)
  )
)
