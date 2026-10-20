


library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(tibble)

celltypes <- c(
  "Astrocytes", 
  "Radial glial cells", 
  "Oligodendrocyte precursor cells",
  "Oligodendrocytes", 
  "Immature neurons", 
  "Mature neurons"
)

# =========================================
# 1) Wygenerowanie vektora genów UNION
# =========================================

genes2heatmap <- combined_summaryEdgerChi2_res0.3RNAnoScaledBrain[celltypes] %>% 
  lapply(function(x){
    x %>%
      filter(meanCPM_CORT > 50 | meanCPM_DEX > 50) %>%
      filter(log2FC_edgerCORT > 1 | log2FC_edgerDEX > 1) %>%
      filter(PValue_edgerCORT < 0.01 | PValue_edgerDEX < 0.01) %>%
      .$gene_symbol
  }) %>% 
  unlist() %>% 
  unique()

length(genes2heatmap)   # kontrolnie ile jest genów

# =========================================
# 2) Wyciąganie QN z każdego klastra dla WYBRANYCH genów
# =========================================

df_list <- lapply(celltypes, function(cl) {
  df <- combined_summaryEdgerChi2_res0.3RNAnoScaledBrain[[cl]] %>%
    filter(gene_symbol %in% genes2heatmap) %>%
    select(gene_symbol, paste0("QN_A", 1:12))
  
  colnames(df)[-1] <- paste0(gsub(" ", "_", cl), "_", colnames(df)[-1])
  df
})

names(df_list) <- celltypes

# =========================================
# 3) Full join po gene_symbol
# =========================================

merged <- purrr::reduce(df_list, full_join, by = "gene_symbol") %>%
  column_to_rownames("gene_symbol")

mat <- as.matrix(merged)

# =========================================
# 4) Z-score per gene + clamp –2 do 2 + CLEAN NA
# =========================================

z_mat <- t(scale(t(mat)))
z_mat[!is.finite(z_mat)] <- 0           # usuń NaN / NA / Inf
z_mat <- pmax(pmin(z_mat, 2), -2)       # clamp
z_mat[!is.finite(z_mat)] <- 0

# =========================================
# 5) Kolejność kolumn CTRL → CORT → DEX
# =========================================

ctrl_cols <- grep("A([1-4])$", colnames(z_mat), value = TRUE)
cort_cols <- grep("A(5|6|7|8)$", colnames(z_mat), value = TRUE)
dex_cols  <- grep("A(9|1[0-2])$", colnames(z_mat), value = TRUE)

z_mat <- z_mat[, c(ctrl_cols, cort_cols, dex_cols)]

# annotation
groups <- ifelse(colnames(z_mat) %in% ctrl_cols, "CTRL",
                 ifelse(colnames(z_mat) %in% cort_cols, "CORT", "DEX"))

ha <- HeatmapAnnotation(
  df = data.frame(group = groups),
  col = list(group = c("CTRL"="grey40", "CORT"="firebrick2", "DEX"="steelblue3"))
)

# =========================================
# 6) HEATMAPA
# =========================================

Heatmap(
  z_mat,
  name = "z-score",
  top_annotation = ha,
  show_row_names = TRUE,
  show_column_names = TRUE,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_title = paste0("Genes (n=", nrow(z_mat), ")"),
  column_title = "QN Z-score across cell types",
  col = colorRamp2(c(-2,0,2), c("navy", "white", "firebrick2"))
)
