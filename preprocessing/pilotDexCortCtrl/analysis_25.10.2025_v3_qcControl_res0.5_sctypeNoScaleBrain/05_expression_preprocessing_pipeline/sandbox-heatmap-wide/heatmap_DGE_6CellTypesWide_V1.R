library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(tibble)

# ===============================================
# 0) Klastery do użycia
# ===============================================

clusters_to_use <- c(
  "Astrocytes", "Radial glial cells", "Oligodendrocyte precursor cells",
  "Oligodendrocytes", "Immature neurons", "Mature neurons"
)

# filtr genów jak wcześniej
filter_genes <- function(df) {
  df %>%
    dplyr::filter(meanCPM_CORT > 50 | meanCPM_DEX > 50) %>%
    dplyr::filter(log2FC_edgerCORT > 1 | log2FC_edgerDEX > 1) %>%
    dplyr::filter(PValue_edgerCORT < 0.01 | PValue_edgerDEX < 0.01)
}

# ===============================================
# 1) Filtrowanie w każdym klastrze + QN_A1–QN_A12
# ===============================================

df_list <- lapply(clusters_to_use, function(cl_name) {
  tmp <- combined_summaryEdgerChi2_res0.3RNAnoScaledBrain[[cl_name]] %>%
    filter_genes() %>%
    dplyr::select(gene_symbol, paste0("QN_A", 1:12))
  
  # prefix z nazwą klastra w kolumnach
  colnames(tmp)[-1] <- paste0(gsub(" ", "_", cl_name), "_", colnames(tmp)[-1])
  tmp
})

names(df_list) <- clusters_to_use

# ===============================================
# 2) Full join po gene_symbol
# ===============================================

merged <- purrr::reduce(df_list, dplyr::full_join, by = "gene_symbol") %>%
  tibble::column_to_rownames("gene_symbol")

mat <- as.matrix(merged)

# ===============================================
# 3) Z-score po wierszach, czyszczenie NA/NaN/Inf, clamp [-2, 2]
# ===============================================

z_mat <- t(scale(t(mat)))          # z-score w wierszach

# kluczowa linia: usuwamy NA, NaN i Inf (wszystko, co nie jest skończone)
z_mat[!is.finite(z_mat)] <- 0      # NA, NaN, Inf -> 0 (średnia w z-score)

# dopiero teraz przycinamy zakres
z_mat <- pmax(pmin(z_mat, 2), -2)

# na wszelki wypadek, jeśli jakiś dziad się ostał:
# (nie powinno być potrzebne, ale nie zaszkodzi)
z_mat[!is.finite(z_mat)] <- 0

# ===============================================
# 4) Kolejność próbek: CTRL (A1–A4), CORT (A5–A8), DEX (A9–A12)
# ===============================================

ctrl_cols <- grep("A([1-4])$", colnames(z_mat), value = TRUE)
cort_cols <- grep("A(5|6|7|8)$", colnames(z_mat), value = TRUE)
dex_cols  <- grep("A(9|1[0-2])$", colnames(z_mat), value = TRUE)

z_mat <- z_mat[, c(ctrl_cols, cort_cols, dex_cols)]

# ===============================================
# 5) Annotation: CTRL / CORT / DEX na górze
# ===============================================

groups <- ifelse(colnames(z_mat) %in% ctrl_cols, "CTRL",
                 ifelse(colnames(z_mat) %in% cort_cols, "CORT", "DEX"))

condition <- data.frame(group = groups)
rownames(condition) <- colnames(z_mat)

ha <- ComplexHeatmap::HeatmapAnnotation(
  df = condition,
  col = list(group = c("CTRL"="grey40", "CORT"="firebrick2", "DEX"="steelblue3")),
  annotation_name_side = "left"
)

# ===============================================
# 6) HEATMAPA
# ===============================================

ComplexHeatmap::Heatmap(
  z_mat,
  name = "z-score",
  top_annotation = ha,
  show_row_names = TRUE,
  show_column_names = TRUE,
  cluster_rows = TRUE,        # klastrujemy geny
  cluster_columns = FALSE,    # kolumny bez klastrowania
  row_title = paste0("Genes (n=", nrow(z_mat), ")"),
  column_title = "Combined QN Z-score across 6 neural cell types",
  col = circlize::colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick2"))
)
