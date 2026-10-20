# ============================================
# 🔹 HEATMAPA CramerV_signed dla CORT i DEX
#    (tylko geny z przynajmniej jednym FDR < 0.05)
# ============================================

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ComplexHeatmap)
library(circlize)

# --------------------------------------------
# 1️⃣ Połączenie wszystkich wyników w jedną ramkę danych
# --------------------------------------------

cell_order <- c("Astrocytes", "Oligodendrocytes", "Mature neurons")

heatmap_df <- results_list %>%
  purrr::imap_dfr(~ mutate(.x, cluster = .y)) %>%
  select(cluster, gene_symbol, treatment, CramerV_signed, p_chi2) %>%
  filter(cluster %in% cell_order)

# --------------------------------------------
# 2️⃣ Obliczanie FDR w obrębie każdej kombinacji
# --------------------------------------------

heatmap_df <- heatmap_df %>%
  group_by(cluster, treatment) %>%
  mutate(FDR = p.adjust(p_chi2, method = "fdr")) %>%
  ungroup()

# --------------------------------------------
# 3️⃣ Czyszczenie danych i filtrowanie genów istotnych po FDR
# --------------------------------------------

# Zamień NA na 0 i wyzeruj nieistotne efekty
heatmap_df <- heatmap_df %>%
  mutate(
    CramerV_signed = ifelse(is.na(CramerV_signed), 0, CramerV_signed),
    CramerV_signed = ifelse(FDR < 0.05, CramerV_signed, 0)
  )

# Zostaw tylko geny z przynajmniej jednym FDR < 0.05
genes_keep <- heatmap_df %>%
  filter(FDR < 0.05) %>%
  distinct(gene_symbol) %>%
  pull(gene_symbol)

heatmap_df <- heatmap_df %>% filter(gene_symbol %in% genes_keep)

# --------------------------------------------
# 4️⃣ Przekształcenie do formatu szerokiego
# --------------------------------------------

heatmap_df_wide <- heatmap_df %>%
  unite(col = "col_id", treatment, cluster, sep = "_") %>%
  select(gene_symbol, col_id, CramerV_signed) %>%
  pivot_wider(names_from = col_id, values_from = CramerV_signed) %>%
  tibble::column_to_rownames("gene_symbol") %>%
  as.matrix()

# --------------------------------------------
# 5️⃣ Wymuszenie kolejności kolumn
# --------------------------------------------

desired_col_order <- c(
  paste0("CORT_", cell_order),
  paste0("DEX_", cell_order)
)

desired_col_order <- desired_col_order[desired_col_order %in% colnames(heatmap_df_wide)]
heatmap_df_wide <- heatmap_df_wide[, desired_col_order, drop = FALSE]

# --------------------------------------------
# 6️⃣ Przygotowanie adnotacji kolumn
# --------------------------------------------

treatment_labels <- ifelse(grepl("^CORT_", colnames(heatmap_df_wide)), "CORT", "DEX")
treatment_colors <- c("CORT" = "#3E96FF", "DEX" = "#FF4E4E")

cluster_labels <- gsub("^(CORT_|DEX_)", "", colnames(heatmap_df_wide))
cluster_colors <- c(
  "Astrocytes" = "#66c2a5",
  "Oligodendrocytes" = "#fc8d62",
  "Mature neurons" = "#8da0cb"
)

top_annot <- HeatmapAnnotation(
  Treatment = treatment_labels,
  Cluster = cluster_labels,
  col = list(
    Treatment = treatment_colors,
    Cluster = cluster_colors
  ),
  annotation_name_gp = gpar(fontsize = 9),
  simple_anno_size = unit(4, "mm")
)

# --------------------------------------------
# 7️⃣ Definicja kolorów dla wartości
# --------------------------------------------

col_fun <- colorRamp2(c(-0.5, 0, 0.5),
                      c("dodgerblue3", "white", "firebrick3"))

# --------------------------------------------
# 8️⃣ Rysowanie heatmapy
# --------------------------------------------

Heatmap(
  heatmap_df_wide,
  name = "Cramér’s V (signed)",
  col = col_fun,
  top_annotation = top_annot,
  cluster_rows = TRUE,
  cluster_columns = FALSE,  # zachowaj kolejność CORT→DEX
  row_names_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 8),
  column_names_rot = 45,
  border = TRUE,
  heatmap_legend_param = list(
    title = "Cramér’s V",
    at = c(-0.5, 0, 0.5),
    labels = c("↓", "0", "↑")
  )
)
