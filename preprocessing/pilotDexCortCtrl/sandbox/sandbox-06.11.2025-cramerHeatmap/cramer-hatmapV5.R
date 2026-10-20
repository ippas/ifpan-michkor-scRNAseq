# ============================================
# 🔹 HEATMAPA Cramér’s V (signed)
#     z FDR < 0.10, winsoryzacją 5–95%,
#     i wartościami w komórkach
# ============================================

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)

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
# 3️⃣ Czyszczenie danych i filtrowanie genów po FDR < 0.10
# --------------------------------------------

heatmap_df <- heatmap_df %>%
  mutate(
    CramerV_signed = ifelse(is.na(CramerV_signed), 0, CramerV_signed),
    CramerV_signed = ifelse(FDR < 0.10, CramerV_signed, 0)
  )

genes_keep <- heatmap_df %>%
  filter(FDR < 0.10) %>%
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
# 5️⃣ Ustal kolejność kolumn
# --------------------------------------------

desired_col_order <- c(
  paste0("CORT_", cell_order),
  paste0("DEX_", cell_order)
)
desired_col_order <- desired_col_order[desired_col_order %in% colnames(heatmap_df_wide)]
heatmap_df_wide <- heatmap_df_wide[, desired_col_order, drop = FALSE]

# --------------------------------------------
# 6️⃣ Winsoryzacja 5–95%
# --------------------------------------------

vals <- as.numeric(heatmap_df_wide)
vals <- vals[!is.na(vals)]
q_low  <- quantile(vals, 0.05, na.rm = TRUE)
q_high <- quantile(vals, 0.95, na.rm = TRUE)

heatmap_df_wide_trimmed <- pmin(pmax(heatmap_df_wide, q_low), q_high)

# --------------------------------------------
# 7️⃣ Kolory i adnotacje
# --------------------------------------------

max_abs <- max(abs(c(q_low, q_high)))
col_fun <- colorRamp2(c(-max_abs, 0, max_abs),
                      c("navy", "white", "firebrick"))

treatment_labels <- ifelse(grepl("^CORT_", colnames(heatmap_df_wide_trimmed)), "CORT", "DEX")
treatment_colors <- c("CORT" = "#3E96FF", "DEX" = "#FF4E4E")

cluster_labels <- gsub("^(CORT_|DEX_)", "", colnames(heatmap_df_wide_trimmed))
cluster_colors <- c(
  "Astrocytes" = "#66c2a5",
  "Oligodendrocytes" = "#fc8d62",
  "Mature neurons"   = "#8da0cb"
)

top_annot <- HeatmapAnnotation(
  Treatment = treatment_labels,
  Cluster = cluster_labels,
  col = list(Treatment = treatment_colors, Cluster = cluster_colors),
  annotation_name_gp = gpar(fontsize = 9),
  simple_anno_size = unit(4, "mm")
)

# --------------------------------------------
# 8️⃣ Rysowanie heatmapy z wartościami
# --------------------------------------------

Heatmap(
  heatmap_df_wide_trimmed,
  name = "Cramér’s V (signed)",
  col = col_fun,
  top_annotation = top_annot,
  cluster_rows = TRUE,
  cluster_columns = FALSE,     # zachowaj kolejność CORT→DEX
  row_names_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 8),
  column_names_rot = 45,
  border = TRUE,
  
  # 🔸 rysowanie wartości w komórkach
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(
      sprintf("%.2f", heatmap_df_wide_trimmed[i, j]),
      x, y,
      gp = gpar(fontsize = 6.5, col = "black")
    )
  },
  
  heatmap_legend_param = list(
    title = "Cramér’s V",
    at = c(-max_abs, 0, max_abs),
    labels = c("↓", "0", "↑")
  )
)
