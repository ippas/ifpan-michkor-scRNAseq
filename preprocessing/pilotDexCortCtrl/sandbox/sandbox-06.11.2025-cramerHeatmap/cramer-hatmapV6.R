# ============================================
# 🔹 HEATMAPA Cramér’s V (signed)
#     FDR < 0.10, winsoryzacja KOLORU 5–95%
#     (etykiety pokazują ORYGINALNE wartości)
# ============================================

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)

# --------------------------------------------
# 1) Zbierz i przygotuj dane
# --------------------------------------------
cell_order <- c("Astrocytes", "Oligodendrocytes", "Mature neurons")

df <- results_list %>%
  purrr::imap_dfr(~ dplyr::mutate(.x, cluster = .y)) %>%
  select(cluster, gene_symbol, treatment, CramerV_signed, p_chi2) %>%
  filter(cluster %in% cell_order) %>%
  group_by(cluster, treatment) %>%
  mutate(FDR = p.adjust(p_chi2, method = "fdr")) %>%
  ungroup()

# wyzeruj nieistotne i zostaw tylko geny z ≥1 istotnym wynikiem (FDR < 0.10)
df <- df %>%
  mutate(CramerV_signed = ifelse(is.na(CramerV_signed), 0, CramerV_signed),
         CramerV_signed = ifelse(FDR < 0.10, CramerV_signed, 0))

genes_keep <- df %>% filter(FDR < 0.10) %>% distinct(gene_symbol) %>% pull(gene_symbol)
df <- df %>% filter(gene_symbol %in% genes_keep)

# --------------------------------------------
# 2) Szeroka macierz (oryginalne wartości do etykiet)
# --------------------------------------------
mat_text <- df %>%
  unite(col = "col_id", treatment, cluster, sep = "_") %>%
  select(gene_symbol, col_id, CramerV_signed) %>%
  pivot_wider(names_from = col_id, values_from = CramerV_signed) %>%
  tibble::column_to_rownames("gene_symbol") %>%
  as.matrix()

# wymuś kolejność kolumn: CORT → DEX oraz kolejność klastrów
desired_col_order <- c(paste0("CORT_", cell_order), paste0("DEX_", cell_order))
desired_col_order <- desired_col_order[desired_col_order %in% colnames(mat_text)]
mat_text <- mat_text[, desired_col_order, drop = FALSE]

# --------------------------------------------
# 3) Zrób macierz TYLKO do kolorowania (winsoryzacja)
# --------------------------------------------
vals <- as.numeric(mat_text); vals <- vals[!is.na(vals)]
q_low  <- quantile(vals, 0.05, na.rm = TRUE)
q_high <- quantile(vals, 0.95, na.rm = TRUE)

# przycięte wartości tylko do koloru
mat_color <- pmin(pmax(mat_text, q_low), q_high)

# skala symetryczna (tylko dla koloru/legendy)
max_abs <- max(abs(c(q_low, q_high)))
col_fun <- colorRamp2(c(-max_abs, 0, max_abs), c("navy", "white", "firebrick"))

# --------------------------------------------
# 4) Adnotacje kolumn
# --------------------------------------------
treatment_labels <- ifelse(grepl("^CORT_", colnames(mat_text)), "CORT", "DEX")
cluster_labels   <- gsub("^(CORT_|DEX_)", "", colnames(mat_text))

top_annot <- HeatmapAnnotation(
  Treatment = treatment_labels,
  Cluster   = cluster_labels,
  col = list(
    Treatment = c("CORT"="#3E96FF","DEX"="#FF4E4E"),
    Cluster   = c("Astrocytes"="#66c2a5","Oligodendrocytes"="#fc8d62","Mature neurons"="#8da0cb")
  ),
  annotation_name_gp = gpar(fontsize = 9),
  simple_anno_size = unit(4, "mm")
)

# --------------------------------------------
# 5) Heatmapa: kolor z mat_color, tekst z mat_text
# --------------------------------------------
Heatmap(
  mat_color,                              # ⬅️ kolory po winsoryzacji
  name = "Cramér’s V (signed)",
  col  = col_fun,
  top_annotation   = top_annot,
  cluster_rows     = TRUE,
  cluster_columns  = FALSE,               # zachowaj kolejność CORT → DEX
  row_names_gp     = gpar(fontsize = 7),
  column_names_gp  = gpar(fontsize = 8),
  column_names_rot = 45,
  border = TRUE,
  heatmap_legend_param = list(
    title  = "Cramér’s V",
    at     = c(-max_abs, 0, max_abs),
    labels = c("↓", "0", "↑")
  ),
  # 🔸 wartości w komórkach = ORYGINALNE (z mat_text)
  cell_fun = function(j, i, x, y, w, h, fill) {
    val <- mat_text[i, j]
    lab <- sprintf("%.2f", val)
    # opcjonalnie: lepszy kontrast napisu
    # jasny tekst na ciemnym tle i odwrotnie
    txt_col <- ifelse(abs(mat_color[i, j]) > (0.6*max_abs), "white", "black")
    grid.text(lab, x, y, gp = gpar(fontsize = 6.5, col = txt_col))
  }
)
