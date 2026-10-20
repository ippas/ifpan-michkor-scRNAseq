# ============================================
# 🔹 HEATMAPA — 3 typy komórek (poprawne FDR i jedna linia legend)
# ============================================

library(dplyr)
library(tidyr)
library(purrr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# --------------------------------------------
# 1️⃣ Definicje klastrów
# --------------------------------------------

# pełny zestaw klastrów do liczenia FDR
clusters_fdr_universe <- c(
  "Radial glial cells",
  "Oligodendrocyte precursor cells",
  "Immature neurons",
  "Astrocytes",
  "Mature neurons",
  "Oligodendrocytes"
)

# tylko 3 do wizualizacji
display_clusters <- c("Astrocytes", "Mature neurons", "Oligodendrocytes")

# --------------------------------------------
# 2️⃣ Łączenie danych i obliczanie FDR na całym zestawie
# --------------------------------------------

df_all <- results_list[clusters_fdr_universe] %>%
  purrr::imap_dfr(~ mutate(.x, cluster = .y)) %>%
  select(cluster, gene_symbol, treatment, CramerV_signed, p_chi2) %>%
  group_by(cluster, treatment) %>%
  mutate(FDR = p.adjust(p_chi2, method = "fdr")) %>%
  ungroup()

# --------------------------------------------
# 3️⃣ Wybór genów, które przeszły FDR < 0.2 w dowolnym klastrze
# --------------------------------------------

genes_keep <- df_all %>%
  filter(FDR < 0.20) %>%
  distinct(gene_symbol) %>%
  pull(gene_symbol)

# --------------------------------------------
# 4️⃣ Przygotowanie danych tylko dla 3 typów komórek (z zachowaniem spójnej listy genów)
# --------------------------------------------

df <- df_all %>%
  filter(cluster %in% display_clusters,
         gene_symbol %in% genes_keep) %>%
  mutate(
    CramerV_signed = ifelse(is.na(CramerV_signed), 0, CramerV_signed),
    CramerV_signed = ifelse(FDR < 0.20, CramerV_signed, 0)
  )

# --------------------------------------------
# 5️⃣ Przygotowanie macierzy do heatmapy
# --------------------------------------------

cell_order <- display_clusters

mat_text <- df %>%
  unite(col = "col_id", treatment, cluster, sep = "_") %>%
  select(gene_symbol, col_id, CramerV_signed) %>%
  pivot_wider(names_from = col_id, values_from = CramerV_signed) %>%
  tibble::column_to_rownames("gene_symbol") %>%
  as.matrix()

mat_fdr <- df %>%
  unite(col = "col_id", treatment, cluster, sep = "_") %>%
  select(gene_symbol, col_id, FDR) %>%
  pivot_wider(names_from = col_id, values_from = FDR) %>%
  tibble::column_to_rownames("gene_symbol") %>%
  as.matrix()

desired_col_order <- c(paste0("CORT_", cell_order), paste0("DEX_", cell_order))
desired_col_order <- desired_col_order[desired_col_order %in% colnames(mat_text)]

mat_text <- mat_text[, desired_col_order, drop = FALSE]
mat_fdr  <- mat_fdr[, desired_col_order, drop = FALSE]

# --------------------------------------------
# 6️⃣ Winsoryzacja kolorów (10–90%)
# --------------------------------------------

vals <- as.numeric(mat_text)
vals <- vals[!is.na(vals)]
q_low  <- quantile(vals, 0.10, na.rm = TRUE)
q_high <- quantile(vals, 0.90, na.rm = TRUE)

mat_color <- pmin(pmax(mat_text, q_low), q_high)

max_abs <- max(abs(c(q_low, q_high)))
col_fun <- colorRamp2(c(-max_abs, 0, max_abs),
                      c("navy", "white", "firebrick"))

# --------------------------------------------
# 7️⃣ Adnotacje i kolory
# --------------------------------------------

treatment_labels <- ifelse(grepl("^CORT_", colnames(mat_text)), "CORT", "DEX")
cluster_labels   <- gsub("^(CORT_|DEX_)", "", colnames(mat_text))

cluster_colors <- c(
  "Astrocytes"       = "#e5c494",
  "Mature neurons"   = "#b3b3b3",
  "Oligodendrocytes" = "#1f78b4"
)

top_annot <- HeatmapAnnotation(
  Treatment = treatment_labels,
  Cluster   = cluster_labels,
  col = list(
    Treatment = c("CORT" = "#3E96FF", "DEX" = "#FF4E4E"),
    Cluster   = cluster_colors
  ),
  annotation_name_gp = gpar(fontsize = 10),
  simple_anno_size = unit(5, "mm")
)

# --------------------------------------------
# 8️⃣ Funkcja gwiazdek dla FDR < 0.2 / 0.1 / 0.05
# --------------------------------------------

get_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.05) return("***")
  else if (p < 0.10) return("**")
  else if (p < 0.20) return("*")
  else return("")
}

# --------------------------------------------
# 9️⃣ Tworzenie heatmapy
# --------------------------------------------

ht <- Heatmap(
  mat_color,
  name = "Cramér’s V (signed)",
  col = col_fun,
  top_annotation = top_annot,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_gp = gpar(fontsize = 9),
  column_names_rot = 45,
  border = TRUE,
  
  cell_fun = function(j, i, x, y, w, h, fill) {
    val <- mat_text[i, j]
    fdr <- mat_fdr[i, j]
    stars <- get_stars(fdr)
    label <- sprintf("%.2f%s", val, stars)
    txt_col <- ifelse(abs(mat_color[i, j]) > (0.6 * max_abs), "white", "black")
    grid.text(label, x, y, gp = gpar(fontsize = 8.5, fontface = "bold", col = txt_col))
  }
)

# --------------------------------------------
# 🔟 Rysowanie legend w jednym wierszu (na górze)
# --------------------------------------------

# Główna legenda (skala)
lgd_heatmap <- Legend(
  title = "Cramér’s V",
  col_fun = col_fun,
  at = c(-max_abs, 0, max_abs),
  labels = c("↓", "0", "↑"),
  direction = "horizontal",
  title_position = "topcenter",
  title_gp = gpar(fontsize = 10, fontface = "bold"),
  labels_gp = gpar(fontsize = 9)
)

# Legenda dla treatmentów
lgd_treatment <- Legend(
  title = "Treatment",
  legend_gp = gpar(fill = c("CORT" = "#3E96FF", "DEX" = "#FF4E4E")),
  labels = c("CORT", "DEX"),
  nrow = 1
)

# Legenda dla klastrów
lgd_clusters <- Legend(
  title = "Cluster",
  legend_gp = gpar(fill = cluster_colors),
  labels = names(cluster_colors),
  nrow = 1
)

# Połączenie legend w jednym wierszu
combined_lgd <- packLegend(
  lgd_heatmap, lgd_treatment, lgd_clusters,
  direction = "horizontal"
)

# Rysowanie całości
draw(ht, heatmap_legend_list = combined_lgd, newpage = TRUE)
