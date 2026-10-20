# ============================================
# 🔹 HEATMAPA Cramér’s V (signed)
#     FDR < 0.10 z gwiazdkami istotności
#     (winsoryzacja koloru 10–90%)
#     + legenda pod spodem i większy tekst
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

df <- results_list %>%
  purrr::imap_dfr(~ mutate(.x, cluster = .y)) %>%
  select(cluster, gene_symbol, treatment, CramerV_signed, p_chi2) %>%
  filter(cluster %in% cell_order)

# --------------------------------------------
# 2️⃣ Obliczanie FDR w obrębie każdej kombinacji
# --------------------------------------------

df <- df %>%
  group_by(cluster, treatment) %>%
  mutate(FDR = p.adjust(p_chi2, method = "fdr")) %>%
  ungroup()

# --------------------------------------------
# 3️⃣ Czyszczenie danych i filtracja po FDR < 0.10
# --------------------------------------------

df <- df %>%
  mutate(
    CramerV_signed = ifelse(is.na(CramerV_signed), 0, CramerV_signed),
    CramerV_signed = ifelse(FDR < 0.10, CramerV_signed, 0)
  )

genes_keep <- df %>% filter(FDR < 0.10) %>% distinct(gene_symbol) %>% pull(gene_symbol)
df <- df %>% filter(gene_symbol %in% genes_keep)

# --------------------------------------------
# 4️⃣ Przygotowanie macierzy
# --------------------------------------------

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
# 5️⃣ Winsoryzacja kolorów (10–90%)
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
# 6️⃣ Adnotacje
# --------------------------------------------

treatment_labels <- ifelse(grepl("^CORT_", colnames(mat_text)), "CORT", "DEX")
cluster_labels   <- gsub("^(CORT_|DEX_)", "", colnames(mat_text))

top_annot <- HeatmapAnnotation(
  Treatment = treatment_labels,
  Cluster   = cluster_labels,
  col = list(
    Treatment = c("CORT" = "#3E96FF", "DEX" = "#FF4E4E"),
    Cluster   = c("Astrocytes" = "#66c2a5",
                  "Oligodendrocytes" = "#fc8d62",
                  "Mature neurons" = "#8da0cb")
  ),
  annotation_name_gp = gpar(fontsize = 10),
  simple_anno_size = unit(5, "mm")
)

# --------------------------------------------
# 7️⃣ Funkcja gwiazdek
# --------------------------------------------

get_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  else if (p < 0.05) return("**")
  else if (p < 0.10) return("*")
  else return("")
}

# --------------------------------------------
# 8️⃣ Heatmapa
# --------------------------------------------

ht <- Heatmap(
  mat_color,
  name = "Cramér’s V (signed)",
  col = col_fun,
  top_annotation = top_annot,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),  # 🔹 większy font nazw genów
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 45,
  border = TRUE,
  
  heatmap_legend_param = list(
    title = "Cramér’s V",
    at = c(-max_abs, 0, max_abs),
    labels = c("↓", "0", "↑"),
    legend_direction = "horizontal",  # 🔹 pozioma legenda
    title_position = "topcenter",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9)
  ),
  
  # 🔸 Teksty i gwiazdki
  cell_fun = function(j, i, x, y, w, h, fill) {
    val <- mat_text[i, j]
    fdr <- mat_fdr[i, j]
    stars <- get_stars(fdr)
    label <- sprintf("%.2f%s", val, stars)
    txt_col <- ifelse(abs(mat_color[i, j]) > (0.6 * max_abs), "white", "black")
    grid.text(label, x, y, gp = gpar(fontsize = 8.5, fontface = "bold", col = txt_col)) # 🔹 większy font w kwadracie
  }
)

# --------------------------------------------
# 9️⃣ Rysowanie z legendą pod spodem
# --------------------------------------------

draw(ht, heatmap_legend_side = "bottom")
