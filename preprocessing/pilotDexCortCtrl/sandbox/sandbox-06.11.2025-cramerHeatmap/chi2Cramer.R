# ##############################################################################
# ---- ususes data ----
# ##############################################################################

cluster_summary_list_GlobalNeural$Astrocytes


# ============================================
# 🔹 Funkcja: testy chi2, OR i Cramér’s V
# ============================================



# ============================================
# 🔹 Analiza chi2, OR, Cramér’s V i kierunek efektu
# ============================================

# ============================================
# 🔹 Przykład użycia:
# ============================================

# Załóżmy że masz ramkę danych jak na zrzucie ekranu
# astro_df <- read.csv("twoj_plik.csv", sep="\t")

results_df <- analyze_gene_detection_signed(cluster_summary_list_GlobalNeural$`Mature neurons`, return_long = TRUE)

results_list <- cluster_summary_list_GlobalNeural %>% lapply(., function(x){
  analyze_gene_detection_signed(x, return_long = TRUE)})
# Wyniki


results_df %>% 
  filter(p_chi2 < 0.05) %>% 
  filter(treatment == "DEX") %>% 
  filter(CramerV > 0.05) %>% as.data.frame()



library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(purrr)
library(tibble)

# ============================================
# 🔹 1. Połączenie wszystkich wyników w jedną ramkę danych
# ============================================

cell_order <- c("Astrocytes", "Oligodendrocytes", "Mature neurons")

heatmap_df <- results_list %>%
  purrr::imap_dfr(~ mutate(.x, cluster = .y)) %>%
  select(cluster, gene_symbol, treatment, CramerV_signed, p_chi2) %>%
  filter(cluster %in% cell_order)

# ============================================
# 🔹 2. Wyczyść i przekształć dane
# ============================================

heatmap_df <- heatmap_df %>%
  mutate(
    CramerV_signed = ifelse(is.na(CramerV_signed), 0, CramerV_signed),
    CramerV_signed = ifelse(p_chi2 < 0.05, CramerV_signed, 0)
  )

# scal treatment i cluster w jedną nazwę kolumny
heatmap_df_wide <- heatmap_df %>%
  unite(col = "col_id", treatment, cluster, sep = "_") %>%
  select(gene_symbol, col_id, CramerV_signed) %>%
  pivot_wider(names_from = col_id, values_from = CramerV_signed) %>%
  tibble::column_to_rownames("gene_symbol") %>%
  as.matrix()

# ============================================
# 🔹 3. Wymuszenie kolejności kolumn
# ============================================

# kolejność: najpierw Astrocytes, potem Oligo, potem Neurons
desired_col_order <- c(
  paste0("CORT_", cell_order),
  paste0("DEX_", cell_order)
)

# zachowaj tylko kolumny, które faktycznie istnieją
desired_col_order <- desired_col_order[desired_col_order %in% colnames(heatmap_df_wide)]
heatmap_df_wide <- heatmap_df_wide[, desired_col_order, drop = FALSE]

# ============================================
# 🔹 4. Przygotowanie adnotacji kolumn
# ============================================

# kolory treatmentów (góra)
treatment_labels <- ifelse(grepl("^CORT_", colnames(heatmap_df_wide)), "CORT", "DEX")
treatment_colors <- c("CORT" = "#3E96FF", "DEX" = "#FF4E4E")

# kolory typów komórek (drugi poziom adnotacji)
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

# ============================================
# 🔹 5. Definicja kolorów dla wartości
# ============================================

col_fun <- colorRamp2(c(-0.5, 0, 0.5),
                      c("dodgerblue3", "white", "firebrick3"))

# ============================================
# 🔹 6. Rysowanie heatmapy
# ============================================

Heatmap(
  heatmap_df_wide,
  name = "Cramér’s V (signed)",
  col = col_fun,
  top_annotation = top_annot,
  cluster_rows = TRUE,
  cluster_columns = FALSE,  # nie zmieniać kolejności kolumn
  row_names_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 8),
  column_names_rot = 45,
  heatmap_legend_param = list(
    title = "Cramér’s V",
    at = c(-0.5, 0, 0.5),
    labels = c("↓", "0", "↑")
  )
)

