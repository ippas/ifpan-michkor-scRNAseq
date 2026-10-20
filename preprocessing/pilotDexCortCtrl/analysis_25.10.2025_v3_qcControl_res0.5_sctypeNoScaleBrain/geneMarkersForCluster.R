Idents(integrated_pilotDexCortCtrl_filtered) <- "res0.3_RNA_noScaled_brain"


Idents(integrated_pilotDexCortCtrl_filtered) <- "res0.3_RNA_noScaled_brain"

markers_res03 <- FindAllMarkers(
  object = integrated_pilotDexCortCtrl_filtered,
  assay = "RNA",
  slot = "counts",       # <--- tu zmiana
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

markers_res03 %>% 
  mutate(cluster = as.character(cluster)) %>% .$cluster %>% unique

markers_res03 %>% 
  group_by(cluster) %>% 
  nest %>%
  mutate(n_genes)



split(markers_res03, markers_res03$cluster) %>% 
  lapply(., head, 100) %>% .$Astrocytes %>% 
  filter(grepl("AQP4", gene))


markers_res03 %>% 
  filter(cluster == "Astrocytes") %>% 
  # head(50) %>%
  # head(310) %>% 
  # tail(10)
  filter(grepl("EPHA6|SLC1A2|AQP4", gene))
# ##############################################################################
# ---- save to xlsx ----
# ##############################################################################
# ============================================================
# 📘 Eksport wyników FindAllMarkers do Excela
# - split na klastry
# - top 50 genów per cluster (posortowane po avg_log2FC malejąco)
# - kolumny: name_id, ensembl_id, gene_symbol
# - format p_val i p_val_adj (notacja naukowa, 3 cyfry)
# - zablokowany pierwszy wiersz w każdym arkuszu
# ============================================================

# ============================================================
# 📘 Eksport wyników FindAllMarkers do Excela (Top 50 per cluster)
# ============================================================

# 📦 Wymagane pakiety
library(dplyr)
library(openxlsx)

# 📁 Ścieżka wyjściowa
output_dir <- "/home/mateusz/projects/ifpan-michkor-scRNAseq/results/results_analysis_02.11.2025_v3_qcControl_res0.3_scTypeNoScaleBrain11/geneMarkers_FinMarkersSeurat/"

# 📄 Nazwa pliku wyjściowego
output_file <- file.path(
  output_dir,
  "seuratFindMarkersTop100_res0.3_clusters16_brainCell11_minPct0.25_log2FC0.25.xlsx"
)

# ============================================================
# 🔹 1. Podziel dane na klastry i wybierz top 50 genów
# ============================================================
cluster_list <- split(markers_res03, markers_res03$cluster) %>%
  lapply(function(df) {
    df %>%
      arrange(desc(avg_log2FC)) %>%
      head(50)
  })

# ============================================================
# 🔹 2. Utwórz workbook
# ============================================================
wb <- createWorkbook()

# ============================================================
# 🔹 3. Dla każdego klastra – przygotowanie tabeli i zapis
# ============================================================
for (cl_name in names(cluster_list)) {
  df <- cluster_list[[cl_name]]
  
  # --- name_id z kolumny gene
  df <- df %>%
    mutate(
      name_id = gene,
      ensembl_id = sub("-.*", "", gene),
      gene_symbol = sub(".*-", "", gene)
    )
  
  # --- usuń kolumnę gene
  df <- df %>% select(-any_of("gene"))
  
  # --- format p_val i p_val_adj
  df <- df %>%
    mutate(
      p_val = formatC(p_val, format = "e", digits = 3),
      p_val_adj = formatC(p_val_adj, format = "e", digits = 3),
      avg_log2FC = round(avg_log2FC, 3)   # ✅ trzy miejsca po przecinku
    )
  
  # --- kolejność kolumn
  df <- df %>%
    select(name_id, ensembl_id, gene_symbol, everything())
  
  # --- nazwa arkusza (bez znaków niedozwolonych)
  safe_name <- gsub("[^A-Za-z0-9_\\-]", "_", cl_name)
  
  # --- dodaj arkusz
  addWorksheet(wb, safe_name)
  
  # --- zapisz dane + autofiltr
  writeData(wb, sheet = safe_name, x = df, withFilter = TRUE)
  
  # --- zamroź pierwszy wiersz
  freezePane(wb, sheet = safe_name, firstActiveRow = 2)
  
  # --- dostosuj szerokość kolumn ensembl_id i gene_symbol
  setColWidths(wb, sheet = safe_name, cols = which(colnames(df) %in% c("ensembl_id", "gene_symbol")), widths = "auto")
}

# ============================================================
# 💾 4. Zapisz workbook
# ============================================================
saveWorkbook(wb, output_file, overwrite = TRUE)

cat("✅ Zapisano plik Excel:\n", output_file, "\n")
