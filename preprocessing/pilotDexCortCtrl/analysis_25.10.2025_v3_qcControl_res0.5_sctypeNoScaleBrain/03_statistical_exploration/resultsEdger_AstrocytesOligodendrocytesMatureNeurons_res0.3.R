edgeR_anovaOneFactor_treatment_pseudobulk_res0.4_clusters16_celltypes11_results$`cluster_Mature neurons`$f_test_global %>% 
  filter(FDR < 0.01) %>% 
  filter(abs(logFC.treatmentDEX)  > 1 | abs(logFC.treatmentCORT) > 1) %>% 
  .$gene_symbol


edgeR_anovaOneFactor_treatment_pseudobulk_res0.4_clusters16_celltypes11_results$`cluster_Mature neurons`$f_test_global %>% 
  filter(FDR < 0.01) %>% 
  filter(abs(logFC.treatmentDEX)  > 1 | abs(logFC.treatmentCORT) > 1) %>% 
  .$gene_symbol



library(dplyr)
library(readr)
library(openxlsx)
library(purrr)

# 📁 Katalog i ścieżka do pliku wynikowego
out_dir <- "results/results_analysis_02.11.2025_v3_qcControl_res0.3_scTypeNoScaleBrain11"
out_xlsx <- file.path(out_dir, "05.11.2025_edgeRresults_clusters.xlsx")

# 🔹 Funkcja pomocnicza do pobierania i łączenia wyników
get_cluster_results <- function(cluster_name, results_list) {
  dex <- results_list[[cluster_name]]$DEX_vs_CTRL %>%
    filter(FDR < 0.01, abs(logFC) > 1) %>%
    mutate(regulationBy = "DEX")
  
  cort <- results_list[[cluster_name]]$CORT_vs_CTRL %>%
    filter(FDR < 0.01, abs(logFC) > 1) %>%
    mutate(regulationBy = "CORT")
  
  bind_rows(dex, cort)
}

# 🔹 Lista klastrów i ich przyjazne nazwy arkuszy
cluster_sheets <- list(
  "Mature neurons" = "cluster_Mature neurons",
  "Astrocytes" = "cluster_Astrocytes",
  "Oligodendrocytes" = "cluster_Oligodendrocytes"
)

# 🔹 Utwórz workbook
wb <- createWorkbook()

# 🔹 Iteracja po klastrach
for (sheet_name in names(cluster_sheets)) {
  cluster_key <- cluster_sheets[[sheet_name]]
  
  df <- get_cluster_results(cluster_key, edgeR_anovaOneFactor_treatment_pseudobulk_res0.4_clusters16_celltypes11_results)
  
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, df)
  freezePane(wb, sheet = sheet_name, firstRow = TRUE)
  setColWidths(wb, sheet_name, cols = 1:ncol(df), widths = "auto")
}

# 🔹 Zapisz plik
saveWorkbook(wb, out_xlsx, overwrite = TRUE)

cat("\n✅ Zapisano plik Excel z trzema zakładkami:\n", out_xlsx, "\n")
