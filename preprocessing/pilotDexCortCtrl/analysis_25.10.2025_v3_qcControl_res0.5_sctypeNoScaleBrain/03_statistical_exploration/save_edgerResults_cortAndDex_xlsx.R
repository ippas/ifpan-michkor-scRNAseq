# ##############################################################################
# ---- function ----
# ##############################################################################
library(dplyr)
library(purrr)
library(openxlsx)
library(stringr)

write_edger_results_to_xlsx <- function(
    res_list,
    contrast = "CORT_vs_CTRL",
    output_path,
    fdr_threshold = 0.05,
    logfc_threshold = 1,
    splitRegulation = FALSE,
    overwrite = TRUE
) {
  # Workbook i style
  wb <- createWorkbook()
  style_num3   <- createStyle(numFmt = "0.000")      # 3 miejsca po przecinku
  style_sci3   <- createStyle(numFmt = "0.000E+00")  # notacja naukowa
  style_header <- createStyle(textDecoration = "bold", halign = "center", valign = "center")
  
  # Czyszczenie nazw arkuszy + alias OPC
  clean_sheet_name <- function(name) {
    name <- str_remove(name, "^cluster_")
    name <- str_trim(name)
    name <- str_replace_all(name, " ", "_")
    if (identical(name, "Oligodendrocyte_precursor_cells")) name <- "OPC"
    substr(name, 1, 31)
  }
  
  # Iteracja po klastrach
  walk(names(res_list), function(cluster_name) {
    if (!contrast %in% names(res_list[[cluster_name]])) return(NULL)
    
    df <- res_list[[cluster_name]][[contrast]] %>%
      filter(FDR < fdr_threshold, abs(logFC) > logfc_threshold) %>%
      arrange(FDR, desc(abs(logFC)))
    
    if (nrow(df) == 0) return(NULL)
    
    base_name <- clean_sheet_name(cluster_name)
    
    # Helper do zapisu jednego arkusza
    write_one_sheet <- function(dat, sheet_nm) {
      addWorksheet(wb, sheet_nm)
      writeData(wb, sheet = sheet_nm, dat, headerStyle = style_header)
      freezePane(wb, sheet = sheet_nm, firstActiveRow = 2)
      
      cn <- colnames(dat)
      cols_num3 <- which(cn %in% c("logFC", "logCPM", "F"))
      cols_sci3 <- which(cn %in% c("PValue", "FDR"))
      
      if (length(cols_num3) > 0) {
        addStyle(wb, sheet = sheet_nm, style = style_num3,
                 cols = cols_num3, rows = 2:(nrow(dat) + 1), gridExpand = TRUE)
      }
      if (length(cols_sci3) > 0) {
        addStyle(wb, sheet = sheet_nm, style = style_sci3,
                 cols = cols_sci3, rows = 2:(nrow(dat) + 1), gridExpand = TRUE)
      }
    }
    
    if (!isTRUE(splitRegulation)) {
      write_one_sheet(df, base_name)
    } else {
      df_up   <- dplyr::filter(df, logFC > 0)
      df_down <- dplyr::filter(df, logFC < 0)
      if (nrow(df_up)   > 0) write_one_sheet(df_up,   substr(paste0(base_name, "_Up"),   1, 31))
      if (nrow(df_down) > 0) write_one_sheet(df_down, substr(paste0(base_name, "_Down"), 1, 31))
    }
  })
  
  # Zapis pliku
  saveWorkbook(wb, output_path, overwrite = overwrite)
  message("✅ Zapisano wyniki (", contrast, ") do: ", output_path)
}



write_edger_results_to_xlsx(
  res_list = edgeR_anovaOneFactor_treatment_pseudobulk_res0.5_clusters13_results,
  contrast = "CORT_vs_CTRL",
  output_path = "/home/mateusz/projects/ifpan-michkor-scRNAseq/results/results_analysis_25.10.2025_v3_qcControl_res0.5_sctypeNoScaleBrain/edgerResults/edgeR_cortVsCtrl_SumUpDown_res0.5_celltypes13_fdr0.05log2FC1.xlsx"
)


write_edger_results_to_xlsx(
  res_list = edgeR_anovaOneFactor_treatment_pseudobulk_res0.5_clusters13_results,
  contrast = "DEX_vs_CTRL",
  output_path = "/home/mateusz/projects/ifpan-michkor-scRNAseq/results/results_analysis_25.10.2025_v3_qcControl_res0.5_sctypeNoScaleBrain/edgerResults/edgeR_dexVsCtrl_SumUpDown_res0.5_celltypes13_fdr0.05log2FC1.xlsx"
)

write_edger_results_to_xlsx(
  res_list = edgeR_anovaOneFactor_treatment_pseudobulk_res0.5_clusters13_results,
  contrast = "CORT_vs_CTRL",
  splitRegulation = TRUE,
  output_path = "/home/mateusz/projects/ifpan-michkor-scRNAseq/results/results_analysis_25.10.2025_v3_qcControl_res0.5_sctypeNoScaleBrain/edgerResults/edgeR_cortVsCtrl_SeparateUpDown_res0.5_celltypes13_fdr0.05log2FC1.xlsx"
)

write_edger_results_to_xlsx(
  res_list = edgeR_anovaOneFactor_treatment_pseudobulk_res0.5_clusters13_results,
  contrast = "DEX_vs_CTRL",
  splitRegulation = TRUE,
  output_path = "/home/mateusz/projects/ifpan-michkor-scRNAseq/results/results_analysis_25.10.2025_v3_qcControl_res0.5_sctypeNoScaleBrain/edgerResults/edgeR_dexVsCtrl_SeparateUpDown_res0.5_celltypes13_fdr0.05log2FC1.xlsx"
)



