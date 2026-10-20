rbind(summary_all_clusters_DEX, summary_all_clusters_CORT) %>% 
  filter(
    cluster %in% c("Astrocytes", "Radial glial cells", "Oligodendrocyte precursor cells",
                   "Oligodendrocytes", "Immature neurons", "Mature neurons")
  ) %>% 
  filter(logFC_threshold == 1) %>% 
  filter(class == "all") %>% 
  filter(cpm_threshold == 50) %>% 
  filter(p_value_threshold == 0.01) %>% 
  mutate(celltype = case_when(
    cluster == "Oligodendrocyte precursor cells" ~ "OPC",
    cluster == "Radial glial cells" ~ "RGS",
    TRUE ~ cluster
  )) %>% 
  mutate(label = paste0(celltype, "_", contrast, "_", direction,  "_edgerCPM", cpm_threshold, "P", p_value_threshold, "log2FC", logFC_threshold)) %>% 
  select(label, genes) %>% 
  split(.$label) %>% 
  lapply(., function(x){
    x %>% 
      .$genes %>% .[[1]]
  }) %>% 
  saveRDS("results/results_GRsignatures_30.11.2025/grSignatures_6cellTypes_edgerCPM50p0.01log2FC1.rds")


library(dplyr)
library(openxlsx)

library(dplyr)
library(openxlsx)

rbind(summary_all_clusters_DEX, summary_all_clusters_CORT) %>% 
  filter(
    cluster %in% c("Astrocytes", "Radial glial cells",
                   "Oligodendrocyte precursor cells", "Oligodendrocytes",
                   "Immature neurons", "Mature neurons"),
    logFC_threshold == 1,
    class == "all",
    cpm_threshold == 50,
    p_value_threshold == 0.01
  ) %>%
  mutate(celltype = case_when(
    cluster == "Oligodendrocyte precursor cells" ~ "OPC",
    cluster == "Radial glial cells" ~ "RGS",
    TRUE ~ cluster
  )) %>%
  mutate(label = paste0(
    celltype, "_", contrast, "_", direction,
    "_edgerCPM", cpm_threshold,
    "P", p_value_threshold,
    "log2FC", logFC_threshold
  )) %>%
  select(label, genes) %>%
  split(.$label) %>%
  lapply(function(x) x$genes[[1]]) %>%      # <── TU była przyczyna błędu
  {
    DEX_UP    <- unlist(. [grepl("DEX.*up",   names(.), ignore.case = TRUE)], use.names = FALSE)
    DEX_DOWN  <- unlist(. [grepl("DEX.*down", names(.), ignore.case = TRUE)], use.names = FALSE)
    CORT_UP   <- unlist(. [grepl("CORT.*up",  names(.), ignore.case = TRUE)], use.names = FALSE)
    CORT_DOWN <- unlist(. [grepl("CORT.*down",names(.), ignore.case = TRUE)], use.names = FALSE)
    
    # (opcjonalnie) unikalne i posortowane
    DEX_UP    <- sort(unique(DEX_UP))
    DEX_DOWN  <- sort(unique(DEX_DOWN))
    CORT_UP   <- sort(unique(CORT_UP))
    CORT_DOWN <- sort(unique(CORT_DOWN))
    
    wb <- createWorkbook()
    
    write_sheet <- function(sheet, vec) {
      addWorksheet(wb, sheet)
      writeData(wb, sheet, data.frame(GENES = vec))
      freezePane(wb, sheet, firstRow = TRUE)
      setColWidths(wb, sheet, cols = 1, widths = "auto")
    }
    
    write_sheet("DEX_UP",    DEX_UP)
    write_sheet("DEX_DOWN",  DEX_DOWN)
    write_sheet("CORT_UP",   CORT_UP)
    write_sheet("CORT_DOWN", CORT_DOWN)
    
    saveWorkbook(
      wb,
      "results/results_GRsignatures_30.11.2025/grSignatures_6cellTypes_edgerCPM50p0.01log2FC1.xlsx",
      overwrite = TRUE
    )
  }
