################################################################################
# 1. Przygotowanie obiektu listowego dla wyników połączonych
################################################################################

cluster_names <- names(log2FC_normCPM_res0.3RNAnoScaledBrain)
combined_summaryEdgerChi2_res0.3RNAnoScaledBrain <- vector("list", length(cluster_names))
names(combined_summaryEdgerChi2_res0.3RNAnoScaledBrain) <- cluster_names


################################################################################
# 2. Pętla łącząca dane dla każdego klastra
################################################################################

for (i in seq_along(cluster_names)) {
  
  cl <- cluster_names[i]
  cat("🔗 Łączenie danych dla klastra:", cl, "\n")
  
  # =====================================================================
  # A) Pobranie surowych CPM z normCPM_res0.3_RNA_noScaledBrain
  #    Zamiana macierzy dgeMatrix -> dataframe, rownames -> gene_symbol,
  #    kolumny A1-A12 na swoim miejscu (bez pivotowania / sumowania).
  # =====================================================================
  
  df_Acols <- normCPM_res0.3_RNA_noScaledBrain[[cl]] %>%
    as.matrix() %>%                      # konwersja DGE -> zwykła macierz
    as.data.frame() %>%                  # zamiana na dataframe
    tibble::rownames_to_column("gene_symbol")  # rownames -> gene_symbol
  
  
  # =====================================================================
  # B) CPM + log2ratio (średnie CPM per grupa + log2FC na podstawie CPM)
  # =====================================================================
  
  df_cpm <- log2FC_normCPM_res0.3RNAnoScaledBrain[[cl]] %>%
    dplyr::rename(
      meanCPM_CTRL     = "mean_CTRL",
      meanCPM_CORT     = "mean_CORT",
      meanCPM_DEX      = "mean_DEX",
      log2FC_CORT_CPM  = "log2FC_CORT",
      log2FC_DEX_CPM   = "log2FC_DEX"
    ) %>%
    dplyr::mutate(max_CPM = pmax(meanCPM_CTRL, meanCPM_CORT, meanCPM_DEX, na.rm = TRUE))
  
  
  # =====================================================================
  # C) Wyniki edgeR CORT_vs_CTRL
  # =====================================================================
  
  df_edgerCORT <- edgerResults_res0.3RNAnoScaledBrain[[cl]]$CORT_vs_CTRL %>%
    select(-c(gene_symbol, ensembl_id)) %>% 
    dplyr::rename(
      gene_symbol       = "name_id",
      log2FC_edgerCORT  = "logFC",
      log2CPM_edgerCORT = "logCPM",
      F_edgerCORT       = "F",
      PValue_edgerCORT  = "PValue",
      FDR_edgerCORT     = "FDR"
    )
  
  
  # =====================================================================
  # D) Wyniki edgeR DEX_vs_CTRL
  # =====================================================================
  
  df_edgerDEX <- edgerResults_res0.3RNAnoScaledBrain[[cl]]$DEX_vs_CTRL %>%
    select(-c(gene_symbol, ensembl_id)) %>% 
    dplyr::rename(
      gene_symbol       = "name_id",
      log2FC_edgerDEX   = "logFC",
      log2CPM_edgerDEX  = "logCPM",
      F_edgerDEX        = "F",
      PValue_edgerDEX   = "PValue",
      FDR_edgerDEX      = "FDR"
    )
  
  
  # =====================================================================
  # E) Wyniki chi-kwadrat detekcji (z wyliczeniem max_percent)
  # =====================================================================
  
  df_chi2 <- chi2Results_res0.3RNAnoScaledBrain[[cl]] %>%
    dplyr::select(-cluster) %>%
    dplyr::mutate(max_percent = pmax(CTRL_percent, CORT_percent, DEX_percent, na.rm = TRUE))
  
  
  # =====================================================================
  # F) Finalne łączenie (kolejność: A1–A12 → CPM+logratio → edgeR → chi2)
  # =====================================================================
  
  combined_summaryEdgerChi2_res0.3RNAnoScaledBrain[[cl]] <- df_Acols %>%
    dplyr::left_join(df_cpm,       by = "gene_symbol") %>%
    dplyr::left_join(df_edgerCORT, by = "gene_symbol") %>%
    dplyr::left_join(df_edgerDEX,  by = "gene_symbol") %>%
    dplyr::left_join(df_chi2,      by = "gene_symbol")
}

cat("✅ Wszystkie klastry zostały połączone.\n")


################################################################################
# 3. Zapis RDS na dysku
################################################################################

saveRDS(
  combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
  file = "data/pilotDexCortCtrl/r_outputs/statistics_results/combine_summary_edger_chi2/combined_summaryData_edgeRstatistcs_chi2statistics_res0.3RNAnoScaledBrain_25.11.2025.RDS"
)

cat("💾 Plik został zapisany pomyślnie.\n")
