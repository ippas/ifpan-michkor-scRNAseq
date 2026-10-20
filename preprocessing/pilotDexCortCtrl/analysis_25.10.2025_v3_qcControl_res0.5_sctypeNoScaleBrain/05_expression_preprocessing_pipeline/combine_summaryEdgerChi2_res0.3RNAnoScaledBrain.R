cluster_names <- names(log2FC_normCPM_res0.3RNAnoScaledBrain)
combined_summaryEdgerChi2_res0.3RNAnoScaledBrain <- vector("list", length(cluster_names))
names(combined_summaryEdgerChi2_res0.3RNAnoScaledBrain) <- cluster_names

for (i in seq_along(cluster_names)) {
  cl <- cluster_names[i]
  
  cat("🔗 Łączenie danych dla klastra:", cl, "\n")
  
  # ---- CPM + log2ratio ----
  df_cpm <- log2FC_normCPM_res0.3RNAnoScaledBrain[[cl]] %>%
    dplyr::rename(
      meanCPM_CTRL = "mean_CTRL",
      meanCPM_CORT = "mean_CORT",
      meanCPM_DEX  = "mean_DEX",
      log2FC_CORT_CPM = "log2FC_CORT",
      log2FC_DEX_CPM  = "log2FC_DEX"
    ) %>%
    dplyr::mutate(max_CPM = pmax(meanCPM_CTRL, meanCPM_CORT, meanCPM_DEX, na.rm = TRUE))
  
  # ---- edgeR CORT ----
  df_edgerCORT <- edgerResults_res0.3RNAnoScaledBrain[[cl]]$CORT_vs_CTRL %>%
    dplyr::select(-c(ensembl_id, gene_symbol)) %>%
    dplyr::rename(
      gene_symbol       = "name_id",
      log2FC_edgerCORT  = "logFC",
      log2CPM_edgerCORT = "logCPM",
      F_edgerCORT       = "F",
      PValue_edgerCORT  = "PValue",
      FDR_edgerCORT     = "FDR"
    )
  
  # ---- edgeR DEX ----
  df_edgerDEX <- edgerResults_res0.3RNAnoScaledBrain[[cl]]$DEX_vs_CTRL %>%
    dplyr::select(-c(ensembl_id, gene_symbol)) %>%
    dplyr::rename(
      gene_symbol       = "name_id",
      log2FC_edgerDEX   = "logFC",
      log2CPM_edgerDEX  = "logCPM",
      F_edgerDEX        = "F",
      PValue_edgerDEX   = "PValue",
      FDR_edgerDEX      = "FDR"
    )
  
  # ---- chi2 binary ----
  df_chi2 <- chi2Results_res0.3RNAnoScaledBrain[[cl]] %>%
    select(-cluster) %>%
    dplyr::mutate(max_percent = pmax(CTRL_percent, CORT_percent, DEX_percent, na.rm = TRUE))
  
  # ---- JOIN ----
  combined_summaryEdgerChi2_res0.3RNAnoScaledBrain [[cl]] <- df_cpm %>%
    left_join(df_edgerCORT, by = "gene_symbol") %>%
    left_join(df_edgerDEX , by = "gene_symbol") %>%
    left_join(df_chi2,      by = "gene_symbol")
}

cat("✅ Wszystkie klastry zostały połączone.\n")



saveRDS(combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
        file = "data/pilotDexCortCtrl/r_outputs/statistics_results/combine_summary_edger_chi2/combined_summaryData_edgeRstatistcs_chi2statistics_res0.3RNAnoScaledBrain_24.11.2025.RDS")
