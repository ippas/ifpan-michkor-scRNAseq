# ##############################################################################
# ---- prepare data CPM, log2ratio ----
# ##############################################################################
log2FC_normCPM_res0.3RNAnoScaledBrain$Astrocytes %>% 
  dplyr::rename(meanCPM_CTRL = "mean_CTRL") %>% 
  dplyr::rename(meanCPM_CORT = "mean_CORT") %>% 
  dplyr::rename(meanCPM_DEX = "mean_DEX") %>%
  dplyr::rename(log2FC_CORT_CPM = "log2FC_CORT") %>% 
  dplyr::rename(log2FC_DEX_CPM = "log2FC_DEX") %>% 
  dplyr::mutate(max_CPM = pmax(meanCPM_CTRL, meanCPM_CORT, meanCPM_DEX, na.rm = TRUE)) 


# ##############################################################################
# ---- prepare edger data ----
# ##############################################################################
edgerResults_res0.3RNAnoScaledBrain$Astrocytes$CORT_vs_CTRL %>% 
  dplyr::select(-c(ensembl_id, gene_symbol)) %>% 
  dplyr::rename(gene_symbol = "name_id") %>% 
  dplyr::rename(log2FC_edgerCORT = "logFC") %>% 
  dplyr::rename(log2CPM_edgerCORT = "logCPM") %>% 
  dplyr::rename(F_edgerCORT = "F") %>% 
  dplyr::rename(PValue_edgerCORT = "PValue") %>% 
  dplyr::rename(FDR_edgerCORT = "FDR")
  
edgerResults_res0.3RNAnoScaledBrain$Astrocytes$DEX_vs_CTRL %>%
  dplyr::select(-c(ensembl_id, gene_symbol)) %>% 
  dplyr::rename(gene_symbol = "name_id") %>% 
  dplyr::rename(log2FC_edgerDEX = "logFC") %>% 
  dplyr::rename(log2CPM_edgerDEX = "logCPM") %>% 
  dplyr::rename(F_edgerDEX = "F") %>% 
  dplyr::rename(PValue_edgerDEX = "PValue") %>% 
  dplyr::rename(FDR_edgerDEX = "FDR") %>% 
  filter(gene_symbol == "STON1-GTF2A1L")

# ##############################################################################
# ---- prepare binary data----
# ##############################################################################

chi2Results_res0.3RNAnoScaledBrain$Astrocytes %>% 
  select(-cluster) %>% 
  dplyr::mutate(max_percent = pmax(CTRL_percent, CORT_percent, DEX_percent, na.rm = TRUE)) 



# ##############################################################################
# ---- combine data ----
# ##############################################################################
log2FC_normCPM_res0.3RNAnoScaledBrain$Astrocytes %>% 
  dplyr::rename(meanCPM_CTRL = "mean_CTRL") %>% 
  dplyr::rename(meanCPM_CORT = "mean_CORT") %>% 
  dplyr::rename(meanCPM_DEX = "mean_DEX") %>%
  dplyr::rename(log2FC_CORT_CPM = "log2FC_CORT") %>% 
  dplyr::rename(log2FC_DEX_CPM = "log2FC_DEX") %>% 
  dplyr::mutate(max_CPM = pmax(meanCPM_CTRL, meanCPM_CORT, meanCPM_DEX, na.rm = TRUE)) %>% 
  left_join(.,
            edgerResults_res0.3RNAnoScaledBrain$Astrocytes$CORT_vs_CTRL %>% 
              dplyr::select(-c(ensembl_id, gene_symbol)) %>% 
              dplyr::rename(gene_symbol = "name_id") %>% 
              dplyr::rename(log2FC_edgerCORT = "logFC") %>% 
              dplyr::rename(log2CPM_edgerCORT = "logCPM") %>% 
              dplyr::rename(F_edgerCORT = "F") %>% 
              dplyr::rename(PValue_edgerCORT = "PValue") %>% 
              dplyr::rename(FDR_edgerCORT = "FDR"),
            by = "gene_symbol") %>% 
  left_join(.,
            edgerResults_res0.3RNAnoScaledBrain$Astrocytes$DEX_vs_CTRL %>%
              dplyr::select(-c(ensembl_id, gene_symbol)) %>% 
              dplyr::rename(gene_symbol = "name_id") %>% 
              dplyr::rename(log2FC_edgerDEX = "logFC") %>% 
              dplyr::rename(log2CPM_edgerDEX = "logCPM") %>% 
              dplyr::rename(F_edgerDEX = "F") %>% 
              dplyr::rename(PValue_edgerDEX = "PValue") %>% 
              dplyr::rename(FDR_edgerDEX = "FDR"),
            by = "gene_symbol") %>% 
  left_join(.,
            chi2Results_res0.3RNAnoScaledBrain$Astrocytes %>% 
              select(-cluster) %>% 
              dplyr::mutate(max_percent = pmax(CTRL_percent, CORT_percent, DEX_percent, na.rm = TRUE)),
            by = "gene_symbol")
