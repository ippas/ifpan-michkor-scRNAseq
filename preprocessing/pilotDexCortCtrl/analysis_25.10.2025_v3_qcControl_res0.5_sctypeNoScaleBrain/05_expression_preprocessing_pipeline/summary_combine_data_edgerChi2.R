combined_summaryEdgerChi2_res0.3RNAnoScaledBrain$`Radial glial cells` %>% 
  filter(max_percent > 5) %>% 
  filter(log2FC_CORT_CPM > 1)


combined_summaryEdgerChi2_res0.3RNAnoScaledBrain$`Radial glial cells` %>% 
  mutate(
    category_percent = case_when(
      max_percent < 5 ~ "percent_<5",
      max_percent >= 5 & max_percent < 10  ~ "percent_5-10",
      max_percent >= 10 & max_percent < 20  ~ "percent_10-20",
      max_percent >= 20 & max_percent < 30  ~ "percent_20-30",
      max_percent >= 30 & max_percent < 50  ~ "percent_30-50",
      max_percent >= 50 ~ "percent_>50",
      TRUE         ~ "other"    # domyślna wartość podobna do else
    )
  ) %>% 
  mutate(
    category_percent = factor(
      category_percent,
      levels = c(
        "percent_<5",
        "percent_5-10",
        "percent_10-20",
        "percent_20-30",
        "percent_30-50",
        "percent_>50"
      )
    )
  ) %>% 
  mutate(gr_signature = ifelse(gene_symbol %in% c(genes_global, genes_neural), TRUE, FALSE)) %>% 
  filter(gr_signature) %>% 
  group_by(category_percent) %>% 
  nest() %>% 
  mutate(n_genes = map(data, ~ .x %>% .$gene_symbol %>% unique %>% length)) %>% 
  mutate(log2FC_CPM_CORTgreater2 = map(data, ~ .x %>% filter(log2FC_CORT_CPM > 2))) %>% 
  mutate(log2FC_CPM_DEXgreater2 = map(data, ~ .x %>% filter(log2FC_DEX_CPM > 2))) %>% 
  mutate(log2FC_CPM_CORTgreater1 = map(data, ~ .x %>% filter(log2FC_CORT_CPM > 1))) %>% 
  mutate(log2FC_CPM_DEXgreater1 = map(data, ~ .x %>% filter(log2FC_DEX_CPM > 1))) %>% 
  mutate(log2FC_CPM_CORTgreater0.5 = map(data, ~ .x %>% filter(log2FC_CORT_CPM > 0.5))) %>% 
  mutate(log2FC_CPM_DEXgreater0.5 = map(data, ~ .x %>% filter(log2FC_DEX_CPM > 0.5))) %>% 
  # edgeR log2FC
  mutate(edger_CORT_log2FCgreater0.5 = map(data, ~ .x %>% filter(log2FC_edgerCORT > 0.5))) %>% 
  mutate(edger_CORT_log2FCgreater1 = map(data, ~ .x %>% filter(log2FC_edgerCORT > 1))) %>% 
  mutate(edger_CORT_log2FCgreater2 = map(data, ~ .x %>% filter(log2FC_edgerCORT > 2))) %>% 
  mutate(edger_DEX_log2FCgreater0.5 = map(data, ~ .x %>% filter(log2FC_edgerDEX > 0.5))) %>% 
  mutate(edger_DEX_log2FCgreater1 = map(data, ~ .x %>% filter(log2FC_edgerDEX > 1))) %>% 
  mutate(edger_DEX_log2FCgreater2 = map(data, ~ .x %>% filter(log2FC_edgerDEX > 2))) %>% 
  # edgeR; lof2FC, FDR
  mutate(edger_CORT_log2FCgreater0.5FDR0.05 = map(data, ~ .x %>% filter(log2FC_edgerCORT > 0.5 & FDR_edgerCORT < 0.05))) %>% 
  mutate(edger_CORT_log2FCgreater0.5FDR0.01 = map(data, ~ .x %>% filter(log2FC_edgerCORT > 0.5 & FDR_edgerCORT < 0.01))) %>% 
  mutate(edger_CORT_log2FCgreater0.5FDR0.001 = map(data, ~ .x %>% filter(log2FC_edgerCORT > 0.5 & FDR_edgerCORT < 0.001))) %>% 
  mutate(edger_DEX_log2FCgreater0.5FDR0.05 = map(data, ~ .x %>% filter(log2FC_edgerDEX > 0.5 & FDR_edgerCORT < 0.05))) %>% 
  mutate(edger_DEX_log2FCgreater0.5FDR0.01 = map(data, ~ .x %>% filter(log2FC_edgerDEX > 0.5 & FDR_edgerCORT < 0.01))) %>% 
  mutate(edger_DEX_log2FCgreater0.5FDR0.001 = map(data, ~ .x %>% filter(log2FC_edgerDEX > 0.5 & FDR_edgerCORT < 0.001))) %>% 
  # chi2
  mutate(chi2_CORT_p0.05 = map(data, ~ .x %>% filter(p_chi2_CTRL_vs_CORT < 0.05, CramerV_signed_CTRL_vs_CORT > 0))) %>% 
  mutate(chi2_CORT_p0.01 = map(data, ~ .x %>% filter(p_chi2_CTRL_vs_CORT < 0.01, CramerV_signed_CTRL_vs_CORT > 0))) %>% 
  mutate(chi2_CORT_p0.001 = map(data, ~ .x %>% filter(p_chi2_CTRL_vs_CORT < 0.001, CramerV_signed_CTRL_vs_CORT > 0))) %>% 
  mutate(chi2_DEX_p0.05 = map(data, ~ .x %>% filter(p_chi2_CTRL_vs_DEX < 0.05, CramerV_signed_CTRL_vs_DEX > 0))) %>% 
  mutate(chi2_DEX_p0.01 = map(data, ~ .x %>% filter(p_chi2_CTRL_vs_DEX < 0.01, CramerV_signed_CTRL_vs_DEX > 0))) %>% 
  mutate(chi2_DEX_p0.001 = map(data, ~ .x %>% filter(p_chi2_CTRL_vs_DEX < 0.001, CramerV_signed_CTRL_vs_DEX > 0))) %>% 
  select(c(chi2_DEX_p0.05, chi2_DEX_p0.01, chi2_DEX_p0.001,
           edger_DEX_log2FCgreater0.5FDR0.05, edger_DEX_log2FCgreater0.5FDR0.01, edger_DEX_log2FCgreater0.5FDR0.001,
           log2FC_CPM_CORTgreater0.5, edger_DEX_log2FCgreater0.5)) %>% 
  arrange(category_percent) %>% 
  select(category_percent, chi2_DEX_p0.01) %>% 
  unnest %>% 
  as.data.frame() %>% 
  filter(log2FC_DEX_CPM > 0.5)
  
  





  
