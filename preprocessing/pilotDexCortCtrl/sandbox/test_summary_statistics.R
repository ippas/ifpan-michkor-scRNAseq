res_global$combined %>% 
  filter(g_FDR < 0.05) %>% 
  filter(DEX_CTRL_logFC > 1 | DEX_CTRL_logFC < -1 | CORT_CTRL_logFC > 1 | CORT_CTRL_logFC < -1) %>% 
  head
  

# 1. Definiujemy kryteria istotności
sig_genes <- res_global$combined %>%
  filter(g_FDR < 0.05) %>% 
  filter(DEX_CTRL_logFC > 1 | DEX_CTRL_logFC < -1 | CORT_CTRL_logFC > 1 | CORT_CTRL_logFC < -1) %>% 
  mutate(
    sig_DEX_CTRL  = DEX_CTRL_FDR  < 0.05 & abs(DEX_CTRL_logFC)  > 1,
    sig_CORT_CTRL = CORT_CTRL_FDR < 0.05 & abs(CORT_CTRL_logFC) > 1,
    sig_CORT_DEX  = CORT_DEX_FDR  < 0.05 & abs(CORT_DEX_logFC)  > 1
  )



sig_genes %>%
  dplyr::count(sig_DEX_CTRL, sig_CORT_CTRL, sig_CORT_DEX) %>%
  arrange(desc(n))


# 2. Policz ile genów w każdej kombinacji
venn_table <- sig_genes %>%
  count(sig_DEX_CTRL, sig_CORT_CTRL, sig_CORT_DEX) %>%
  arrange(desc(n))

sig_genes %>%
  group_by(sig_DEX_CTRL, sig_CORT_CTRL, sig_CORT_DEX) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n))
