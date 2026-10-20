results_list_bulk <- bulk_summary_list_AllGenes %>% lapply(., function(x){
  analyze_gene_detection_signed(x, return_long = TRUE)})



bulk_summary_list_AllGenes$`all-well` %>% 
  filter(!(CTRL_percent < 5 & DEX_percent < 5)) %>% 
  analyze_gene_detection_signed(return_long = TRUE) %>% 
  filter(treatment == "DEX") %>% 
  mutate(fdr = p.adjust(p_chi2, method = "fdr")) %>% 
  filter(p_chi2 < 0.05)
  

results_list_bulk$`all-well` %>% 
  filter(treatment == "DEX") %>% 
  filter(p_chi2 < 0.05)



results_list %>% 
  lapply(., function(x){
    x %>%
      filter(treatment == "DEX") %>% 
      mutate(fdr = p.adjust(p_chi2, method = "fdr")) %>% 
      filter(fdr < 0.1) %>% 
      filter(direction == "up") %>% nrow
    
  }) %>% unlist %>% sort
