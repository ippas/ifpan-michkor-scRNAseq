run_enrichr <- function(gene_list, database) {
  results <- enrichR::enrichr(gene_list, database)[[1]]
  
  results %>% 
    dplyr::select(!c("Old.P.value", "Old.Adjusted.P.value")) %>% 
    mutate(n_genes = as.integer(str_split_fixed(Overlap, "/", 2)[, 1]),
           set_size = as.integer(str_split_fixed(Overlap, "/", 2)[, 2])) %>%
    set_colnames(c("term", "overlap", "pvalue", "fdr", "odds_ratio", "combined_score", "genes", "n_genes", "set_size")) %>% 
    dplyr::select(c(term, overlap, n_genes, set_size, everything())) -> results
  
  return(results)
}
