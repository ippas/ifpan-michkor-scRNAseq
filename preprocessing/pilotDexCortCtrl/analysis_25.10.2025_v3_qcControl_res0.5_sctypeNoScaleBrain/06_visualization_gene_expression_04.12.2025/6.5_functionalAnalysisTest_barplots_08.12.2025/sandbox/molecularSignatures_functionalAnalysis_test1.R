clusters_to_run <- paste0("cluster_", 1:5)

databases_to_run <- c(
  "GO_Biological_Process_2025",
  "Reactome_Pathways_2024",
  "MGI_Mammalian_Phenotype_Level_4_2024",
  "DisGeNET"
)

# Główna lista
enrichr_results <- list()

for (db in databases_to_run) {
  
  # podlista dla danej bazy
  enrichr_results[[db]] <- list()
  
  for (clust in clusters_to_run) {
    gene_list <- res_heatmap$genes_per_cluster[[clust]]
    
    df_res <- run_enrichr(
      gene_list = gene_list,
      database = db
    ) %>%
      dplyr::mutate(
        cluster = clust,
        database = db
      )
    
    # zapis do: enrichr_results[[baza]][[klaster]]
    enrichr_results[[db]][[clust]] <- df_res
  }
}




enrichr_results$Reactome_Pathways_2024 %>% 
  lapply(., function(x){
    x %>% filter(n_genes > 2,
                 pvalue < 1) %>% head(5)
  })


enrichr_results$MGI_Mammalian_Phenotype_Level_4_2024 %>% 
  lapply(., function(x){
    x %>% filter(n_genes > 2,
                 pvalue < 0.05) %>% head(10)
  })


enrichr_results$DisGeNET%>% 
  lapply(., function(x){
    x %>% filter(n_genes > 0,
                 pvalue < 1) %>% filter(term %in% c("Bipolar Disorder", "Unipolar Depression", "Schizophrenia"))
  }) %>% 
  lapply(., function(x){x$genes %>% strsplit(";")}) %>% 
  unname %>% unlist %>% unique() %>% length()

enrichr_results$GO_Biological_Process_2025 %>% 
  lapply(., function(x){
    x %>% filter(n_genes > 2,
                 pvalue < 1) %>% head(5)
  })



enrichr_results$DisGeNET%>% 
  lapply(., function(x){
    x %>% filter(n_genes > 0,
                 pvalue < 1) %>% filter(term %in% c("Bipolar Disorder", "Unipolar Depression", "Schizophrenia"))
  })
