BiocManager::install("fgsea")

library(fgsea)


# ##############################################################################
# ---- read data ----
# ##############################################################################
load("data/pilotDexCortCtrl/r_outputs/cluster_markers/clusterMarkersAndPlotsTopAndEnrichr_seurat_clusters_res0.13_fdr0.01logFC1minpct0.5.RData")

cluster_markers_enrichr_results$markers_up_list$cluster_0 %>% head


# ##############################################################################
# ---- prepare data ----
# ##############################################################################

cluster_markers_enrichr_results$markers_up_list %>%
  lapply(function(df) {
    df %>%
      select(name_id, gene_symbol) %>%
      mutate(rank = rev(seq_len(nrow(df)))) %>% 
      select(gene_symbol, rank) %>% 
      with(setNames(rank, gene_symbol))
  })-> cluster_markers_test

cluster_markers_test$cluster_0

panglaoDB_all_markers %>% 
  select(c(official.gene.symbol, cell.type)) %>% 
  mutate(cell.type = str_replace_all(cell.type, " ", "_")) %>% 
  group_by(cell.type) %>% 
  nest() %>% 
  mutate(n_genes = map(data, ~ .x %>% nrow)) %>% 
  unnest(n_genes) %>% 
  filter(n_genes >= 10) %>% 
  select(-n_genes) %>% 
  unnest(data) %>% 
  split(.$cell.type) %>% 
  lapply(function(x) unique(x$official.gene.symbol)) -> panglaoDB_list_AllMarkers


panglaoDB_all_markers %>% 
  filter(organ == "Brain") %>% 
  select(c(official.gene.symbol, cell.type)) %>% 
  mutate(cell.type = str_replace_all(cell.type, " ", "_")) %>% 
  group_by(cell.type) %>% 
  nest() %>% 
  mutate(n_genes = map(data, ~ .x %>% nrow)) %>% 
  unnest(n_genes) %>% 
  filter(n_genes >= 10) %>% 
  select(-n_genes) %>% 
  unnest(data) %>% 
  split(.$cell.type) %>% 
  lapply(function(x) unique(x$official.gene.symbol)) -> panglaoDB_list_BrainMarkers


panglaoDB_list_AllMarkers


res <- fgsea(
  pathways = panglaoDB_list_BrainMarkers,   # lista sygnatur
  stats    = cluster_markers_test$cluster_0,      # wektor nazwany
  nperm    = 1000        # liczba permutacji (więcej = dokładniej)
)
res <- fgsea(
  pathways = panglaoDB_list_BrainMarkers,   # lista sygnatur
  stats    = cluster_markers_test$cluster_7,      # wektor nazwany
  nperm    = 1000        # liczba permutacji (więcej = dokładniej)
)


res <- fgsea(
  pathways = panglaoDB_list_BrainMarkers,   # lista sygnatur
  stats    = cluster_markers_test$cluster_0,      # wektor nazwany
  nperm    = 1000        # liczba permutacji (więcej = dokładniej)
)




res %>% 
  arrange(desc(NES))

