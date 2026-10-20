raw_counts_per_cluster <- extract_raw_by_cluster_cells(
  seurat_obj = integrated_pilotDexCortCtrl,
  cluster_col = "seurat_clusters_res0.13",
  assay      = "RNA",
  layer      = "counts"
)


raw_counts_per_cluster$cluster_0 %>% dim

boxplot_geneSymbolExpression_by_cluster(
  data_list = raw_counts_per_cluster,          # lista z extract_raw_by_cluster_cells()
  gene_symbol = "SYN3",
  reorder_by_median = F,
  point_size = 0.2
)


boxplot_multiGeneSymbols_by_cluster(
  data_list    = agg_list_mean,
  gene_symbols = c("S100B","GFAP","AQP4","ALDH1L1", "SYN3", "PDGFRA"),
  n_rows = 3,                 # układ 2 wiersze
  collect_guides = F,
  point_size = 0.2
)


boxplot_multiGeneSymbols_by_cluster(
  data_list    = agg_list_mean,
  gene_symbols = c("VWA3A", "POSTN", "TTC6", "MOG", "GALNTL6", "RARB", "RGS3", "LINC02742"),
  n_rows = 3,                 # układ 2 wiersze
  collect_guides = F,
  point_size = 0.2
)

boxplot_multiGeneSymbols_by_cluster(
  data_list    = agg_list_mean,
  gene_symbols = c("MOG", "RARB", "NEK2", "TTC6", "NTN4", "POSTN", "DNAAF1"),
  n_rows = 3,                 # układ 2 wiersze
  collect_guides = F,
  point_size = 0.2
)



tmp <- select_high_abundant_low_discrete(
  raw_counts_list = raw_counts_per_cluster,
  top_n = 100,
  min_pct = 0.7,
  max_gini = 0.5,
  min_mean = 0.,
  return_long = TRUE
)

markers_specific <- find_specific_markers(
  raw_counts_list = raw_counts_per_cluster,
  min_mean        = 0.1,
  min_log2fc      = 0.5,
  min_tau         = 0.5,
  verbose         = TRUE
)


markers_specific %>%  
  separate(gene, into = c("ensembl_id", "gene_symbol"), sep = "-", extra = "merge") %>% 
  filter(cluster == "cluster_8")
