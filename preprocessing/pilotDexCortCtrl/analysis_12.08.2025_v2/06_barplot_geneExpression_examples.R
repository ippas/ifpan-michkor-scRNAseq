agg_list_mean <- aggregate_counts_by_cluster_sample(
  seurat_obj  = integrated_pilotDexCortCtrl,
  cluster_col = "seurat_clusters_res0.13",
  sample_col  = "sample_id",
  assay       = "RNA",
  layer       = "counts",  # <- Seurat v5
  agg         = "mean",
  return      = "list",
  progress    = TRUE
)


barplot_geneExpression_by_cluster(
  data_list = agg_list_mean,
  gene = "ENSG00000160307-S100B",
  reorder_by_mean = TRUE
)


barplot_geneSymbolExpression_by_cluster(
  data_list = agg_list_mean,
  gene_symbol = "COL15A1",
  reorder_by_mean = F
)
