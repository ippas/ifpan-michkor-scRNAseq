prepare_proteincoding_avg_expr <- function(
    seurat_obj,
    cluster_col = "seurat_clusters_res0.13",
    sample_col  = "sample_id",
    assay       = "RNA",
    layer       = "counts",
    biomart_version = "114",
    select_samples  = NULL,
    verbose = TRUE
) {
  if (verbose) message("Step 1: Downloading genes from Ensembl v", biomart_version, "...")
  genes_all <- download_biomart_genes(version = biomart_version)
  
  genes_proteinCoding <- genes_all %>%
    dplyr::filter(
      gene_biotype == "protein_coding",
      chromosome_name %in% c(1:22, "X", "Y"),
      external_gene_name != ""
    )
  
  if (verbose) message("Step 2: Extracting Ensembl IDs from Seurat object...")
  genes_df_seuratObj <- data.frame(
    name_id = rownames(seurat_obj[[assay]]),
    stringsAsFactors = FALSE
  ) %>%
    tidyr::separate(
      name_id, into = c("ensembl_id", "gene_symbol"),
      sep = "-", extra = "drop"
    ) %>%
    dplyr::select(ensembl_id)
  
  genes_df_proteinCoding_seuratObj <- genes_df_seuratObj %>%
    dplyr::left_join(genes_proteinCoding,
                     by = c("ensembl_id" = "ensembl_gene_id")) %>%
    tidyr::drop_na() %>%
    dplyr::rename(gene_symbol = external_gene_name)
  
  if (verbose) message("Step 3: Aggregating mean expression per cluster and sample...")
  agg_list_mean <- aggregate_counts_by_cluster_sample(
    seurat_obj  = seurat_obj,
    cluster_col = cluster_col,
    sample_col  = sample_col,
    assay       = assay,
    layer       = layer,
    agg         = "mean",
    return      = "list",
    select_samples = select_samples,
    progress    = verbose
  )
  
  if (verbose) message("Step 4: Filtering expression list by Ensembl IDs...")
  agg_list_mean <- filter_list_by_ensemblID(
    gene_expression_list   = agg_list_mean,
    ensembl_ids            = genes_df_proteinCoding_seuratObj$ensembl_id,
    return_only_geneSymbol = TRUE
  )
  
  if (verbose) message("Step 5: Computing mean expression across clusters...")
  agg_mean_expr <- sapply(agg_list_mean, Matrix::rowMeans)
  
  return(agg_mean_expr)
}


