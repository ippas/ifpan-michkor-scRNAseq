#' Prepare gene × cluster average expression using your sample-then-cluster pipeline
#'
#' Replicates your original steps:
#'  1) Ensembl v114 protein-coding filter (chr 1:22,X,Y; non-empty symbol)
#'  2) Parse Ensembl IDs from Seurat feature names (assumes "ENSG-..."; trims version)
#'  3) aggregate_counts_by_cluster_sample(..., return="list")
#'  4) filter_list_by_ensemblID(..., return_only_geneSymbol=TRUE)
#'  5) rowMeans across samples per cluster → genes × clusters matrix (rownames = HGNC symbols)
#'
#' Requires your helpers: download_biomart_genes(), aggregate_counts_by_cluster_sample(),
#' filter_list_by_ensemblID().
prepare_avg_expr_by_sample_cluster <- function(
    seurat_obj,
    cluster_col   = "seurat_clusters_res0.13",
    sample_col    = "sample_id",
    assay         = "RNA",
    layer         = "counts",     # your code used "counts"
    agg           = "mean",
    biomart_version = "114",
    chromosomes   = c(as.character(1:22), "X", "Y"),
    verbose       = TRUE
) {
  # --- sanity & deps ---
  stopifnot(requireNamespace("dplyr", quietly = TRUE))
  stopifnot(requireNamespace("tidyr", quietly = TRUE))
  stopifnot(requireNamespace("Matrix", quietly = TRUE))
  if (!exists("download_biomart_genes")) {
    stop("download_biomart_genes() not found. Please source your GitHub script first.")
  }
  if (!exists("aggregate_counts_by_cluster_sample")) {
    stop("aggregate_counts_by_cluster_sample() not found (needed to match your pipeline).")
  }
  if (!exists("filter_list_by_ensemblID")) {
    stop("filter_list_by_ensemblID() not found (needed to match your pipeline).")
  }
  
  # --- 1) Ensembl v114 protein-coding table ---
  if (verbose) message("Downloading Ensembl v", biomart_version, " annotations…")
  genes_all <- download_biomart_genes(version = biomart_version)
  genes_protein <- genes_all %>%
    dplyr::filter(
      gene_biotype == "protein_coding",
      chromosome_name %in% chromosomes,
      external_gene_name != ""
    ) %>%
    dplyr::select(ensembl_gene_id, external_gene_name) %>%
    dplyr::distinct()
  
  # --- 2) Ensembl IDs from Seurat feature names ("ENSG-..."; strip version) ---
  rn <- rownames(seurat_obj[[assay]])
  ensembl_id <- sub("-.*$", "", rn)          # "ENSG-TP53" → "ENSG…"
  ensembl_id <- sub("\\.\\d+$", "", ensembl_id)  # remove ".xx" versions
  
  genes_df_seurat <- dplyr::tibble(name_id = rn, ensembl_id = ensembl_id)
  
  # Keep only those features that are protein-coding and have a symbol in v114
  genes_df_pc <- genes_df_seurat %>%
    dplyr::left_join(genes_protein, by = c("ensembl_id" = "ensembl_gene_id")) %>%
    tidyr::drop_na(external_gene_name) %>%
    dplyr::rename(gene_symbol = external_gene_name)
  
  if (!nrow(genes_df_pc)) {
    stop("No protein-coding features matched between Seurat and Ensembl v114.")
  }
  
  # --- 3) Aggregate per sample within clusters (your helper) ---
  if (verbose) message("Aggregating counts by cluster/sample using your helper…")
  agg_list <- aggregate_counts_by_cluster_sample(
    seurat_obj  = seurat_obj,
    cluster_col = cluster_col,
    sample_col  = sample_col,
    assay       = assay,
    layer       = layer,
    agg         = agg,
    return      = "list",
    progress    = verbose
  )
  
  # --- 4) Filter list to protein-coding Ensembl IDs and return gene symbols ---
  agg_list <- filter_list_by_ensemblID(
    gene_expression_list   = agg_list,
    ensembl_ids            = genes_df_pc$ensembl_id,
    return_only_geneSymbol = TRUE
  )
  # At this point, each element is a matrix/data.frame of genes (symbols) × samples for a given cluster
  
  # --- 5) Mean across samples → one column per cluster ---
  if (verbose) message("Collapsing to mean per cluster (rowMeans over samples)…")
  avg_mat <- sapply(agg_list, Matrix::rowMeans)
  
  # Safety: ensure numeric and finite
  avg_mat <- as.matrix(avg_mat)
  mode(avg_mat) <- "numeric"
  avg_mat[!is.finite(avg_mat)] <- 0
  
  if (verbose) {
    message("Done: ", nrow(avg_mat), " genes × ", ncol(avg_mat), " clusters.")
  }
  avg_mat  # genes (symbols) × clusters (as in your names of agg_list)
}


prepare_avg_expr_by_sample_cluster <- function(
    seurat_obj,
    cluster_col   = "seurat_clusters_res0.13",
    sample_col    = "sample_id",
    assay         = "RNA",
    layer         = "counts",     # your code used "counts"
    agg           = "mean",
    biomart_version = "114",
    chromosomes   = c(as.character(1:22), "X", "Y"),
    select_samples = NULL,        # <--- NEW
    verbose       = TRUE
) {
  # --- sanity & deps ---
  stopifnot(requireNamespace("dplyr", quietly = TRUE))
  stopifnot(requireNamespace("tidyr", quietly = TRUE))
  stopifnot(requireNamespace("Matrix", quietly = TRUE))
  if (!exists("download_biomart_genes")) {
    stop("download_biomart_genes() not found. Please source your GitHub script first.")
  }
  if (!exists("aggregate_counts_by_cluster_sample")) {
    stop("aggregate_counts_by_cluster_sample() not found (needed to match your pipeline).")
  }
  if (!exists("filter_list_by_ensemblID")) {
    stop("filter_list_by_ensemblID() not found (needed to match your pipeline).")
  }
  
  # --- 1) Ensembl v114 protein-coding table ---
  if (verbose) message("Downloading Ensembl v", biomart_version, " annotations…")
  genes_all <- download_biomart_genes(version = biomart_version)
  genes_protein <- genes_all %>%
    dplyr::filter(
      gene_biotype == "protein_coding",
      chromosome_name %in% chromosomes,
      external_gene_name != ""
    ) %>%
    dplyr::select(ensembl_gene_id, external_gene_name) %>%
    dplyr::distinct()
  
  # --- 2) Ensembl IDs from Seurat feature names ("ENSG-..."; strip version) ---
  rn <- rownames(seurat_obj[[assay]])
  ensembl_id <- sub("-.*$", "", rn)          # "ENSG-TP53" → "ENSG…"
  ensembl_id <- sub("\\.\\d+$", "", ensembl_id)  # remove ".xx" versions
  
  genes_df_seurat <- dplyr::tibble(name_id = rn, ensembl_id = ensembl_id)
  
  # Keep only those features that are protein-coding and have a symbol in v114
  genes_df_pc <- genes_df_seurat %>%
    dplyr::left_join(genes_protein, by = c("ensembl_id" = "ensembl_gene_id")) %>%
    tidyr::drop_na(external_gene_name) %>%
    dplyr::rename(gene_symbol = external_gene_name)
  
  if (!nrow(genes_df_pc)) {
    stop("No protein-coding features matched between Seurat and Ensembl v114.")
  }
  
  # --- 3) Aggregate per sample within clusters (your helper) ---
  if (verbose) message("Aggregating counts by cluster/sample using your helper…")
  agg_list <- aggregate_counts_by_cluster_sample(
    seurat_obj  = seurat_obj,
    cluster_col = cluster_col,
    sample_col  = sample_col,
    assay       = assay,
    layer       = layer,
    agg         = agg,
    return      = "list",
    progress    = verbose
  )
  
  # --- 3a) Optionally subset to select_samples ---
  if (!is.null(select_samples)) {
    if (verbose) message("Subsetting to selected samples…")
    agg_list <- lapply(agg_list, function(m) {
      keep <- intersect(colnames(m), select_samples)
      if (length(keep) == 0) {
        warning("No matching samples found in a cluster; returning empty matrix.")
        return(m[, FALSE, drop = FALSE])
      }
      m[, keep, drop = FALSE]
    })
  }
  
  # --- 4) Filter list to protein-coding Ensembl IDs and return gene symbols ---
  agg_list <- filter_list_by_ensemblID(
    gene_expression_list   = agg_list,
    ensembl_ids            = genes_df_pc$ensembl_id,
    return_only_geneSymbol = TRUE
  )
  
  # --- 5) Mean across samples → one column per cluster ---
  if (verbose) message("Collapsing to mean per cluster (rowMeans over samples)…")
  avg_mat <- sapply(agg_list, Matrix::rowMeans)
  
  # Safety: ensure numeric and finite
  avg_mat <- as.matrix(avg_mat)
  mode(avg_mat) <- "numeric"
  avg_mat[!is.finite(avg_mat)] <- 0
  
  if (verbose) {
    message("Done: ", nrow(avg_mat), " genes × ", ncol(avg_mat), " clusters.")
  }
  avg_mat  # genes (symbols) × clusters
}
