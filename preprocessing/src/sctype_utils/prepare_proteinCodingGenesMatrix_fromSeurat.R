prepare_proteinCodingGenesMatrix_fromSeurat <- function(
    seurat_obj,
    assay = "RNA",
    layer = "data",
    biomart_version = "v112",
    method = c("sum", "mean")
) {
  method <- match.arg(method)
  
  message("Step 1: Extracting expression matrix from Seurat object...")
  exprMat <- as.matrix(LayerData(seurat_obj[[assay]], layer = layer))
  
  message("Step 2: Normalizing rownames (extracting pure Ensembl IDs)...")
  ensembl_ids <- sub("-.*", "", rownames(exprMat))
  rownames(exprMat) <- ensembl_ids
  
  message("Step 3: Downloading gene annotations from Ensembl (BioMart version: ", biomart_version, ")...")
  bm <- download_biomart_genes(version = biomart_version)
  
  symbol_col <- if ("hgnc_symbol" %in% colnames(bm)) {
    "hgnc_symbol"
  } else if ("external_gene_name" %in% colnames(bm)) {
    "external_gene_name"
  } else {
    stop("No gene symbol column found (expected 'hgnc_symbol' or 'external_gene_name').")
  }
  
  message("Step 4: Filtering for protein-coding genes on autosomes and sex chromosomes...")
  gene_info <- bm %>%
    dplyr::filter(
      chromosome_name %in% c(as.character(1:22), "X", "Y"),
      gene_biotype == "protein_coding"
    ) %>%
    dplyr::select(
      ensembl_gene_id,
      hgnc_symbol = !!sym(symbol_col),
      gene_biotype,
      chromosome_name
    ) %>%
    dplyr::distinct()
  message("  ▶ Number of protein-coding genes in BioMart: ", nrow(gene_info))
  
  message("Step 5: Matching Ensembl IDs between expression matrix and BioMart annotation...")
  common_ids <- intersect(rownames(exprMat), gene_info$ensembl_gene_id)
  if (length(common_ids) == 0) {
    stop("No common Ensembl IDs found between exprMat and BioMart.")
  }
  message("  ▶ Number of genes matched with expression matrix: ", length(common_ids))
  
  exprMat <- exprMat[common_ids, , drop = FALSE]
  gene_info <- gene_info[match(common_ids, gene_info$ensembl_gene_id), ]
  
  message("Step 6: Replacing Ensembl IDs with gene symbols...")
  rownames(exprMat) <- gene_info$hgnc_symbol
  
  message("Step 7: Aggregating duplicate gene symbols using method = '", method, "' ...")
  if (method == "sum") {
    exprMat <- rowsum(exprMat, group = rownames(exprMat))
  } else if (method == "mean") {
    exprMat <- rowsum(exprMat, group = rownames(exprMat)) / as.vector(table(rownames(exprMat)))
  }
  message("  ▶ Number of unique genes after aggregation: ", nrow(exprMat))
  
  message("✅ Done! Final matrix has ", nrow(exprMat), " genes and ", ncol(exprMat), " cells/samples.")
  
  # clean up temporary variables to save memory
  rm(bm, gene_info, common_ids, ensembl_ids, symbol_col)
  gc()
  
  return(exprMat)
}
