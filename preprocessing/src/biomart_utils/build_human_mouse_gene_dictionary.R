build_human_mouse_gene_dictionary <- function(
    ensembl_version,
    save_tsv = NULL,
    save_rds = NULL
) {
  # --- create mart ---
  human_mart <- biomaRt::useEnsembl(
    biomart = "ensembl",
    dataset = "hsapiens_gene_ensembl",
    version = ensembl_version
  )
  
  # --- get human gene info ---
  gene_info <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
    mart = human_mart
  )
  
  # --- get human → mouse ortholog info ---
  homolog_info <- biomaRt::getBM(
    attributes = c(
      "ensembl_gene_id",
      "mmusculus_homolog_ensembl_gene",
      "mmusculus_homolog_associated_gene_name",
      "mmusculus_homolog_orthology_type",
      "mmusculus_homolog_orthology_confidence"
    ),
    mart = human_mart
  )
  
  # --- merge ---
  out <- merge(gene_info, homolog_info, by = "ensembl_gene_id", all.x = TRUE)
  
  # --- rename ---
  colnames(out) <- c(
    "human_ensembl", "human_symbol", "human_biotype",
    "mouse_ensembl", "mouse_symbol",
    "orthology_type", "orthology_confidence"
  )
  
  # --- save if requested ---
  if (!is.null(save_tsv)) {
    utils::write.table(out, file = save_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  }
  if (!is.null(save_rds)) {
    saveRDS(out, file = save_rds)
  }
  
  out
}
