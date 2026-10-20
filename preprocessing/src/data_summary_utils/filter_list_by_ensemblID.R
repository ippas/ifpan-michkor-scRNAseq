#' Filter gene_expression_list by Ensembl IDs
#'
#' @param gene_expression_list lista obiektów dgeMatrix (po jednym na klaster)
#' @param ensembl_ids wektor Ensembl ID (np. "ENSG00000141510")
#' @param return_only_geneSymbol logic, czy zwrócić tylko symbole genów w rownames
#' @param make_unique logic, czy wymuszać unikalność nazw genów po zamianie
#' @param drop_empty logic, czy usuwać klastry bez trafień
#' @param keep_order logic, czy zachować kolejność ensembl_ids w wynikach
#'
#' @return lista przefiltrowanych macierzy ekspresji genów
filter_list_by_ensemblID <- function(
    gene_expression_list,
    ensembl_ids,
    return_only_geneSymbol = FALSE,
    make_unique = TRUE,
    drop_empty = TRUE,
    keep_order = TRUE
){
  stopifnot(is.list(gene_expression_list))
  ens_set <- unique(as.character(ensembl_ids))
  
  extract_parts <- function(rn){
    ensg   <- sub("-.*$", "", rn)          # część przed "-"
    symbol <- sub("^[^-]+-", "", rn)       # część po "-"
    list(ensg = ensg, symbol = symbol)
  }
  
  out <- lapply(gene_expression_list, function(mat){
    rn <- rownames(mat)
    parts <- extract_parts(rn)
    
    # trafienia po Ensembl ID
    hits <- parts$ensg %in% ens_set
    if (!any(hits)) {
      if (drop_empty) return(NULL)
      return(mat[0, , drop = FALSE]) # pusta macierz
    }
    
    mat_sub <- mat[hits, , drop = FALSE]
    
    # kolejność wg ensembl_ids
    if (keep_order) {
      ord <- match(ens_set, parts$ensg[hits])
      ord <- ord[!is.na(ord)]
      if (length(ord)) mat_sub <- mat_sub[ord, , drop = FALSE]
    }
    
    # tylko symbole genów
    if (return_only_geneSymbol) {
      sym <- sub("^[^-]+-", "", rownames(mat_sub))
      if (make_unique) sym <- make.unique(sym, sep = ".")
      rownames(mat_sub) <- sym
    }
    
    mat_sub
  })
  
  if (drop_empty) {
    keep <- !vapply(out, is.null, logical(1))
    out <- out[keep]
  }
  out
}
