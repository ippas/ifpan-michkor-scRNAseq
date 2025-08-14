run_enrichr_for_all_clusters <- function(markers_list, database, gene_col = "gene_symbol") {
  res_list <- lapply(
    X = names(markers_list),
    FUN = function(cl) {
      res <- run_enrichr(
        gene_list = markers_list[[cl]][[gene_col]],
        database  = database
      )
      dplyr::mutate(res, cluster = cl)
    }
  )
  names(res_list) <- names(markers_list)
  return(res_list)
}
