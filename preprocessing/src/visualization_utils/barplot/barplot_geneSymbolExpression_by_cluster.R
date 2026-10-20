# =====================================================================
# Function: barplot_geneSymbolExpression_by_cluster
# ---------------------------------------------------------------------
# Wrapper function that allows specifying only the gene symbol (e.g. "S100B"),
# automatically finds the full ENSG-symbol ID in the dataset, 
# and passes it to barplot_geneExpression_by_cluster.
# =====================================================================
barplot_geneSymbolExpression_by_cluster <- function(
    data_list,
    gene_symbol,
    ...
){
  # Extract row names from the first cluster matrix
  rn <- rownames(data_list[[1]])
  
  # Match ENSG IDs ending with the requested symbol
  match_idx <- grep(paste0("-", gene_symbol, "$"), rn, value = TRUE)
  
  if (length(match_idx) == 0) {
    stop(sprintf("Gene symbol '%s' was not found.", gene_symbol))
  }
  if (length(match_idx) > 1) {
    warning(sprintf("Multiple matches found for symbol '%s': using the first (%s)", 
                    gene_symbol, match_idx[1]))
  }
  
  gene_id <- match_idx[1]
  
  # Call the main function
  barplot_geneExpression_by_cluster(
    data_list = data_list,
    gene = gene_id,
    ...
  )
}

barplot_multiGeneSymbols_by_cluster <- function(
    data_list,
    gene_symbols,
    n_rows = NULL,
    n_cols = NULL,
    collect_guides = TRUE,
    title = NULL,
    ...
){
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required. Install with install.packages('patchwork').")
  }
  
  # Build individual plots; skip genes that fail matching
  plots <- lapply(gene_symbols, function(sym){
    tryCatch({
      p <- barplot_geneSymbolExpression_by_cluster(
        data_list = data_list,
        gene_symbol = sym,
        ...
      )
      p + ggplot2::ggtitle(sym)
    }, error = function(e){
      warning(sprintf("Skipping '%s': %s", sym, e$message))
      NULL
    })
  })
  plots <- Filter(Negate(is.null), plots)
  k <- length(plots)
  if (k == 0) stop("No valid genes to plot after matching.")
  
  # Decide grid layout
  if (is.null(n_rows) && is.null(n_cols)) {
    n_cols <- ceiling(sqrt(k))
    n_rows <- ceiling(k / n_cols)
  } else if (is.null(n_cols)) {
    n_cols <- ceiling(k / n_rows)
  } else if (is.null(n_rows)) {
    n_rows <- ceiling(k / n_cols)
  }
  
  # Combine with patchwork
  combined <- Reduce(`+`, plots)
  combined <- combined + patchwork::plot_layout(
    ncol = n_cols,
    nrow = n_rows,
    guides = if (collect_guides) "collect" else "auto"
  )
  
  if (!is.null(title)) {
    combined <- combined + patchwork::plot_annotation(title = title)
  }
  
  return(combined)
}

