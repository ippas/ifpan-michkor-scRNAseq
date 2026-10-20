#' Boxplot of gene expression by cluster
#'
#' @param data_list List of cluster matrices (from aggregate_counts_by_cluster_sample or extract_raw_by_cluster_cells)
#' @param gene Gene symbol/ID (must be present in rownames of matrices)
#' @param cluster_colors Optional named vector of colors per cluster
#' @param reorder_by_median Logical; reorder clusters by median expression
#' @param cluster_order Character vector specifying cluster order (e.g., paste0("cluster_", 0:9))
#' @param title Plot title
#' @param ylab Y-axis label
#' @param point_size Size of jittered points
#' @param point_alpha Alpha of jittered points
#' @param box_linewidth Line width of boxplot
#'
#' @return A ggplot object
#'
boxplot_geneExpression_by_cluster <- function(
    data_list,
    gene,
    cluster_colors   = NULL,
    reorder_by_median = FALSE,
    cluster_order    = NULL,
    title            = NULL,
    ylab             = "Expression per cluster",
    point_size       = 1.8,
    point_alpha      = 0.7,
    box_linewidth    = 0.6
){
  stopifnot(is.list(data_list), length(data_list) > 0)
  
  # Collect per-cell values
  df_list <- lapply(names(data_list), function(cl){
    mat <- data_list[[cl]]
    rn <- rownames(mat)
    if (is.null(rn) || !gene %in% rn) return(NULL)
    vals <- as.numeric(mat[gene, , drop = TRUE])
    tibble::tibble(
      cluster = cl,
      sample  = colnames(mat),
      value   = vals
    )
  })
  df <- dplyr::bind_rows(df_list)
  if (is.null(df) || nrow(df) == 0) stop(sprintf("Gene '%s' not found in any cluster.", gene))
  
  # Ordering clusters
  cl_order <- .merge_cluster_order(df$cluster, cluster_order)
  df$cluster <- factor(df$cluster, levels = cl_order)
  
  if (reorder_by_median && is.null(cluster_order)) {
    ord <- df |>
      dplyr::group_by(cluster) |>
      dplyr::summarise(median = median(value, na.rm = TRUE), .groups="drop") |>
      dplyr::arrange(dplyr::desc(median)) |>
      dplyr::pull(cluster) |> as.character()
    df$cluster <- factor(df$cluster, levels = ord)
  }
  
  if (is.null(title)) title <- paste0("Boxplot: ", gene)
  
  # Plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = cluster, y = value, color = cluster)) +
    ggplot2::geom_boxplot(width = 0.6, linewidth = box_linewidth, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.15, height = 0, size = point_size, alpha = point_alpha) +
    ggplot2::labs(title = title, x = NULL, y = ylab) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.major.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 10)
    )
  
  if (!is.null(cluster_colors)) {
    p <- p + ggplot2::scale_color_manual(values = cluster_colors)
  }
  
  return(p)
}


boxplot_geneSymbolExpression_by_cluster <- function(
    data_list,
    gene_symbol,
    ...
){
  stopifnot(is.list(data_list), length(data_list) > 0)
  
  # spróbuj z pierwszej macierzy
  rn <- rownames(data_list[[1]])
  if (is.null(rn)) stop("Row names are NULL in the first cluster matrix.")
  
  match_idx <- grep(paste0("-", gene_symbol, "$"), rn, value = TRUE)
  
  # fallback: przeszukaj wszystkie klastry, gdyby pierwszy nie miał tego wiersza
  if (length(match_idx) == 0) {
    all_rn <- unique(unlist(lapply(data_list, rownames)))
    match_idx <- grep(paste0("-", gene_symbol, "$"), all_rn, value = TRUE)
  }
  
  if (length(match_idx) == 0) {
    stop(sprintf("Gene symbol '%s' was not found.", gene_symbol))
  }
  if (length(match_idx) > 1) {
    warning(sprintf("Multiple matches found for symbol '%s': using the first (%s)",
                    gene_symbol, match_idx[1]))
  }
  
  gene_id <- match_idx[1]
  
  boxplot_geneExpression_by_cluster(
    data_list = data_list,
    gene      = gene_id,
    ...
  )
}

boxplot_multiGeneSymbols_by_cluster <- function(
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
      p <- boxplot_geneSymbolExpression_by_cluster(
        data_list   = data_list,
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

