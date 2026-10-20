# ---------- helpers ----------
# Natural numeric order for labels like "cluster_0","cluster_1",...
.cluster_levels_numeric <- function(x) {
  x <- unique(as.character(x))
  nums <- suppressWarnings(as.numeric(sub("^cluster_", "", x)))
  if (any(is.na(nums))) return(sort(x))
  x[order(nums)]
}

# Merge desired order with observed clusters; append missing in numeric order
.merge_cluster_order <- function(observed, desired = NULL) {
  observed <- unique(as.character(observed))
  if (is.null(desired)) return(.cluster_levels_numeric(observed))
  
  desired <- unique(as.character(desired))
  extra <- setdiff(desired, observed)
  if (length(extra) > 0) {
    warning("Clusters in `cluster_order` not found and will be ignored: ",
            paste(extra, collapse = ", "))
  }
  keep_desired <- intersect(desired, observed)
  missing <- setdiff(observed, keep_desired)
  c(keep_desired, .cluster_levels_numeric(missing))
}


# =====================================================================
# Function: barplot_geneExpression_by_cluster
# ---------------------------------------------------------------------
# Barplot for a given gene (full ENSG-symbol ID).
# - Bars = cluster means (transparent fill, colored outline)
# - Error bars = SE
# - Dots = per-sample values
# - Cluster order: by `cluster_order` if provided; else numeric; unless
#   `reorder_by_mean = TRUE` AND no `cluster_order` provided.
# =====================================================================
barplot_geneExpression_by_cluster <- function(
    data_list,
    gene,
    cluster_colors = NULL,
    reorder_by_mean = FALSE,      # ignored if cluster_order is provided
    cluster_order = NULL,         # character vector, e.g., paste0("cluster_", 0:9)
    title = NULL,
    ylab = "Mean expression per cluster",
    point_size = 1.9,
    point_alpha = 0.9,
    bar_width = 0.7,
    bar_linewidth = 0.6
){
  stopifnot(is.list(data_list), length(data_list) > 0)
  
  # Collect sample-level values
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
  if (is.null(df) || nrow(df) == 0) stop(sprintf("Gene '%s' was not found in any cluster.", gene))
  
  # Summaries
  sum_df <- df |>
    dplyr::group_by(cluster) |>
    dplyr::summarise(
      n    = sum(!is.na(value)),
      mean = mean(value, na.rm = TRUE),
      sd   = stats::sd(value, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(se = ifelse(n > 1, sd / sqrt(n), NA_real_))
  
  # Ordering
  cl_order <- .merge_cluster_order(df$cluster, cluster_order)
  df$cluster     <- factor(df$cluster, levels = cl_order)
  sum_df$cluster <- factor(sum_df$cluster, levels = cl_order)
  
  # If user asked to reorder by mean AND didn't supply cluster_order, do it
  if (reorder_by_mean && is.null(cluster_order)) {
    ord <- sum_df |>
      dplyr::arrange(dplyr::desc(mean)) |>
      dplyr::pull(cluster) |> as.character()
    df$cluster     <- factor(df$cluster, levels = ord)
    sum_df$cluster <- factor(sum_df$cluster, levels = ord)
  }
  
  if (is.null(title)) title <- paste0("Barplot: ", gene)
  
  # Plot
  p <- ggplot2::ggplot() +
    ggplot2::geom_col(
      data = sum_df,
      ggplot2::aes(x = cluster, y = mean, color = cluster),
      fill = NA, width = bar_width, linewidth = bar_linewidth
    ) +
    ggplot2::geom_errorbar(
      data = sum_df,
      ggplot2::aes(x = cluster, ymin = mean - se, ymax = mean + se, color = cluster),
      width = 0.15, linewidth = 0.4
    ) +
    ggplot2::geom_jitter(
      data = df,
      ggplot2::aes(x = cluster, y = value, color = cluster),
      width = 0.15, height = 0, size = point_size, alpha = point_alpha
    ) +
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
  p
}

