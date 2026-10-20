plot_detected_transcripts_per_cluster <- function(
    seurat_obj,
    cluster_var      = "seurat_clusters",
    assay            = "RNA",
    slot             = "counts",     
    features         = NULL,         
    min_counts       = 1,            
    min_cells_prop   = 0.05,         
    order_desc       = TRUE,         
    bar_fill         = "grey70",     
    text_size_axis   = 10,
    text_size_title  = 12
) {
  stopifnot(cluster_var %in% colnames(seurat_obj@meta.data))
  mat <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = slot)
  if (!is.null(features)) {
    keep <- intersect(rownames(mat), features)
    if (length(keep) == 0) stop("Żaden z podanych 'features' nie występuje w macierzy.")
    mat <- mat[keep, , drop = FALSE]
  }
  meta <- seurat_obj@meta.data
  common_cells <- intersect(colnames(mat), rownames(meta))
  mat  <- mat[, common_cells, drop = FALSE]
  meta <- meta[common_cells, , drop = FALSE]
  
  cl <- factor(meta[[cluster_var]])
  clusters <- levels(cl)
  
  count_detected_in_cluster <- function(cols_idx) {
    x <- mat[, cols_idx, drop = FALSE]
    detected_cells_per_feature <- Matrix::rowSums(x >= min_counts)
    thr <- ceiling(min_cells_prop * ncol(x))
    sum(detected_cells_per_feature >= thr)
  }
  
  n_detected <- vapply(clusters, function(k) {
    cols_idx <- which(cl == k)
    count_detected_in_cluster(cols_idx)
  }, FUN.VALUE = numeric(1))
  
  res_df <- data.frame(
    cluster = clusters,
    n_detected_transcripts = as.integer(n_detected),
    n_cells = as.integer(vapply(clusters, function(k) sum(cl == k), 0)),
    stringsAsFactors = FALSE
  )
  
  if (order_desc) {
    res_df <- res_df[order(res_df$n_detected_transcripts, decreasing = TRUE), ]
    res_df$cluster <- factor(res_df$cluster, levels = res_df$cluster)
  } else {
    res_df$cluster <- factor(res_df$cluster, levels = clusters)
  }
  
  p <- ggplot2::ggplot(res_df, ggplot2::aes(x = cluster, y = n_detected_transcripts)) +
    ggplot2::geom_col(fill = bar_fill) +
    ggplot2::geom_text(ggplot2::aes(label = n_detected_transcripts), vjust = -0.4, size = 3) +
    ggplot2::labs(
      x = "Cluster",
      y = "Detected transcripts",
      title = sprintf("Detected transcripts per cluster (≥%d counts in ≥%.0f%% cells)",
                      min_counts, 100 * min_cells_prop)
    ) +
    ggplot2::theme_classic(base_size = 12)
  
  print(p)
  invisible(res_df)
}


plot_detected_transcripts_per_cluster <- function(
    seurat_obj,
    cluster_var    = "seurat_clusters",
    assay          = "RNA",
    slot           = "counts",
    min_counts     = 1,
    min_cells_prop = 0.05,
    use_seurat_colors = TRUE,   # <- nowe
    order_numeric     = TRUE    # <- nowe
) {
  stopifnot(cluster_var %in% colnames(seurat_obj@meta.data))
  
  mat  <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = slot)
  meta <- seurat_obj@meta.data
  common_cells <- intersect(colnames(mat), rownames(meta))
  mat  <- mat[, common_cells, drop = FALSE]
  meta <- meta[common_cells, , drop = FALSE]
  
  cl <- factor(meta[[cluster_var]])
  clusters <- levels(cl)
  
  # kolejność numeryczna
  if (order_numeric) {
    clusters <- sort(as.numeric(as.character(clusters)))
    cl <- factor(cl, levels = clusters)
  }
  
  # liczba wykrytych transkryptów w klastrze
  count_detected_in_cluster <- function(cols_idx) {
    x <- mat[, cols_idx, drop = FALSE]
    detected_cells <- Matrix::rowSums(x >= min_counts)
    thr <- ceiling(min_cells_prop * ncol(x))
    sum(detected_cells >= thr)
  }
  
  n_detected <- vapply(levels(cl), function(k) {
    cols_idx <- which(cl == k)
    count_detected_in_cluster(cols_idx)
  }, FUN.VALUE = numeric(1))
  
  res_df <- data.frame(
    cluster = factor(levels(cl), levels = levels(cl)),
    n_detected_transcripts = as.integer(n_detected)
  )
  
  # kolory jak w Seurat
  if (use_seurat_colors) {
    n_clusters <- length(levels(res_df$cluster))
    cluster_colors <- Seurat:::DiscretePalette(n_clusters, palette = "alphabet")
    names(cluster_colors) <- levels(res_df$cluster)
    fill_scale <- ggplot2::scale_fill_manual(values = cluster_colors)
  } else {
    fill_scale <- ggplot2::scale_fill_manual(values = "grey70")
  }
  
  p <- ggplot2::ggplot(res_df, ggplot2::aes(x = cluster, y = n_detected_transcripts, fill = cluster)) +
    ggplot2::geom_col() +
    ggplot2::geom_text(ggplot2::aes(label = n_detected_transcripts), vjust = -0.4, size = 3) +
    ggplot2::labs(x = "Cluster", y = "Detected transcripts") +
    fill_scale +
    ggplot2::theme_classic()
  
  print(p)
  invisible(res_df)
}

plot_detected_transcripts_per_cluster <- function(
    seurat_obj,
    cluster_var    = "seurat_clusters",
    assay          = "RNA",
    slot           = "counts",
    min_counts     = 1,
    min_cells_prop = 0.05,
    order_numeric  = TRUE,
    match_seurat_colors = TRUE
) {
  mat  <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = slot)
  meta <- seurat_obj@meta.data
  common_cells <- intersect(colnames(mat), rownames(meta))
  mat  <- mat[, common_cells, drop = FALSE]
  meta <- meta[common_cells, , drop = FALSE]
  
  cl <- factor(meta[[cluster_var]])
  if (order_numeric) {
    cl_levels <- sort(as.numeric(as.character(levels(cl))))
    cl <- factor(cl, levels = cl_levels)
  }
  
  count_detected_in_cluster <- function(cols_idx) {
    x <- mat[, cols_idx, drop = FALSE]
    detected_cells <- Matrix::rowSums(x >= min_counts)
    thr <- ceiling(min_cells_prop * ncol(x))
    sum(detected_cells >= thr)
  }
  
  n_detected <- vapply(levels(cl), function(k) {
    cols_idx <- which(cl == k)
    count_detected_in_cluster(cols_idx)
  }, FUN.VALUE = numeric(1))
  
  res_df <- data.frame(
    cluster = factor(levels(cl), levels = levels(cl)),
    n_detected_transcripts = as.integer(n_detected)
  )
  
  p <- ggplot2::ggplot(res_df, ggplot2::aes(x = cluster, y = n_detected_transcripts, fill = cluster)) +
    ggplot2::geom_col() +
    ggplot2::geom_text(ggplot2::aes(label = n_detected_transcripts), vjust = -0.4, size = 3) +
    ggplot2::labs(x = "Cluster", y = "Detected transcripts") +
    ggplot2::theme_classic()
  
  if (match_seurat_colors) {
    # używa tego samego generatora co DimPlot
    p <- p + ggplot2::scale_fill_hue()
  }
  
  print(p)
  invisible(res_df)
}
