pseudobulk_heatmap_global_complexHeatmap <- function(
    agg,
    sample_meta,
    sample_col   = "sample_id",
    group_col    = "treatment",
    group_levels = c("CTRL","DEX","CORT"),
    group_palette = c(CTRL="#1b9e77", DEX="#d95f02", CORT="#7570b3"),
    genes        = NULL,
    top_n_var    = 50,
    transform    = c("logCPM","log2","none"),
    log2_offset  = 1,
    scale_rows   = TRUE,
    z_score_limits = NULL,           # clamp z-scores, e.g. c(-4,4)
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    legend_title = "Expression (z)",
    show_colnames = TRUE,
    show_rownames = TRUE,
    verbose      = TRUE,
    return_matrix = TRUE,
    quantile_normalize = FALSE,
    heatmap_colors = c("blue","white","red")   # <--- NEW
) {
  .msg <- function(...) if (isTRUE(verbose)) message(...)
  
  # ------------------------------
  # Simple quantile normalization in base R
  # ------------------------------
  .qn_base <- function(m) {
    m <- as.matrix(m)
    n <- nrow(m); p <- ncol(m)
    ord <- apply(m, 2, order, na.last = TRUE)
    m_sorted <- m
    for (j in seq_len(p)) m_sorted[, j] <- m[ord[, j], j]
    avg <- rowMeans(m_sorted, na.rm = TRUE)
    m_qn <- matrix(NA_real_, n, p, dimnames = dimnames(m))
    for (j in seq_len(p)) {
      inv <- integer(n); inv[ord[, j]] <- seq_len(n)
      m_qn[, j] <- avg[inv]
    }
    m_qn
  }
  
  transform <- match.arg(transform)
  
  # ------------------------------
  # Helpers to merge matrices
  # ------------------------------
  .sum_by_rownames <- function(A, B) {
    A <- as.matrix(A); B <- as.matrix(B)
    all_genes <- union(rownames(A), rownames(B))
    common_samples <- intersect(colnames(A), colnames(B))
    if (!length(common_samples)) stop("No overlapping sample columns between matrices.")
    A2 <- matrix(0, nrow = length(all_genes), ncol = length(common_samples),
                 dimnames = list(all_genes, common_samples))
    B2 <- A2
    A2[rownames(A), intersect(colnames(A), common_samples)] <-
      A[, intersect(colnames(A), common_samples), drop = FALSE]
    B2[rownames(B), intersect(colnames(B), common_samples)] <-
      B[, intersect(colnames(B), common_samples), drop = FALSE]
    A2 + B2
  }
  
  .as_total_counts <- function(agg) {
    if (is.list(agg)) {
      mats <- Filter(function(m) (is.matrix(m) || inherits(m, "Matrix")) &&
                       nrow(m) > 0 && ncol(m) > 0, agg)
      if (!length(mats)) stop("Input list is empty or contains empty matrices.")
      .msg(sprintf("[heatmap] Merging %d cluster matrices into one global matrix...", length(mats)))
      mats <- lapply(mats, as.matrix)
      Reduce(.sum_by_rownames, mats)
    } else if (is.matrix(agg) || inherits(agg, "Matrix")) {
      as.matrix(agg)
    } else {
      stop("'agg' must be a matrix, sparse Matrix, or a list of such.")
    }
  }
  
  # ------------------------------
  # Assemble counts and align metadata
  # ------------------------------
  counts <- .as_total_counts(agg)
  stopifnot(all(c(sample_col, group_col) %in% colnames(sample_meta)))
  smeta <- as.data.frame(sample_meta)
  if (any(duplicated(smeta[[sample_col]]))) {
    .msg("[heatmap] Deduplicating sample metadata by ", sample_col)
    smeta <- smeta[!duplicated(smeta[[sample_col]]), , drop = FALSE]
  }
  keep_cols <- intersect(colnames(counts), smeta[[sample_col]])
  if (!length(keep_cols)) stop("No overlap between counts colnames and sample_meta[[sample_col]].")
  counts <- counts[, keep_cols, drop = FALSE]
  smeta  <- smeta[match(colnames(counts), smeta[[sample_col]]), , drop = FALSE]
  if (!is.null(group_levels)) smeta[[group_col]] <- factor(smeta[[group_col]], levels = group_levels)
  
  .msg(sprintf("[heatmap] Matrix: %d genes × %d samples", nrow(counts), ncol(counts)))
  .msg(sprintf("[heatmap] Groups: %s", paste(levels(smeta[[group_col]]), collapse = ", ")))
  
  # ------------------------------
  # Quantile normalization (optional)
  # ------------------------------
  if (isTRUE(quantile_normalize)) {
    .msg("[heatmap] Quantile normalization: custom base R implementation")
    counts <- .qn_base(counts)
  }
  
  # ------------------------------
  # Transform
  # ------------------------------
  mat <- counts
  if (transform == "logCPM") {
    .msg("[heatmap] Transform: edgeR::cpm(log=TRUE)")
    mat <- edgeR::cpm(mat, log = TRUE, prior.count = 1)
  } else if (transform == "log2") {
    .msg("[heatmap] Transform: log2(counts + offset), offset = ", log2_offset)
    mat <- log2(mat + log2_offset)
  } else {
    .msg("[heatmap] Transform: none (raw scale)")
  }
  
  # ------------------------------
  # Gene selection
  # ------------------------------
  if (!is.null(genes)) {
    genes <- intersect(genes, rownames(mat))
    if (!length(genes)) stop("Provided 'genes' not found in matrix rownames.")
    mat <- mat[genes, , drop = FALSE]
    .msg(sprintf("[heatmap] Using provided genes: %d", length(genes)))
  } else {
    v <- matrixStats::rowVars(mat)
    ord <- order(v, decreasing = TRUE)
    n <- min(top_n_var, nrow(mat))
    mat <- mat[ord[seq_len(n)], , drop = FALSE]
    .msg(sprintf("[heatmap] Using top %d most variable genes", n))
  }
  
  # ------------------------------
  # Row scaling (z-score) + clamping
  # ------------------------------
  if (isTRUE(scale_rows)) {
    .msg("[heatmap] Scaling rows (z-score)")
    row_means <- rowMeans(mat, na.rm = TRUE)
    row_sds   <- matrixStats::rowSds(mat, na.rm = TRUE)
    row_sds[row_sds == 0] <- 1
    mat <- sweep(mat, 1, row_means, "-")
    mat <- sweep(mat, 1, row_sds,   "/")
    legend_title <- "Expression (z)"
  }
  if (!is.null(z_score_limits)) {
    .msg(sprintf("[heatmap] Clamping z-scores to range [%.1f, %.1f]",
                 z_score_limits[1], z_score_limits[2]))
    mat[mat < z_score_limits[1]] <- z_score_limits[1]
    mat[mat > z_score_limits[2]] <- z_score_limits[2]
  }
  
  # ------------------------------
  # Order columns by group
  # ------------------------------
  ord_cols <- order(smeta[[group_col]], smeta[[sample_col]])
  mat   <- mat[, ord_cols, drop = FALSE]
  smeta <- smeta[ord_cols, , drop = FALSE]
  
  # ------------------------------
  # Annotation and colors
  # ------------------------------
  if (!all(names(group_palette) %in% levels(smeta[[group_col]]))) {
    group_palette <- setNames(group_palette[seq_along(levels(smeta[[group_col]]))],
                              levels(smeta[[group_col]]))
  }
  anno_df <- setNames(data.frame(smeta[[group_col]]), group_col)
  ha <- ComplexHeatmap::HeatmapAnnotation(
    df  = anno_df,
    col = setNames(list(group_palette), group_col),
    annotation_legend_param = list(title = group_col)
  )
  
  # heatmap colors
  if (is.null(z_score_limits)) {
    col_fun <- circlize::colorRamp2(c(-2, 0, 2), heatmap_colors)
  } else {
    col_fun <- circlize::colorRamp2(c(z_score_limits[1], 0, z_score_limits[2]), heatmap_colors)
  }
  
  # ------------------------------
  # Draw heatmap
  # ------------------------------
  ht <- ComplexHeatmap::Heatmap(
    mat,
    name = legend_title,
    col = col_fun,
    top_annotation = ha,
    show_row_names = show_rownames,
    show_column_names = show_colnames,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_cols,
    row_names_gp = grid::gpar(fontsize = 8),
    column_names_gp = grid::gpar(fontsize = 8),
    heatmap_legend_param = list(direction = "horizontal")
  )
  
  .msg("[heatmap] Drawing (legend bottom)...")
  ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
  
  if (return_matrix) {
    return(list(heatmap = ht, matrix = mat, annotation = smeta[, c(sample_col, group_col), drop = FALSE]))
  } else {
    invisible(ht)
  }
}



pseudobulk_heatmap_global_complexHeatmap <- function(
    agg,
    sample_meta,
    sample_col   = "sample_id",
    group_col    = "treatment",
    group_levels = c("CTRL","DEX","CORT"),
    group_palette = c(CTRL="#1b9e77", DEX="#d95f02", CORT="#7570b3"),
    genes        = NULL,
    top_n_var    = 50,
    transform    = c("logCPM","log2","none"),
    log2_offset  = 1,
    scale_rows   = TRUE,
    z_score_limits = NULL,           
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    legend_title = "Expression (z)",
    show_colnames = TRUE,
    show_rownames = TRUE,
    verbose      = TRUE,
    return_matrix = TRUE,
    quantile_normalize = FALSE,
    heatmap_colors = c("blue","white","red")
) {
  .msg <- function(...) if (isTRUE(verbose)) message(...)
  
  # --- simple quantile normalization ---
  .qn_base <- function(m) {
    m <- as.matrix(m)
    n <- nrow(m); p <- ncol(m)
    ord <- apply(m, 2, order, na.last = TRUE)
    m_sorted <- m
    for (j in seq_len(p)) m_sorted[, j] <- m[ord[, j], j]
    avg <- rowMeans(m_sorted, na.rm = TRUE)
    m_qn <- matrix(NA_real_, n, p, dimnames = dimnames(m))
    for (j in seq_len(p)) {
      inv <- integer(n); inv[ord[, j]] <- seq_len(n)
      m_qn[, j] <- avg[inv]
    }
    m_qn
  }
  
  transform <- match.arg(transform)
  
  # --- merge matrices ---
  .sum_by_rownames <- function(A, B) {
    A <- as.matrix(A); B <- as.matrix(B)
    all_genes <- union(rownames(A), rownames(B))
    common_samples <- intersect(colnames(A), colnames(B))
    if (!length(common_samples)) stop("No overlapping sample columns between matrices.")
    A2 <- matrix(0, nrow = length(all_genes), ncol = length(common_samples),
                 dimnames = list(all_genes, common_samples))
    B2 <- A2
    A2[rownames(A), common_samples] <- A[, common_samples, drop = FALSE]
    B2[rownames(B), common_samples] <- B[, common_samples, drop = FALSE]
    A2 + B2
  }
  .as_total_counts <- function(agg) {
    if (is.list(agg)) {
      mats <- Filter(function(m) (is.matrix(m) || inherits(m, "Matrix")) &&
                       nrow(m) > 0 && ncol(m) > 0, agg)
      if (!length(mats)) stop("Input list is empty or contains empty matrices.")
      .msg(sprintf("[heatmap] Merging %d cluster matrices into one global matrix...", length(mats)))
      mats <- lapply(mats, as.matrix)
      Reduce(.sum_by_rownames, mats)
    } else if (is.matrix(agg) || inherits(agg, "Matrix")) {
      as.matrix(agg)
    } else {
      stop("'agg' must be a matrix, sparse Matrix, or a list of such.")
    }
  }
  
  # --- counts + meta ---
  counts <- .as_total_counts(agg)
  stopifnot(all(c(sample_col, group_col) %in% colnames(sample_meta)))
  smeta <- as.data.frame(sample_meta)
  if (any(duplicated(smeta[[sample_col]]))) {
    .msg("[heatmap] Deduplicating sample metadata by ", sample_col)
    smeta <- smeta[!duplicated(smeta[[sample_col]]), , drop = FALSE]
  }
  keep_cols <- intersect(colnames(counts), smeta[[sample_col]])
  if (!length(keep_cols)) stop("No overlap between counts colnames and sample_meta[[sample_col]].")
  counts <- counts[, keep_cols, drop = FALSE]
  smeta  <- smeta[match(colnames(counts), smeta[[sample_col]]), , drop = FALSE]
  if (!is.null(group_levels)) smeta[[group_col]] <- factor(smeta[[group_col]], levels = group_levels)
  
  .msg(sprintf("[heatmap] Matrix: %d genes × %d samples", nrow(counts), ncol(counts)))
  .msg(sprintf("[heatmap] Groups: %s", paste(levels(smeta[[group_col]]), collapse = ", ")))
  
  # --- quantile normalize ---
  if (isTRUE(quantile_normalize)) {
    .msg("[heatmap] Quantile normalization: custom base R implementation")
    counts <- .qn_base(counts)
  }
  
  # --- transform ---
  mat <- counts
  if (transform == "logCPM") {
    .msg("[heatmap] Transform: edgeR::cpm(log=TRUE)")
    mat <- edgeR::cpm(mat, log = TRUE, prior.count = 1)
  } else if (transform == "log2") {
    .msg("[heatmap] Transform: log2(counts + offset), offset = ", log2_offset)
    mat <- log2(mat + log2_offset)
  } else {
    .msg("[heatmap] Transform: none (raw scale)")
  }
  
  # --- gene selection ---
  if (!is.null(genes)) {
    genes <- intersect(genes, rownames(mat))
    if (!length(genes)) stop("Provided 'genes' not found in matrix rownames.")
    mat <- mat[genes, , drop = FALSE]
    .msg(sprintf("[heatmap] Using provided genes: %d", length(genes)))
  } else {
    v <- matrixStats::rowVars(mat)
    ord <- order(v, decreasing = TRUE)
    n <- min(top_n_var, nrow(mat))
    mat <- mat[ord[seq_len(n)], , drop = FALSE]
    .msg(sprintf("[heatmap] Using top %d most variable genes", n))
  }
  
  # --- scale rows ---
  if (isTRUE(scale_rows)) {
    .msg("[heatmap] Scaling rows (z-score)")
    row_means <- rowMeans(mat, na.rm = TRUE)
    row_sds   <- matrixStats::rowSds(mat, na.rm = TRUE)
    row_sds[row_sds == 0] <- 1
    mat <- sweep(mat, 1, row_means, "-")
    mat <- sweep(mat, 1, row_sds,   "/")
    legend_title <- "Expression (z)"
  }
  if (!is.null(z_score_limits)) {
    .msg(sprintf("[heatmap] Clamping z-scores to range [%.1f, %.1f]",
                 z_score_limits[1], z_score_limits[2]))
    mat[mat < z_score_limits[1]] <- z_score_limits[1]
    mat[mat > z_score_limits[2]] <- z_score_limits[2]
  }
  
  # --- order cols ---
  ord_cols <- order(smeta[[group_col]], smeta[[sample_col]])
  mat   <- mat[, ord_cols, drop = FALSE]
  smeta <- smeta[ord_cols, , drop = FALSE]
  
  # --- annotation ---
  if (!all(names(group_palette) %in% levels(smeta[[group_col]]))) {
    group_palette <- setNames(group_palette[seq_along(levels(smeta[[group_col]]))],
                              levels(smeta[[group_col]]))
  }
  anno_df <- setNames(data.frame(smeta[[group_col]]), group_col)
  ha <- ComplexHeatmap::HeatmapAnnotation(
    df  = anno_df,
    col = setNames(list(group_palette), group_col),
    annotation_legend_param = list(title = group_col)
  )
  
  # --- colors ---
  if (is.null(z_score_limits)) {
    col_fun <- circlize::colorRamp2(c(-2, 0, 2), heatmap_colors)
  } else {
    col_fun <- circlize::colorRamp2(c(z_score_limits[1], 0, z_score_limits[2]), heatmap_colors)
  }
  
  # --- heatmap ---
  ht <- ComplexHeatmap::Heatmap(
    mat,
    name = legend_title,
    col = col_fun,
    top_annotation = ha,
    show_row_names = show_rownames,
    show_column_names = show_colnames,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_cols,
    row_names_side = "right",
    column_names_side = "top",   # <--- column names on top
    row_names_gp = grid::gpar(fontsize = 8),
    column_names_gp = grid::gpar(fontsize = 8),
    heatmap_legend_param = list(direction = "horizontal")
  )
  
  .msg("[heatmap] Drawing (legend bottom, colnames top)...")
  # ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
  
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "bottom", 
    annotation_legend_side = "bottom", 
    merge_legend = TRUE
  )
  
  # if (return_matrix) {
  #   return(list(heatmap = ht, matrix = mat, annotation = smeta[, c(sample_col, group_col), drop = FALSE]))
  # } else {
  #   invisible(ht)
  # }
}
