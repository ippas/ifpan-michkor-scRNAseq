fgsea_heatmap_from_results <- function(
    fgsea_results,
    path_col       = "predict_cell_type",
    value_col      = "NES",
    padj_max       = NULL,
    fill_na        = NA,
    zscore_row     = FALSE,
    zscore_limit   = NULL,
    axis_x         = "celltype",
    axis_y         = "cluster",
    numbers        = TRUE,
    numbers_digits = 2,
    numbers_thresh = 1.5,
    numbers_bold   = FALSE,
    outline_max    = TRUE,
    outline_col    = "darkgreen",
    outline_lwd    = 2,
    column_names_side = "top",
    column_names_rot  = 45,
    outer_border_col  = "black",
    outer_border_lwd  = 2,
    title_prefix      = "",
    row_dend_side_space = grid::unit(5, "mm"),
    main           = NULL,
    celltype_rename = NULL,
    color_palette  = list(low = "darkblue", mid = "white", high = "darkred"),
    return_matrix  = TRUE
){
  stopifnot(requireNamespace("dplyr", quietly = TRUE))
  stopifnot(requireNamespace("tidyr", quietly = TRUE))
  stopifnot(requireNamespace("purrr", quietly = TRUE))
  stopifnot(requireNamespace("ComplexHeatmap", quietly = TRUE))
  stopifnot(requireNamespace("circlize", quietly = TRUE))
  stopifnot(requireNamespace("grid", quietly = TRUE))
  
  if (is.null(fgsea_results$results)) stop("fgsea_results must contain $results")
  if (is.null(main)) main <- paste0(title_prefix, "Normalized Enrichment Score heatmap")
  
  # --- sizes of gene sets per cell type ---
  set_sizes <- NULL
  if (!is.null(fgsea_results$gene_lists_used)) {
    set_sizes <- vapply(fgsea_results$gene_lists_used, length, integer(1))
  }
  
  # --- long -> wide ---
  all_clusters  <- names(fgsea_results$results)
  all_celltypes <- unique(unlist(purrr::map(fgsea_results$results, ~ .x[[path_col]])))
  
  val_long <- purrr::imap_dfr(fgsea_results$results, function(df, cl){
    if (!nrow(df)) return(NULL)
    if (!value_col %in% names(df)) stop("No column ", value_col, " in results.")
    if (!is.null(padj_max)) {
      df <- dplyr::filter(df, !is.na(.data$padj) & .data$padj <= padj_max)
    }
    dplyr::transmute(df, cluster = cl, celltype = .data[[path_col]], value = .data[[value_col]])
  })
  
  if (axis_x == "celltype" && axis_y == "cluster") {
    wide <- tidyr::pivot_wider(val_long, names_from = celltype, values_from = value)
    mat  <- as.data.frame(wide, check.names = FALSE)
    rownames(mat) <- mat$cluster; mat$cluster <- NULL
    mat  <- as.matrix(mat)
    mat  <- mat[all_clusters, all_celltypes, drop = FALSE]
  } else {
    stop("Na razie wspieram tylko celltype (x) × cluster (y)")
  }
  
  if (!is.na(fill_na)) mat[is.na(mat)] <- fill_na
  
  # --- nazwy kolumn z n= ---
  orig_names <- all_celltypes
  pretty_names <- orig_names
  if (!is.null(celltype_rename)) {
    if (is.list(celltype_rename)) celltype_rename <- unlist(celltype_rename)
    stopifnot(is.character(celltype_rename))
    pretty_names <- dplyr::recode(orig_names, !!!celltype_rename, .default = orig_names)
  }
  if (!is.null(set_sizes)) {
    n_vals <- set_sizes[orig_names]
    colnames(mat) <- paste0(pretty_names, "\n(n=", n_vals, ")")
  } else {
    colnames(mat) <- pretty_names
  }
  
  # --- kolorowanie ---
  mm <- mat
  if (zscore_row) {
    row_center <- rowMeans(mm, na.rm = TRUE)
    row_sd     <- apply(mm, 1, sd, na.rm = TRUE)
    mm <- sweep(mm, 1, row_center, "-")
    mm <- sweep(mm, 1, ifelse(row_sd == 0 | is.na(row_sd), 1, row_sd), "/")
    if (!is.null(zscore_limit)) {
      mm <- t(apply(mm, 1, function(x){
        rng <- range(x, na.rm = TRUE)
        if (diff(rng) == 0) return(rep(0, length(x)))
        scales::rescale(x, to = zscore_limit, from = rng)
      }))
      rownames(mm) <- rownames(mat); colnames(mm) <- colnames(mat)
    }
  }
  
  rng  <- range(mm, na.rm = TRUE)
  vmax <- max(abs(rng))
  col_fun <- circlize::colorRamp2(
    c(-vmax, 0, vmax),
    c(color_palette$low, color_palette$mid, color_palette$high)
  )
  
  num_mat <- round(mat, numbers_digits)
  max_pos <- apply(mat, 1, function(x) if (all(is.na(x))) NA_integer_ else which.max(x))
  
  ht <- ComplexHeatmap::Heatmap(
    mm, name = "Normalized Enrichment Score", col = col_fun,
    cluster_rows = TRUE, cluster_columns = TRUE,
    show_row_names = TRUE, show_column_names = TRUE,
    column_title = main,
    column_names_side = column_names_side,
    column_names_rot  = column_names_rot,
    heatmap_legend_param = list(
      direction = "horizontal",
      legend_width = grid::unit(6, "cm"),
      title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 10)
    ),
    row_dend_side = "left",
    row_dend_width = row_dend_side_space,
    cell_fun = function(j, i, x, y, w, h, col) {
      if (numbers && is.finite(num_mat[i, j])) {
        grid::grid.text(
          as.character(num_mat[i, j]), x, y,
          gp = grid::gpar(
            col = if (abs(num_mat[i, j]) > numbers_thresh) "white" else "black",
            fontsize = 10,
            fontface = if (numbers_bold) "bold" else "plain"
          )
        )
      }
      if (outline_max && !is.na(max_pos[i]) && j == max_pos[i]) {
        grid::grid.rect(x = x, y = y, width = w, height = h,
                        gp = grid::gpar(col = outline_col, fill = NA, lwd = outline_lwd))
      }
    }
  )
  
  # --- legendy ---
  lgd_color <- ComplexHeatmap::Legend(
    col_fun = col_fun,
    title   = "Normalized Enrichment Score",
    direction = "horizontal",
    legend_width = grid::unit(6, "cm"),
    title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
    labels_gp = grid::gpar(fontsize = 10)
  )
  
  lgd_outline <- ComplexHeatmap::Legend(
    title = NULL,
    labels = "max NES per cluster",
    type = "grid",
    graphics = list(function(x, y, w, h) {
      grid::grid.rect(x, y, w, h, gp = grid::gpar(fill = NA, col = outline_col, lwd = outline_lwd))
    }),
    grid_width  = grid::unit(5, "mm"),
    grid_height = grid::unit(5, "mm"),
    direction   = "horizontal"
  )
  lgd_outline <- ComplexHeatmap::Legend(
    title = NULL,
    labels = "max NES per cluster",
    type = "grid",
    graphics = list(function(x, y, w, h) {
      grid::grid.rect(x, y, w, h, gp = grid::gpar(fill = NA, col = outline_col, lwd = outline_lwd))
    }),
    grid_width  = grid::unit(5, "mm"),
    grid_height = grid::unit(5, "mm"),
    direction   = "horizontal"
  )
  
  # Rysujemy heatmapę z WŁĄCZONĄ legendą kolorów,
  # a legendę z ramką przekazujemy jako annotation legend.
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side     = "bottom",
    annotation_legend_side  = "bottom",
    show_heatmap_legend     = TRUE,              # <- ważne
    show_annotation_legend  = TRUE,
    merge_legends           = TRUE,
    annotation_legend_list  = list(lgd_outline)  # <- dokładamy ramkę
  )
  
  # czarna ramka dookoła heatmapy
  ComplexHeatmap::decorate_heatmap_body("Normalized Enrichment Score", {
    grid::grid.rect(gp = grid::gpar(col = outer_border_col, fill = NA, lwd = outer_border_lwd))
  })
  
  if (return_matrix) {
    invisible(list(matrix_raw = mat, matrix_color = mm, heatmap = ht))
  } else {
    invisible(ht)
  }
}
