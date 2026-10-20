fgsea_heatmap_log10p_from_results <- function(
    fgsea_results,
    path_col       = "predict_cell_type",
    pval_col       = "pval",
    nes_col        = "NES",
    padj_max       = NULL,
    pval_floor     = 1e-300,
    fill_na        = NA,
    axis_x         = "celltype",    # "celltype" lub "cluster"
    axis_y         = "cluster",     # "cluster"  lub "celltype"
    numbers        = TRUE,
    numbers_digits = 3,             # miejsca po przecinku dla -log10(p)
    nes_digits     = 2,             # miejsca po przecinku dla NES
    numbers_size_pval = 9,          # rozmiar tekstu -log10(p)
    numbers_size_nes  = 9,          # rozmiar tekstu NES
    numbers_thresh = 1.3,           # próg koloru napisów wg -log10(p)
    numbers_bold   = FALSE,         # bold dla obu tekstów
    outline_max    = TRUE,          # włącz/wyłącz zielone ramki
    outline_mode   = c("max_nes","min_pval"),
    outline_reduce = c("by_row","by_col"),
    outline_col    = "darkgreen",
    outline_lwd    = 2,
    column_names_side = "top",
    column_names_rot  = 45,
    outer_border_col  = "black",
    outer_border_lwd  = 2,
    title_prefix      = "",
    row_dend_side_space    = grid::unit(20, "mm"),
    column_dend_side_space = grid::unit(20, "mm"),
    celltype_rename   = NULL,
    cluster_rename    = NULL,
    color_palette  = list(low = "white", mid = "lightcoral", high = "darkred"),
    return_matrix  = TRUE,
    main           = NULL,
    cluster_rows    = TRUE,         # << sterowanie dendrogramami przy wywołaniu
    cluster_columns = TRUE
){
  outline_mode   <- match.arg(outline_mode)
  outline_reduce <- match.arg(outline_reduce)
  
  stopifnot(requireNamespace("dplyr", quietly = TRUE))
  stopifnot(requireNamespace("tidyr", quietly = TRUE))
  stopifnot(requireNamespace("purrr", quietly = TRUE))
  stopifnot(requireNamespace("ComplexHeatmap", quietly = TRUE))
  stopifnot(requireNamespace("circlize", quietly = TRUE))
  stopifnot(requireNamespace("grid", quietly = TRUE))
  
  if (is.null(fgsea_results$results)) stop("fgsea_results must contain $results")
  if (is.null(main)) main <- paste0(title_prefix, "-log10(p-value) heatmap")
  
  # --- sizes of gene sets per cell type ---
  set_sizes <- NULL
  if (!is.null(fgsea_results$gene_lists_used)) {
    set_sizes <- vapply(fgsea_results$gene_lists_used, length, integer(1))
  }
  
  # --- long (cluster, celltype, log10p, nes) ---
  val_long <- purrr::imap_dfr(fgsea_results$results, function(df, cl){
    if (!nrow(df)) return(NULL)
    if (!pval_col %in% names(df)) stop("No column ", pval_col, " in results.")
    if (!nes_col  %in% names(df)) stop("No column ", nes_col,  " in results.")
    if (!is.null(padj_max)) {
      if (!"padj" %in% names(df)) stop("padj_max given but no 'padj' column in results.")
      df <- dplyr::filter(df, !is.na(.data$padj) & .data$padj <= padj_max)
    }
    pv <- df[[pval_col]]
    pv[is.na(pv)] <- NA_real_
    pv[pv <= 0]  <- pval_floor
    lpv <- -log10(pv)
    nes <- df[[nes_col]]
    dplyr::transmute(df,
                     cluster  = cl,
                     celltype = .data[[path_col]],
                     log10p   = lpv,
                     nes      = nes
    )
  })
  if (nrow(val_long) == 0) stop("No rows after filtering/transform.")
  
  all_clusters  <- names(fgsea_results$results)
  all_celltypes <- unique(unlist(purrr::map(fgsea_results$results, ~ .x[[path_col]])))
  
  # --- macierze bazowe: rows = clusters, cols = celltypes ---
  wide <- tidyr::pivot_wider(val_long, names_from = celltype, values_from = c(log10p, nes))
  dfw  <- as.data.frame(wide, check.names = FALSE)
  rownames(dfw) <- dfw$cluster; dfw$cluster <- NULL
  
  mat_log10p_base <- as.matrix(dfw[, grepl("^log10p_", colnames(dfw)), drop = FALSE])
  mat_nes_base    <- as.matrix(dfw[, grepl("^nes_",    colnames(dfw)), drop = FALSE])
  
  mat_log10p_base <- mat_log10p_base[all_clusters, paste0("log10p_", all_celltypes), drop = FALSE]
  mat_nes_base    <- mat_nes_base[   all_clusters, paste0("nes_",    all_celltypes), drop = FALSE]
  colnames(mat_log10p_base) <- all_celltypes
  colnames(mat_nes_base)    <- all_celltypes
  rownames(mat_log10p_base) <- all_clusters
  rownames(mat_nes_base)    <- all_clusters
  
  # --- zamiana osi ---
  if (axis_x == "celltype" && axis_y == "cluster") {
    mat_log10p <- mat_log10p_base
    mat_nes    <- mat_nes_base
  } else if (axis_x == "cluster" && axis_y == "celltype") {
    mat_log10p <- t(mat_log10p_base)
    mat_nes    <- t(mat_nes_base)
  } else {
    stop("Supported only: (axis_x='celltype', axis_y='cluster') OR (axis_x='cluster', axis_y='celltype').")
  }
  
  # --- NA fill (opcjonalnie) ---
  if (!is.na(fill_na)) {
    mat_log10p[is.na(mat_log10p)] <- fill_na
    mat_nes[is.na(mat_nes)]       <- fill_na
  }
  
  # --- RENAME + n= (najpierw celltype, potem cluster) ---
  
  ## zbuduj mapę ładnych nazw celltype
  pretty_celltypes <- all_celltypes
  if (!is.null(celltype_rename)) {
    if (is.list(celltype_rename)) celltype_rename <- unlist(celltype_rename)
    stopifnot(is.character(celltype_rename))
    pretty_celltypes <- dplyr::recode(all_celltypes, !!!celltype_rename, .default = all_celltypes)
  }
  ct_map <- setNames(pretty_celltypes, all_celltypes)
  
  if (axis_x == "celltype") {
    # celltype na KOLUMNACH
    ct_ids <- colnames(mat_log10p)
    pretty <- unname(ct_map[ct_ids]); pretty[is.na(pretty)] <- ct_ids
    if (!is.null(set_sizes)) {
      n_vals <- set_sizes[ct_ids]
      colnames(mat_log10p) <- paste0(pretty, "\n(n=", n_vals, ")")
      colnames(mat_nes)    <- colnames(mat_log10p)
    } else {
      colnames(mat_log10p) <- pretty
      colnames(mat_nes)    <- pretty
    }
  } else if (axis_x == "cluster" && axis_y == "celltype") {
    # celltype na WIERSZACH
    ct_ids <- rownames(mat_log10p)
    pretty <- unname(ct_map[ct_ids]); pretty[is.na(pretty)] <- ct_ids
    if (!is.null(set_sizes)) {
      n_vals <- set_sizes[ct_ids]
      rownames(mat_log10p) <- paste0(pretty, "\n(n=", n_vals, ")")
      rownames(mat_nes)    <- rownames(mat_log10p)
    } else {
      rownames(mat_log10p) <- pretty
      rownames(mat_nes)    <- pretty
    }
  }
  
  ## zbuduj mapę ładnych nazw clusterów
  pretty_clusters <- all_clusters
  if (!is.null(cluster_rename)) {
    if (is.list(cluster_rename)) cluster_rename <- unlist(cluster_rename)
    stopifnot(is.character(cluster_rename))
    pretty_clusters <- dplyr::recode(all_clusters, !!!cluster_rename, .default = all_clusters)
  }
  cl_map <- setNames(pretty_clusters, all_clusters)
  
  # zastosuj rename klastrów na wymiarze, gdzie aktualnie leżą
  if (axis_y == "cluster") {
    # clustery na WIERSZACH
    cl_ids <- rownames(mat_log10p)
    cl_ids_clean <- sub("\\n\\(n=.*\\)$", "", cl_ids)  # asekuracja
    pretty <- unname(cl_map[cl_ids_clean]); pretty[is.na(pretty)] <- cl_ids_clean
    rownames(mat_log10p) <- pretty
    rownames(mat_nes)    <- pretty
  } else if (axis_x == "cluster") {
    # clustery na KOLUMNACH
    cl_ids <- colnames(mat_log10p)
    cl_ids_clean <- sub("\\n\\(n=.*\\)$", "", cl_ids)
    pretty <- unname(cl_map[cl_ids_clean]); pretty[is.na(pretty)] <- cl_ids_clean
    colnames(mat_log10p) <- pretty
    colnames(mat_nes)    <- pretty
  }
  
  # --- kolory (-log10 p) ---
  rng  <- range(mat_log10p, na.rm = TRUE)
  vmin <- 0
  vmax <- max(1e-9, rng[2])
  midv <- (vmin + vmax) / 2
  col_fun <- circlize::colorRamp2(
    c(vmin, midv, vmax),
    c(color_palette$low, color_palette$mid, color_palette$high)
  )
  
  # --- numery do wyświetlania ---
  num_mat_log10p_raw <- mat_log10p
  num_mat_log10p_txt <- matrix(
    sprintf(paste0("%.", numbers_digits, "f"), as.vector(mat_log10p)),
    nrow = nrow(mat_log10p), ncol = ncol(mat_log10p),
    dimnames = dimnames(mat_log10p)
  )
  num_mat_nes_txt <- matrix(
    sprintf(paste0("%.", nes_digits, "f"), as.vector(mat_nes)),
    nrow = nrow(mat_nes), ncol = ncol(mat_nes),
    dimnames = dimnames(mat_nes)
  )
  
  # --- które komórki podkreślić ---
  crit_mat <- if (outline_mode == "max_nes") mat_nes else mat_log10p
  pos_row <- pos_col <- NULL
  if (outline_max) {
    if (outline_reduce == "by_row") {
      pos_row <- apply(crit_mat, 1, function(x) if (all(is.na(x))) NA_integer_ else which.max(x))
    } else {
      pos_col <- apply(crit_mat, 2, function(x) if (all(is.na(x))) NA_integer_ else which.max(x))
    }
  }
  
  # --- legenda ramki ---
  lgd_outline <- ComplexHeatmap::Legend(
    title  = NULL,
    labels = if (outline_mode == "max_nes") "max NES (normalized enrichment score)" else "min p-value",
    type   = "grid",
    graphics = list(function(x, y, w, h) {
      grid::grid.rect(x, y, w, h, gp = grid::gpar(fill = NA, col = outline_col, lwd = outline_lwd))
    }),
    grid_width  = grid::unit(10, "mm"),
    grid_height = grid::unit(10, "mm"),
    direction   = "horizontal",
    labels_gp   = grid::gpar(fontsize = 12)
  )
  
  # --- heatmapa ---
  ht <- ComplexHeatmap::Heatmap(
    mat_log10p,
    name = "-log10(p-value)", col = col_fun,
    cluster_rows    = cluster_rows,
    cluster_columns = cluster_columns,
    show_row_names = TRUE, show_column_names = TRUE,
    column_title = main,
    column_names_side = column_names_side,
    column_names_rot  = column_names_rot,
    heatmap_legend_param = list(
      direction       = "horizontal",
      title_position  = "leftcenter",   # tytuł po lewej stronie paska
      legend_width    = grid::unit(6, "cm"),
      title_gp        = grid::gpar(fontsize = 12, fontface = "bold"),
      labels_gp       = grid::gpar(fontsize = 10)
    ),
    # dendrogramy
    row_dend_side      = "left",
    row_dend_width     = row_dend_side_space,
    column_dend_side   = "top",
    column_dend_height = column_dend_side_space,
    # rysowanie komórek
    cell_fun = function(j, i, x, y, w, h, col) {
      if (numbers && is.finite(num_mat_log10p_raw[i, j])) {
        col_txt <- if (num_mat_log10p_raw[i, j] > numbers_thresh) "white" else "black"
        # góra: -log10(pval)
        grid::grid.text(
          num_mat_log10p_txt[i, j], x, y + h*0.2,
          gp = grid::gpar(
            col = col_txt,
            fontsize = numbers_size_pval,
            fontface = if (numbers_bold) "bold" else "plain"
          )
        )
        # dół: NES
        grid::grid.text(
          num_mat_nes_txt[i, j], x, y - h*0.2,
          gp = grid::gpar(
            col = col_txt,
            fontsize = numbers_size_nes,
            fontface = if (numbers_bold) "bold" else "plain"
          )
        )
      }
      if (outline_max) {
        draw_row <- !is.null(pos_row) && !is.na(pos_row[i]) && j == pos_row[i]
        draw_col <- !is.null(pos_col) && !is.na(pos_col[j]) && i == pos_col[j]
        if (draw_row || draw_col) {
          grid::grid.rect(x = x, y = y, width = w, height = h,
                          gp = grid::gpar(col = outline_col, fill = NA, lwd = outline_lwd))
        }
      }
    }
  )
  
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side     = "bottom",
    annotation_legend_side  = "bottom",
    show_heatmap_legend     = TRUE,
    show_annotation_legend  = TRUE,
    merge_legends           = TRUE,
    annotation_legend_list  = list(lgd_outline)
  )
  
  ComplexHeatmap::decorate_heatmap_body("-log10(p-value)", {
    grid::grid.rect(gp = grid::gpar(col = outer_border_col, fill = NA, lwd = outer_border_lwd))
  })
  
  if (return_matrix) {
    invisible(list(
      matrix_log10p = mat_log10p,
      matrix_nes    = mat_nes,
      heatmap       = ht,
      axis_x        = axis_x,
      axis_y        = axis_y
    ))
  } else {
    invisible(ht)
  }
}
