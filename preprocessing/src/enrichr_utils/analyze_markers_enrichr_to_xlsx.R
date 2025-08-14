analyze_markers_enrichr_to_xlsx <- function(
    seurat_obj,
    xlsx_path,
    cluster_col           = "seurat_clusters",
    assay                 = NULL,
    min.pct               = 0.5,
    test.use              = "wilcox",
    p_adj_thresh          = 0.05,
    logfc_up_thresh       = 0.5,
    logfc_down_thresh     = -0.5,
    gene_id_sep           = "-",
    top_n_featureplots    = 9,
    featureplot_pt_size   = 0.5,
    enrichr_databases     = c("CellMarker_2024", "CellMarker_Augmented_2021",
                              "Human_Gene_Atlas", "ARCHS4_Cell-lines"),
    enr_min_genes         = 2,
    enr_n_top             = 10,
    enr_fill_limits       = c(0, 10),
    enr_squish            = TRUE,
    enr_symmetric_limits  = FALSE,
    axis_text_size        = 12,
    axis_title_size       = 12,
    legend_text_size      = 10,
    legend_title_size     = 12,
    plot_title_size       = 14,
    label_text_size       = 4,
    seed                  = 0,
    verbose               = TRUE
) {
  stopifnot("Seurat" %in% .packages(all.available = TRUE))
  if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
    stop(sprintf("Column '%s' not found in seurat_obj@meta.data.", cluster_col))
  }
  if (is.null(xlsx_path) || !nzchar(xlsx_path)) {
    stop("Provide a valid 'xlsx_path'.")
  }
  
  set.seed(seed)
  
  original_assay <- Seurat::DefaultAssay(seurat_obj)
  original_idents <- Seurat::Idents(seurat_obj)
  
  on.exit({
    Seurat::DefaultAssay(seurat_obj) <<- original_assay
    if (!is.null(original_idents)) Seurat::Idents(seurat_obj) <<- original_idents
  }, add = TRUE)
  
  if (!is.null(assay)) {
    if (!assay %in% names(seurat_obj@assays)) {
      stop(sprintf("Assay '%s' not found in seurat_obj.", assay))
    }
    Seurat::DefaultAssay(seurat_obj) <- assay
  }
  
  Seurat::Idents(seurat_obj) <- seurat_obj@meta.data[[cluster_col]]
  
  # 1/6: FindAllMarkers
  if (verbose) message("[1/6] Running FindAllMarkers for '", cluster_col, "' ...")
  markers_df <- Seurat::FindAllMarkers(
    object  = seurat_obj,
    only.pos = FALSE,
    min.pct = min.pct,
    test.use = test.use
  )
  
  # 2/6: Prepare marker lists
  if (verbose) message("[2/6] Preparing marker lists ...")
  markers_df <- markers_df %>%
    tibble::rownames_to_column(var = "rowname") %>%
    dplyr::select(-rowname) %>%
    dplyr::filter(.data$p_val_adj < p_adj_thresh) %>%
    dplyr::mutate(name_id = .data$gene) %>%
    dplyr::select(name_id, dplyr::everything()) %>%
    tidyr::separate(
      col     = "gene",
      into    = c("gene_id", "gene_symbol"),
      sep     = gene_id_sep,
      remove  = TRUE,
      extra   = "merge",
      fill    = "right"
    ) %>%
    dplyr::mutate(cluster = paste0("cluster_", .data$cluster))
  
  markers_up_list <- markers_df %>%
    dplyr::filter(.data$avg_log2FC >  logfc_up_thresh) %>%
    dplyr::group_by(.data$cluster) %>%
    dplyr::group_split() %>%
    rlang::set_names(unique(markers_df$cluster)[unique(markers_df$cluster) %in% dplyr::group_keys(dplyr::group_by(dplyr::filter(markers_df, avg_log2FC > logfc_up_thresh), cluster))$cluster])
  
  markers_down_list <- markers_df %>%
    dplyr::filter(.data$avg_log2FC <  logfc_down_thresh) %>%
    dplyr::group_by(.data$cluster) %>%
    dplyr::group_split() %>%
    rlang::set_names(unique(markers_df$cluster)[unique(markers_df$cluster) %in% dplyr::group_keys(dplyr::group_by(dplyr::filter(markers_df, avg_log2FC < logfc_down_thresh), cluster))$cluster])
  
  if (length(markers_up_list)) {
    names(markers_up_list) <- purrr::map_chr(markers_up_list, ~ unique(.x$cluster))
  }
  if (length(markers_down_list)) {
    names(markers_down_list) <- purrr::map_chr(markers_down_list, ~ unique(.x$cluster))
  }
  
  # 3/6: FeaturePlots
  if (verbose) message("[3/6] Preparing FeaturePlots for top ", top_n_featureplots, " genes per cluster ...")
  top_markersUp_plots_list <- markers_up_list %>%
    purrr::imap(function(df_cl, cl_name) {
      genes <- df_cl %>%
        dplyr::slice_head(n = top_n_featureplots) %>%
        dplyr::pull(name_id) %>%
        unique()
      if (length(genes) == 0) return(NULL)
      (\(genes_vec) Seurat::FeaturePlot(
        seurat_obj,
        features = genes_vec,
        pt.size  = featureplot_pt_size
      ))(genes)
    })
  
  # 4/6: Enrichr
  if (verbose) message("[4/6] Running Enrichr for databases: ",
                       paste(enrichr_databases, collapse = ", "), " ...")
  enrichrResults_perDatabase_perCluster <- lapply(
    enrichr_databases,
    function(db) {
      run_enrichr_for_all_clusters(
        markers_list = markers_up_list,
        database     = db,
        gene_col     = "gene_symbol"
      )
    }
  )
  names(enrichrResults_perDatabase_perCluster) <- enrichr_databases
  
  # 5/6: Enrichr plots
  if (verbose) message("[5/6] Building Enrichr plots ...")
  clusters <- sort(unique(c(names(markers_up_list), names(markers_down_list))))
  plots_perCluster_perDatabase <- lapply(clusters, function(cl) {
    available_dbs <- enrichr_databases[
      purrr::map_lgl(enrichr_databases, ~ !is.null(enrichrResults_perDatabase_perCluster[[.x]][[cl]]))
    ]
    if (!length(available_dbs)) return(NULL)
    
    patchwork::wrap_plots(
      plotlist = lapply(available_dbs, function(db) {
        barplot_enrichr_results(
          df                 = enrichrResults_perDatabase_perCluster[[db]][[cl]],
          min_genes          = enr_min_genes,
          n_top              = enr_n_top,
          title              = paste("Cluster", cl, "–", db),
          axis_text_size     = axis_text_size,
          axis_title_size    = axis_title_size,
          legend_text_size   = legend_text_size,
          legend_title_size  = legend_title_size,
          plot_title_size    = plot_title_size,
          label_text_size    = label_text_size,
          fill_limits        = enr_fill_limits,
          squish_to_limits   = enr_squish,
          symmetric_limits   = enr_symmetric_limits
        )
      }),
      ncol = 1
    ) +
      patchwork::plot_layout(guides = "collect") &
      ggplot2::theme(legend.position = "bottom")
  })
  names(plots_perCluster_perDatabase) <- clusters
  
  # 6/6: Save XLSX
  if (verbose) message("[6/6] Writing XLSX to: ", xlsx_path)
  write_cluster_listAndTwoPlots_to_xlsx(
    cluster_list = markers_up_list,
    plot_list_a  = top_markersUp_plots_list,
    plot_list_b  = plots_perCluster_perDatabase,
    file_path    = xlsx_path,
    img_width_a  = 11, img_height_a = 9,
    img_width_b  = 12, img_height_b = 13,
    gap_rows     = 45
  )
  
  if (verbose) message("✅ Done.")
  
  invisible(list(
    markers_up_list   = markers_up_list,
    markers_down_list = markers_down_list,
    enrichr_results   = enrichrResults_perDatabase_perCluster,
    feature_plots     = top_markersUp_plots_list,
    enrichr_plots     = plots_perCluster_perDatabase,
    xlsx_path         = xlsx_path
  ))
}

analyze_markers_enrichr_to_xlsx <- function(
    seurat_obj,
    xlsx_path,
    cluster_col           = "seurat_clusters",
    assay                 = NULL,
    only.pos              = FALSE,               # NEW: przekazywane do FindAllMarkers
    min.pct               = 0.5,
    test.use              = "wilcox",
    p_adj_thresh          = 0.05,
    logfc_up_thresh       = 0.5,
    logfc_down_thresh     = -0.5,
    gene_id_sep           = "-",
    top_n_featureplots    = 9,
    featureplot_pt_size   = 0.5,
    enrichr_databases     = c("CellMarker_2024", "CellMarker_Augmented_2021",
                              "Human_Gene_Atlas", "ARCHS4_Cell-lines"),
    enr_min_genes         = 2,
    enr_n_top             = 10,
    enr_fill_limits       = c(0, 10),
    enr_squish            = TRUE,
    enr_symmetric_limits  = FALSE,
    axis_text_size        = 12,
    axis_title_size       = 12,
    legend_text_size      = 10,
    legend_title_size     = 12,
    plot_title_size       = 14,
    label_text_size       = 4,
    seed                  = 0,
    verbose               = TRUE
) {
  stopifnot("Seurat" %in% .packages(all.available = TRUE))
  if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
    stop(sprintf("Column '%s' not found in seurat_obj@meta.data.", cluster_col))
  }
  if (is.null(xlsx_path) || !nzchar(xlsx_path)) {
    stop("Provide a valid 'xlsx_path'.")
  }
  
  set.seed(seed)
  
  # Zapisz stan początkowy
  original_assay  <- Seurat::DefaultAssay(seurat_obj)
  original_idents <- Seurat::Idents(seurat_obj)
  
  # Przywracanie w try() żeby nie wywalało błędu
  on.exit({
    try({
      Seurat::DefaultAssay(seurat_obj) <- original_assay
      if (!is.null(original_idents)) {
        Seurat::Idents(seurat_obj) <- original_idents
      }
    }, silent = TRUE)
  }, add = TRUE)
  
  # Ustaw assay jeśli podano
  if (!is.null(assay)) {
    if (!assay %in% names(seurat_obj@assays)) {
      stop(sprintf("Assay '%s' not found in seurat_obj.", assay))
    }
    Seurat::DefaultAssay(seurat_obj) <- assay
  }
  
  # Ustaw identyfikatory na kolumnę klastrów
  Seurat::Idents(seurat_obj) <- seurat_obj@meta.data[[cluster_col]]
  
  # [1/7] FindAllMarkers
  if (verbose) message("[1/7] Running FindAllMarkers for '", cluster_col, "' ...")
  markers_df <- Seurat::FindAllMarkers(
    object  = seurat_obj,
    only.pos = only.pos,     # NEW
    min.pct = min.pct,
    test.use = test.use
  )
  
  # [2/7] Prepare marker lists
  if (verbose) message("[2/7] Preparing marker lists ...")
  markers_df <- markers_df %>%
    tibble::rownames_to_column(var = "rowname") %>%
    dplyr::select(-rowname) %>%
    dplyr::filter(.data$p_val_adj < p_adj_thresh) %>%
    dplyr::mutate(name_id = .data$gene) %>%
    dplyr::select(name_id, dplyr::everything()) %>%
    tidyr::separate(
      col     = "gene",
      into    = c("gene_id", "gene_symbol"),
      sep     = gene_id_sep,
      remove  = TRUE,
      extra   = "merge",
      fill    = "right"
    ) %>%
    dplyr::mutate(cluster = paste0("cluster_", .data$cluster))
  
  markers_up_list <- markers_df %>%
    dplyr::filter(.data$avg_log2FC >  logfc_up_thresh) %>%
    dplyr::group_by(.data$cluster) %>%
    dplyr::group_split()
  if (length(markers_up_list)) {
    names(markers_up_list) <- purrr::map_chr(markers_up_list, ~ unique(.x$cluster))
  }
  
  markers_down_list <- markers_df %>%
    dplyr::filter(.data$avg_log2FC <  logfc_down_thresh) %>%
    dplyr::group_by(.data$cluster) %>%
    dplyr::group_split()
  if (length(markers_down_list)) {
    names(markers_down_list) <- purrr::map_chr(markers_down_list, ~ unique(.x$cluster))
  }
  
  # [3/7] FeaturePlots
  if (verbose) message("[3/7] Preparing FeaturePlots for top ", top_n_featureplots, " genes per cluster ...")
  top_markersUp_plots_list <- markers_up_list %>%
    purrr::imap(function(df_cl, cl_name) {
      genes <- df_cl %>%
        dplyr::slice_head(n = top_n_featureplots) %>%
        dplyr::pull(name_id) %>%
        unique()
      if (length(genes) == 0) return(NULL)
      (\(genes_vec) Seurat::FeaturePlot(
        seurat_obj,
        features = genes_vec,
        pt.size  = featureplot_pt_size
      ))(genes)
    })
  
  # [4/7] Enrichr
  if (verbose) message("[4/7] Running Enrichr for databases: ",
                       paste(enrichr_databases, collapse = ", "), " ...")
  enrichrResults_perDatabase_perCluster <- lapply(
    enrichr_databases,
    function(db) {
      run_enrichr_for_all_clusters(
        markers_list = markers_up_list,
        database     = db,
        gene_col     = "gene_symbol"
      )
    }
  )
  names(enrichrResults_perDatabase_perCluster) <- enrichr_databases
  
  # [5/7] Enrichr plots
  if (verbose) message("[5/7] Building Enrichr plots ...")
  clusters <- sort(unique(c(names(markers_up_list), names(markers_down_list))))
  plots_perCluster_perDatabase <- lapply(clusters, function(cl) {
    available_dbs <- enrichr_databases[
      purrr::map_lgl(enrichr_databases, ~ !is.null(enrichrResults_perDatabase_perCluster[[.x]][[cl]]))
    ]
    if (!length(available_dbs)) return(NULL)
    
    patchwork::wrap_plots(
      plotlist = lapply(available_dbs, function(db) {
        barplot_enrichr_results(
          df                 = enrichrResults_perDatabase_perCluster[[db]][[cl]],
          min_genes          = enr_min_genes,
          n_top              = enr_n_top,
          title              = paste("Cluster", cl, "–", db),
          axis_text_size     = axis_text_size,
          axis_title_size    = axis_title_size,
          legend_text_size   = legend_text_size,
          legend_title_size  = legend_title_size,
          plot_title_size    = plot_title_size,
          label_text_size    = label_text_size,
          fill_limits        = enr_fill_limits,
          squish_to_limits   = enr_squish,
          symmetric_limits   = enr_symmetric_limits
        )
      }),
      ncol = 1
    ) +
      patchwork::plot_layout(guides = "collect") &
      ggplot2::theme(legend.position = "bottom")
  })
  names(plots_perCluster_perDatabase) <- clusters
  
  # [6/7] Save XLSX
  if (verbose) message("[6/7] Writing XLSX to: ", xlsx_path)
  write_cluster_listAndTwoPlots_to_xlsx(
    cluster_list = markers_up_list,
    plot_list_a  = top_markersUp_plots_list,
    plot_list_b  = plots_perCluster_perDatabase,
    file_path    = xlsx_path,
    img_width_a  = 11, img_height_a = 9,
    img_width_b  = 12, img_height_b = 13,
    gap_rows     = 45
  )
  
  # [7/7] Build per-cluster summaries
  if (verbose) message("[7/7] Building per-cluster summaries ...")
  cluster_summaries <- purrr::map(clusters, function(cl) {
    list(
      markers_up_df    = markers_up_list[[cl]]   %||% NULL,
      markers_down_df  = markers_down_list[[cl]] %||% NULL,
      feature_plot_top = top_markersUp_plots_list[[cl]] %||% NULL,
      enrichr_plot     = plots_perCluster_perDatabase[[cl]] %||% NULL
    )
  })
  names(cluster_summaries) <- clusters
  
  if (verbose) message("✅ Done.")
  
  invisible(list(
    markers_up_list    = markers_up_list,
    markers_down_list  = markers_down_list,
    enrichr_results    = enrichrResults_perDatabase_perCluster,
    feature_plots      = top_markersUp_plots_list,
    enrichr_plots      = plots_perCluster_perDatabase,
    xlsx_path          = xlsx_path,
    cluster_summaries  = cluster_summaries
  ))
}


analyze_markers_enrichr_to_xlsx <- function(
    seurat_obj,
    xlsx_path,
    cluster_col           = "seurat_clusters",
    assay                 = NULL,
    only.pos              = FALSE,               # NEW: przekazywane do FindAllMarkers
    min.pct               = 0.5,
    test.use              = "wilcox",
    p_adj_thresh          = 0.05,
    logfc_up_thresh       = 0.5,
    logfc_down_thresh     = -0.5,
    gene_id_sep           = "-",
    top_n_featureplots    = 9,
    featureplot_pt_size   = 0.5,
    enrichr_databases     = c("CellMarker_2024", "CellMarker_Augmented_2021",
                              "Human_Gene_Atlas", "ARCHS4_Cell-lines"),
    enr_min_genes         = 2,
    enr_n_top             = 10,
    enr_fill_limits       = c(0, 10),
    enr_squish            = TRUE,
    enr_symmetric_limits  = FALSE,
    axis_text_size        = 12,
    axis_title_size       = 12,
    legend_text_size      = 10,
    legend_title_size     = 12,
    plot_title_size       = 14,
    label_text_size       = 4,
    seed                  = 0,
    verbose               = TRUE
) {
  stopifnot("Seurat" %in% .packages(all.available = TRUE))
  if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
    stop(sprintf("Column '%s' not found in seurat_obj@meta.data.", cluster_col))
  }
  if (is.null(xlsx_path) || !nzchar(xlsx_path)) {
    stop("Provide a valid 'xlsx_path'.")
  }
  
  set.seed(seed)
  
  # --- zapisz stan początkowy
  original_assay  <- Seurat::DefaultAssay(seurat_obj)
  original_idents <- Seurat::Idents(seurat_obj)
  
  on.exit({
    try({
      Seurat::DefaultAssay(seurat_obj) <- original_assay
      if (!is.null(original_idents)) {
        Seurat::Idents(seurat_obj) <- original_idents
      }
    }, silent = TRUE)
  }, add = TRUE)
  
  # --- ustaw assay jeśli podano
  if (!is.null(assay)) {
    if (!assay %in% names(seurat_obj@assays)) {
      stop(sprintf("Assay '%s' not found in seurat_obj.", assay))
    }
    Seurat::DefaultAssay(seurat_obj) <- assay
  }
  
  # --- ustaw Idents
  Seurat::Idents(seurat_obj) <- seurat_obj@meta.data[[cluster_col]]
  
  # [1/7] FindAllMarkers
  if (verbose) message("[1/7] Running FindAllMarkers for '", cluster_col, "' ...")
  markers_df <- Seurat::FindAllMarkers(
    object   = seurat_obj,
    only.pos = only.pos,
    min.pct  = min.pct,
    test.use = test.use
  )
  
  # [2/7] Prepare marker lists
  if (verbose) message("[2/7] Preparing marker lists ...")
  markers_df <- markers_df %>%
    tibble::rownames_to_column(var = "rowname") %>%
    dplyr::select(-rowname) %>%
    dplyr::filter(.data$p_val_adj < p_adj_thresh) %>%
    dplyr::mutate(name_id = .data$gene) %>%
    dplyr::select(name_id, dplyr::everything()) %>%
    tidyr::separate(
      col     = "gene",
      into    = c("gene_id", "gene_symbol"),
      sep     = gene_id_sep,
      remove  = TRUE,
      extra   = "merge",
      fill    = "right"
    ) %>%
    dplyr::mutate(cluster = paste0("cluster_", .data$cluster))
  
  markers_up_list <- markers_df %>%
    dplyr::filter(.data$avg_log2FC >  logfc_up_thresh) %>%
    dplyr::group_by(cluster) %>%
    dplyr::group_split()
  names(markers_up_list) <- purrr::map_chr(markers_up_list, ~ unique(.x$cluster))
  
  markers_down_list <- markers_df %>%
    dplyr::filter(.data$avg_log2FC <  logfc_down_thresh) %>%
    dplyr::group_by(cluster) %>%
    dplyr::group_split()
  names(markers_down_list) <- purrr::map_chr(markers_down_list, ~ unique(.x$cluster))
  
  clusters <- sort(unique(c(names(markers_up_list), names(markers_down_list))))
  
  # [3/7] FeaturePlots
  if (verbose) message("[3/7] Preparing FeaturePlots for top ", top_n_featureplots, " genes per cluster ...")
  top_markersUp_plots_list <- purrr::map(markers_up_list, function(df_cl) {
    genes <- df_cl %>%
      dplyr::slice_head(n = top_n_featureplots) %>%
      dplyr::pull(name_id) %>%
      unique()
    if (length(genes) == 0) return(NULL)
    (\(genes_vec) Seurat::FeaturePlot(
      seurat_obj,
      features = genes_vec,
      pt.size  = featureplot_pt_size
    ))(genes)
  })
  
  # [4/7] Enrichr
  if (verbose) message("[4/7] Running Enrichr for databases: ",
                       paste(enrichr_databases, collapse = ", "), " ...")
  enrichrResults_perDatabase_perCluster <- lapply(
    enrichr_databases,
    function(db) {
      run_enrichr_for_all_clusters(
        markers_list = markers_up_list,
        database     = db,
        gene_col     = "gene_symbol"
      )
    }
  )
  names(enrichrResults_perDatabase_perCluster) <- enrichr_databases
  
  # [5/7] Enrichr plots
  if (verbose) message("[5/7] Building Enrichr plots ...")
  plots_perCluster_perDatabase <- lapply(clusters, function(cl) {
    available_dbs <- enrichr_databases[
      purrr::map_lgl(enrichr_databases, ~ !is.null(enrichrResults_perDatabase_perCluster[[.x]][[cl]]))
    ]
    if (!length(available_dbs)) return(NULL)
    
    patchwork::wrap_plots(
      plotlist = lapply(available_dbs, function(db) {
        barplot_enrichr_results(
          df                 = enrichrResults_perDatabase_perCluster[[db]][[cl]],
          min_genes          = enr_min_genes,
          n_top              = enr_n_top,
          title              = paste("Cluster", cl, "–", db),
          axis_text_size     = axis_text_size,
          axis_title_size    = axis_title_size,
          legend_text_size   = legend_text_size,
          legend_title_size  = legend_title_size,
          plot_title_size    = plot_title_size,
          label_text_size    = label_text_size,
          fill_limits        = enr_fill_limits,
          squish_to_limits   = enr_squish,
          symmetric_limits   = enr_symmetric_limits
        )
      }),
      ncol = 1
    ) +
      patchwork::plot_layout(guides = "collect") &
      ggplot2::theme(legend.position = "bottom")
  })
  names(plots_perCluster_perDatabase) <- clusters
  
  # [6/7] Save XLSX
  if (verbose) message("[6/7] Writing XLSX to: ", xlsx_path)
  write_cluster_listAndTwoPlots_to_xlsx(
    cluster_list = markers_up_list,
    plot_list_a  = top_markersUp_plots_list,
    plot_list_b  = plots_perCluster_perDatabase,
    file_path    = xlsx_path,
    img_width_a  = 11, img_height_a = 9,
    img_width_b  = 12, img_height_b = 13,
    gap_rows     = 45
  )
  
  # [7/7] Build per-cluster summaries
  if (verbose) message("[7/7] Building per-cluster summaries ...")
  # cluster_summaries <- purrr::map(clusters, function(cl) {
  #   list(
  #     markers_up_df    = markers_up_list[[cl]]   %||% NULL,
  #     markers_down_df  = markers_down_list[[cl]] %||% NULL,
  #     feature_plot_top = top_markersUp_plots_list[[cl]] %||% NULL,
  #     enrichr_plot     = plots_perCluster_perDatabase[[cl]] %||% NULL
  #   )
  # })
  # names(cluster_summaries) <- clusters
  
  cluster_summaries <- purrr::map(clusters, function(cl) {
    list(
      markers_up_df    = if (cl %in% names(markers_up_list)) markers_up_list[[cl]] else NULL,
      markers_down_df  = if (cl %in% names(markers_down_list)) markers_down_list[[cl]] else NULL,
      feature_plot_top = if (cl %in% names(top_markersUp_plots_list)) top_markersUp_plots_list[[cl]] else NULL,
      enrichr_plot     = if (cl %in% names(plots_perCluster_perDatabase)) plots_perCluster_perDatabase[[cl]] else NULL
    )
  })
  names(cluster_summaries) <- clusters
  
  if (verbose) message("✅ Done.")
  
  invisible(list(
    markers_up_list    = markers_up_list,
    markers_down_list  = markers_down_list,
    enrichr_results    = enrichrResults_perDatabase_perCluster,
    feature_plots      = top_markersUp_plots_list,
    enrichr_plots      = plots_perCluster_perDatabase,
    xlsx_path          = xlsx_path,
    cluster_summaries  = cluster_summaries
  ))
}



analyze_markers_enrichr_to_xlsx <- function(
    seurat_obj,
    xlsx_path,
    cluster_col           = "seurat_clusters",
    assay                 = NULL,
    min.pct               = 0.5,
    test.use              = "wilcox",
    p_adj_thresh          = 0.05,
    logfc_thresh          = 0.5,   # jeden próg, symetryczny
    gene_id_sep           = "-",
    top_n_featureplots    = 9,
    featureplot_pt_size   = 0.5,
    enrichr_databases     = c("CellMarker_2024", "Human_Gene_Atlas"),
    enrichr_direction     = c("up", "down", "both"),  # NOWOŚĆ
    enr_min_genes         = 2,
    enr_n_top             = 10,
    enr_fill_limits       = c(0, 10),
    enr_squish            = TRUE,
    enr_symmetric_limits  = FALSE,
    axis_text_size        = 12,
    axis_title_size       = 12,
    legend_text_size      = 10,
    legend_title_size     = 12,
    plot_title_size       = 14,
    label_text_size       = 4,
    seed                  = 0,
    verbose               = TRUE
) {
  enrichr_direction <- match.arg(enrichr_direction, c("up", "down", "both"))
  
  if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
    stop(sprintf("Column '%s' not found in seurat_obj@meta.data.", cluster_col))
  }
  if (is.null(xlsx_path) || !nzchar(xlsx_path)) {
    stop("Provide a valid 'xlsx_path'.")
  }
  
  set.seed(seed)
  
  # Zachowanie stanu wejściowego
  original_assay  <- Seurat::DefaultAssay(seurat_obj)
  original_idents <- Seurat::Idents(seurat_obj)
  on.exit({
    try({
      Seurat::DefaultAssay(seurat_obj) <- original_assay
      if (!is.null(original_idents)) {
        Seurat::Idents(seurat_obj) <- original_idents
      }
    }, silent = TRUE)
  }, add = TRUE)
  
  # Ustaw assay
  if (!is.null(assay)) {
    if (!assay %in% names(seurat_obj@assays)) {
      stop(sprintf("Assay '%s' not found in seurat_obj.", assay))
    }
    Seurat::DefaultAssay(seurat_obj) <- assay
  }
  
  # Ustaw Idents
  Seurat::Idents(seurat_obj) <- seurat_obj@meta.data[[cluster_col]]
  
  # [1/7] FindAllMarkers (always both directions)
  if (verbose) message("[1/7] Running FindAllMarkers for '", cluster_col, "' ...")
  markers_df <- Seurat::FindAllMarkers(
    object   = seurat_obj,
    only.pos = FALSE,
    min.pct  = min.pct,
    test.use = test.use
  )
  
  # [2/7] Prepare marker lists
  if (verbose) message("[2/7] Preparing marker lists (up/down split) ...")
  markers_df <- markers_df %>%
    tibble::rownames_to_column(var = "rowname") %>%
    dplyr::select(-rowname) %>%
    dplyr::filter(.data$p_val_adj < p_adj_thresh) %>%
    dplyr::mutate(name_id = .data$gene) %>%
    dplyr::select(name_id, dplyr::everything()) %>%
    tidyr::separate(
      col     = "gene",
      into    = c("gene_id", "gene_symbol"),
      sep     = gene_id_sep,
      remove  = TRUE,
      extra   = "merge",
      fill    = "right"
    ) %>%
    dplyr::mutate(cluster = paste0("cluster_", .data$cluster))
  
  markers_up_list <- markers_df %>%
    dplyr::filter(.data$avg_log2FC > logfc_thresh) %>%
    dplyr::group_by(cluster) %>%
    dplyr::group_split()
  names(markers_up_list) <- purrr::map_chr(markers_up_list, ~ unique(.x$cluster))
  
  markers_down_list <- markers_df %>%
    dplyr::filter(.data$avg_log2FC < -logfc_thresh) %>%
    dplyr::group_by(cluster) %>%
    dplyr::group_split()
  names(markers_down_list) <- purrr::map_chr(markers_down_list, ~ unique(.x$cluster))
  
  clusters <- sort(unique(c(names(markers_up_list), names(markers_down_list))))
  
  # Wybór listy do Enrichr i FeaturePlot
  selected_list <- switch(
    enrichr_direction,
    "up"   = markers_up_list,
    "down" = markers_down_list,
    "both" = {
      # Łączymy up i down
      purrr::map(clusters, function(cl) {
        dplyr::bind_rows(
          markers_up_list[[cl]]   %||% tibble(),
          markers_down_list[[cl]] %||% tibble()
        )
      }) %>% rlang::set_names(clusters)
    }
  )
  
  # [3/7] FeaturePlots (top N z wybranego kierunku)
  if (verbose) message("[3/7] Preparing FeaturePlots (direction: ", enrichr_direction, ") ...")
  top_markers_plots_list <- purrr::map(selected_list, function(df_cl) {
    if (is.null(df_cl) || nrow(df_cl) == 0) return(NULL)
    genes <- df_cl %>%
      dplyr::slice_head(n = top_n_featureplots) %>%
      dplyr::pull(name_id) %>%
      unique()
    if (length(genes) == 0) return(NULL)
    Seurat::FeaturePlot(
      seurat_obj,
      features = genes,
      pt.size  = featureplot_pt_size
    )
  })
  
  # [4/7] Enrichr
  if (verbose) message("[4/7] Running Enrichr (direction: ", enrichr_direction, ") ...")
  enrichrResults_perDatabase_perCluster <- lapply(
    enrichr_databases,
    function(db) {
      run_enrichr_for_all_clusters(
        markers_list = selected_list,
        database     = db,
        gene_col     = "gene_symbol"
      )
    }
  )
  names(enrichrResults_perDatabase_perCluster) <- enrichr_databases
  
  # [5/7] Enrichr plots
  if (verbose) message("[5/7] Building Enrichr plots ...")
  plots_perCluster_perDatabase <- lapply(clusters, function(cl) {
    available_dbs <- enrichr_databases[
      purrr::map_lgl(enrichr_databases, ~ !is.null(enrichrResults_perDatabase_perCluster[[.x]][[cl]]))
    ]
    if (!length(available_dbs)) return(NULL)
    patchwork::wrap_plots(
      plotlist = lapply(available_dbs, function(db) {
        barplot_enrichr_results(
          df                 = enrichrResults_perDatabase_perCluster[[db]][[cl]],
          min_genes          = enr_min_genes,
          n_top              = enr_n_top,
          title              = paste("Cluster", cl, "–", db),
          axis_text_size     = axis_text_size,
          axis_title_size    = axis_title_size,
          legend_text_size   = legend_text_size,
          legend_title_size  = legend_title_size,
          plot_title_size    = plot_title_size,
          label_text_size    = label_text_size,
          fill_limits        = enr_fill_limits,
          squish_to_limits   = enr_squish,
          symmetric_limits   = enr_symmetric_limits
        )
      }),
      ncol = 1
    ) +
      patchwork::plot_layout(guides = "collect") &
      ggplot2::theme(legend.position = "bottom")
  })
  names(plots_perCluster_perDatabase) <- clusters
  
  # [6/7] Save XLSX
  if (verbose) message("[6/7] Writing XLSX to: ", xlsx_path)
  write_cluster_listAndTwoPlots_to_xlsx(
    cluster_list = selected_list,
    plot_list_a  = top_markers_plots_list,
    plot_list_b  = plots_perCluster_perDatabase,
    file_path    = xlsx_path,
    img_width_a  = 11, img_height_a = 9,
    img_width_b  = 12, img_height_b = 13,
    gap_rows     = 45
  )
  
  # [7/7] Build per-cluster summaries
  if (verbose) message("[7/7] Building per-cluster summaries ...")
  # cluster_summaries <- purrr::map(clusters, function(cl) {
  #   list(
  #     markers_up_df    = markers_up_list[[cl]]   %||% NULL,
  #     markers_down_df  = markers_down_list[[cl]] %||% NULL,
  #     feature_plot     = top_markers_plots_list[[cl]] %||% NULL,
  #     enrichr_plot     = plots_perCluster_perDatabase[[cl]] %||% NULL
  #   )
  # })
  # names(cluster_summaries) <- clusters
  
  if (verbose) message("[7/7] Building per-cluster summaries ...")
  cluster_summaries <- purrr::map(clusters, function(cl) {
    list(
      markers_up_df    = if (cl %in% names(markers_up_list)) markers_up_list[[cl]] else NULL,
      markers_down_df  = if (cl %in% names(markers_down_list)) markers_down_list[[cl]] else NULL,
      feature_plot     = if (cl %in% names(top_markers_plots_list)) top_markers_plots_list[[cl]] else NULL,
      enrichr_plot     = if (cl %in% names(plots_perCluster_perDatabase)) plots_perCluster_perDatabase[[cl]] else NULL
    )
  })
  names(cluster_summaries) <- clusters
  
  if (verbose) message("✅ Done.")
  
  invisible(list(
    markers_up_list    = markers_up_list,
    markers_down_list  = markers_down_list,
    enrichr_results    = enrichrResults_perDatabase_perCluster,
    feature_plots      = top_markers_plots_list,
    enrichr_plots      = plots_perCluster_perDatabase,
    xlsx_path          = xlsx_path,
    cluster_summaries  = cluster_summaries
  ))
}
