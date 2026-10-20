seurat_featurePlot_binary_clusterSummary <- function(
    seurat_obj,
    gene_pattern,
    assay = "RNA",
    slot = "counts",
    split.by = NULL,
    split_levels = NULL,
    cluster_col = NULL,
    reduction = "umap",
    cols = c("gray90", "firebrick"),
    point_size = 0.8,
    legend_position = "none",
    detection_threshold = 0,
    return_data = FALSE,
    cluster_summary_figure = NULL,   # 🔹 np. c("Astrocytes" = "Astro", "Oligodendrocytes" = "Oligo")
    title_size = 16,
    facet_title_size = 14,
    label_size = 12
) {
  requireNamespace("ggplot2")
  requireNamespace("dplyr")
  
  
  label_size <- label_size/2.84
  
  # --- lokalna kopia i ustawienie assay
  obj_tmp <- seurat_obj
  if (!is.null(assay)) DefaultAssay(obj_tmp) <- assay
  
  # --- wyszukiwanie genu
  all_features <- rownames(GetAssayData(obj_tmp, assay = assay, slot = slot))
  matched <- grep(paste0("-", gene_pattern, "$"), all_features, value = TRUE)
  if (length(matched) == 0) stop(paste0("❌ No gene found matching pattern: ", gene_pattern))
  if (length(matched) > 1) {
    message("⚠️ Multiple matches found, using first: ", paste(matched, collapse = ", "))
    matched <- matched[1]
  } else {
    message("✅ Found feature: ", matched)
  }
  
  # --- embeddings + expression
  emb <- as.data.frame(Embeddings(obj_tmp, reduction = reduction))
  colnames(emb) <- c("dim1", "dim2")[seq_len(ncol(emb))]
  expr_values <- GetAssayData(obj_tmp, assay = assay, slot = slot)[matched, ]
  detected_binary <- as.numeric(expr_values > detection_threshold)
  
  data_df <- emb |>
    dplyr::mutate(expression = expr_values, detected = detected_binary)
  
  # --- metadata
  if (!is.null(split.by)) data_df$treatment <- obj_tmp@meta.data[[split.by]]
  if (!is.null(cluster_col) && cluster_col %in% colnames(obj_tmp@meta.data))
    data_df$cluster <- obj_tmp@meta.data[[cluster_col]]
  
  # --- porządek treatment
  if (!is.null(split.by)) {
    current_values <- unique(data_df$treatment)
    if (is.null(split_levels)) {
      preferred_order <- c("CTRL", "CORT", "DEX")
      split_levels <- intersect(preferred_order, current_values)
      split_levels <- c(split_levels, setdiff(current_values, split_levels))
    }
    data_df$treatment <- factor(data_df$treatment, levels = split_levels)
  }
  
  # --- summary do facetów
  if (!is.null(split.by)) {
    summary_df <- as.data.frame(table(data_df$treatment))
    colnames(summary_df) <- c("treatment", "total_cells")
    det_summary <- as.data.frame(table(data_df$treatment[data_df$detected == 1]))
    colnames(det_summary) <- c("treatment", "detected_count")
    summary_df <- merge(summary_df, det_summary, by = "treatment", all.x = TRUE)
    summary_df$detected_count[is.na(summary_df$detected_count)] <- 0
    summary_df$percent <- round(100 * summary_df$detected_count / summary_df$total_cells, 1)
    data_df <- merge(data_df, summary_df, by = "treatment", all.x = TRUE)
    data_df$facet_label <- paste0(
      data_df$treatment, "\n",
      data_df$detected_count, "/", data_df$total_cells,
      " (", data_df$percent, "%)"
    )
    if (!is.null(split_levels)) {
      data_df$treatment <- factor(data_df$treatment, levels = split_levels)
      facet_labels_ordered <- data_df |>
        dplyr::distinct(treatment, facet_label) |>
        dplyr::arrange(match(treatment, split_levels)) |>
        dplyr::pull(facet_label)
      data_df$facet_label <- factor(data_df$facet_label, levels = facet_labels_ordered)
    }
    title_text <- paste0(matched, " (binary detection)")
  } else {
    n_detected <- sum(data_df$detected)
    n_total <- nrow(data_df)
    percent_detected <- round(100 * n_detected / n_total, 1)
    title_text <- paste0(
      matched, " — detected in ", n_detected, "/", n_total,
      " cells (", percent_detected, "%)"
    )
  }
  
  # --- summary per cluster
  cluster_summary <- NULL
  if (!is.null(cluster_col)) {
    if (!is.null(split.by)) {
      cluster_summary <- data_df |>
        dplyr::group_by(treatment, cluster) |>
        dplyr::summarise(
          n_cells_total = dplyr::n(),
          n_cells_detected = sum(detected, na.rm = TRUE),
          percent_detected = round(100 * n_cells_detected / n_cells_total, 1),
          .groups = "drop"
        )
    } else {
      cluster_summary <- data_df |>
        dplyr::group_by(cluster) |>
        dplyr::summarise(
          n_cells_total = dplyr::n(),
          n_cells_detected = sum(detected, na.rm = TRUE),
          percent_detected = round(100 * n_cells_detected / n_cells_total, 1),
          .groups = "drop"
        )
    }
    message("📊 Cluster summary:")
    print(cluster_summary)
  }
  
  if (return_data) return(list(data = data_df, cluster_summary = cluster_summary))
  
  # --- wykres
  p <- ggplot2::ggplot(data_df, ggplot2::aes(x = dim1, y = dim2, color = factor(detected))) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::scale_color_manual(values = c("0" = cols[1], "1" = cols[2])) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      legend.position = legend_position,
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = title_size, color = "black"),
      strip.text = ggplot2::element_text(size = facet_title_size, color = "black"),
      text = ggplot2::element_text(color = "black")
    ) +
    ggplot2::labs(title = title_text, color = "Detection")
  
  if (!is.null(split.by)) {
    p <- p + ggplot2::facet_wrap(~facet_label)
  }
  
  # --- 🔹 dodanie tekstów tylko dla wybranych klastrów (z podmianą nazw)
  if (!is.null(cluster_summary_figure) && !is.null(cluster_summary)) {
    # weź nazwy klastrów do filtrowania
    clusters_to_show <- names(cluster_summary_figure)
    if (is.null(clusters_to_show)) clusters_to_show <- cluster_summary_figure
    
    selected_clusters <- cluster_summary |>
      dplyr::filter(cluster %in% clusters_to_show)
    
    if (nrow(selected_clusters) > 0) {
      # podmiana nazw na skrócone
      label_map <- cluster_summary_figure
      selected_clusters$cluster_label <- dplyr::recode(selected_clusters$cluster, !!!label_map)
      
      selected_clusters <- selected_clusters |>
        dplyr::mutate(
          text_label = paste0(cluster_label, ": ", n_cells_detected, "/", n_cells_total,
                              " (", percent_detected, "%)")
        )
      
      # współrzędne dolnego rogu osobno dla każdego treatment
      coords <- data_df |>
        dplyr::group_by(treatment, facet_label) |>
        dplyr::summarise(
          x = min(dim1, na.rm = TRUE),
          y = min(dim2, na.rm = TRUE),
          .groups = "drop"
        )
      
      # przypisz facet_label do tekstów (żeby każdy panel miał tylko swoje)
      selected_clusters <- dplyr::left_join(selected_clusters, coords, by = "treatment") |>
        dplyr::mutate(facet_label = factor(facet_label, levels = unique(coords$facet_label)))
      
      text_blocks <- selected_clusters |>
        dplyr::group_by(facet_label, x, y) |>
        dplyr::summarise(
          text_label = paste(text_label, collapse = "\n"),
          .groups = "drop"
        )
      
      p <- p + ggplot2::geom_text(
        data = text_blocks,
        ggplot2::aes(x = x, y = y, label = text_label),
        hjust = 0, vjust = 0, color = "black", size = label_size,
        inherit.aes = FALSE
      )
    }
  }
  
  return(p)
}

seurat_featurePlot_binary_clusterSummary <- function(
    seurat_obj,
    gene_pattern,
    assay = "RNA",
    slot = "counts",
    split.by = NULL,
    split_levels = NULL,
    cluster_col = NULL,
    reduction = "umap",
    cols = c("gray90", "firebrick"),
    point_size = 0.8,
    legend_position = "none",
    detection_threshold = 0,
    return_data = FALSE,
    cluster_summary_figure = NULL,   # np. c("Astrocytes" = "Astro", "Oligodendrocytes" = "Oligo")
    title_size = 16,
    facet_title_size = 14,
    label_size = 12
) {
  requireNamespace("ggplot2")
  requireNamespace("dplyr")
  
  label_size <- label_size / 2.84
  
  # --- lokalna kopia i ustawienie assay
  obj_tmp <- seurat_obj
  if (!is.null(assay)) DefaultAssay(obj_tmp) <- assay
  
  # --- wyszukiwanie genu
  all_features <- rownames(GetAssayData(obj_tmp, assay = assay, slot = slot))
  matched <- grep(paste0("-", gene_pattern, "$"), all_features, value = TRUE)
  if (length(matched) == 0) stop(paste0("❌ No gene found matching pattern: ", gene_pattern))
  if (length(matched) > 1) {
    message("⚠️ Multiple matches found, using first: ", paste(matched, collapse = ", "))
    matched <- matched[1]
  } else {
    message("✅ Found feature: ", matched)
  }
  
  # --- embeddings + expression
  emb <- as.data.frame(Embeddings(obj_tmp, reduction = reduction))
  colnames(emb) <- c("dim1", "dim2")[seq_len(ncol(emb))]
  expr_values <- GetAssayData(obj_tmp, assay = assay, slot = slot)[matched, ]
  detected_binary <- as.numeric(expr_values > detection_threshold)
  
  data_df <- emb |>
    dplyr::mutate(expression = expr_values, detected = detected_binary)
  
  # --- metadata
  if (!is.null(split.by)) data_df$treatment <- obj_tmp@meta.data[[split.by]]
  if (!is.null(cluster_col) && cluster_col %in% colnames(obj_tmp@meta.data))
    data_df$cluster <- obj_tmp@meta.data[[cluster_col]]
  
  # --- porządek treatment
  if (!is.null(split.by)) {
    current_values <- unique(data_df$treatment)
    if (is.null(split_levels)) {
      preferred_order <- c("CTRL", "CORT", "DEX")
      split_levels <- intersect(preferred_order, current_values)
      split_levels <- c(split_levels, setdiff(current_values, split_levels))
    }
    data_df$treatment <- factor(data_df$treatment, levels = split_levels)
  }
  
  # --- summary do facetów
  if (!is.null(split.by)) {
    summary_df <- as.data.frame(table(data_df$treatment))
    colnames(summary_df) <- c("treatment", "total_cells")
    det_summary <- as.data.frame(table(data_df$treatment[data_df$detected == 1]))
    colnames(det_summary) <- c("treatment", "detected_count")
    summary_df <- merge(summary_df, det_summary, by = "treatment", all.x = TRUE)
    summary_df$detected_count[is.na(summary_df$detected_count)] <- 0
    summary_df$percent <- round(100 * summary_df$detected_count / summary_df$total_cells, 1)
    data_df <- merge(data_df, summary_df, by = "treatment", all.x = TRUE)
    data_df$facet_label <- paste0(
      data_df$treatment, "\n",
      data_df$detected_count, "/", data_df$total_cells,
      " (", data_df$percent, "%)"
    )
    if (!is.null(split_levels)) {
      data_df$treatment <- factor(data_df$treatment, levels = split_levels)
      facet_labels_ordered <- data_df |>
        dplyr::distinct(treatment, facet_label) |>
        dplyr::arrange(match(treatment, split_levels)) |>
        dplyr::pull(facet_label)
      data_df$facet_label <- factor(data_df$facet_label, levels = facet_labels_ordered)
    }
    title_text <- paste0(matched, " (binary detection)")
  } else {
    n_detected <- sum(data_df$detected)
    n_total <- nrow(data_df)
    percent_detected <- round(100 * n_detected / n_total, 1)
    title_text <- paste0(
      matched, " — detected in ", n_detected, "/", n_total,
      " cells (", percent_detected, "%)"
    )
  }
  
  # --- summary per cluster
  cluster_summary <- NULL
  if (!is.null(cluster_col)) {
    if (!is.null(split.by)) {
      cluster_summary <- data_df |>
        dplyr::group_by(treatment, cluster) |>
        dplyr::summarise(
          n_cells_total = dplyr::n(),
          n_cells_detected = sum(detected, na.rm = TRUE),
          percent_detected = round(100 * n_cells_detected / n_cells_total, 1),
          .groups = "drop"
        )
    } else {
      cluster_summary <- data_df |>
        dplyr::group_by(cluster) |>
        dplyr::summarise(
          n_cells_total = dplyr::n(),
          n_cells_detected = sum(detected, na.rm = TRUE),
          percent_detected = round(100 * n_cells_detected / n_cells_total, 1),
          .groups = "drop"
        )
    }
    message("📊 Cluster summary:")
    print(cluster_summary)
  }
  
  if (return_data) return(list(data = data_df, cluster_summary = cluster_summary))
  
  # --- wykres
  p <- ggplot2::ggplot(data_df, ggplot2::aes(x = dim1, y = dim2, color = factor(detected))) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::scale_color_manual(values = c("0" = cols[1], "1" = cols[2])) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      legend.position = legend_position,
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = title_size, color = "black"),
      strip.text = ggplot2::element_text(size = facet_title_size, color = "black"),
      text = ggplot2::element_text(color = "black")
    ) +
    ggplot2::labs(title = title_text, color = "Detection")
  
  if (!is.null(split.by)) {
    p <- p + ggplot2::facet_wrap(~facet_label)
  }
  
  # --- 🔹 dodanie tekstów tylko dla wybranych klastrów (z zachowaniem kolejności)
  if (!is.null(cluster_summary_figure) && !is.null(cluster_summary)) {
    clusters_to_show <- names(cluster_summary_figure)
    if (is.null(clusters_to_show)) clusters_to_show <- cluster_summary_figure
    
    selected_clusters <- cluster_summary |>
      dplyr::filter(cluster %in% clusters_to_show)
    
    if (nrow(selected_clusters) > 0) {
      label_map <- cluster_summary_figure
      selected_clusters$cluster_label <- dplyr::recode(selected_clusters$cluster, !!!label_map)
      
      selected_clusters <- selected_clusters |>
        dplyr::mutate(
          text_label = paste0(cluster_label, ": ", n_cells_detected, "/", n_cells_total,
                              " (", percent_detected, "%)")
        )
      
      # współrzędne dolnego rogu dla każdego treatment
      coords <- data_df |>
        dplyr::group_by(treatment, facet_label) |>
        dplyr::summarise(
          x = min(dim1, na.rm = TRUE),
          y = min(dim2, na.rm = TRUE),
          .groups = "drop"
        )
      
      selected_clusters <- dplyr::left_join(selected_clusters, coords, by = "treatment") |>
        dplyr::mutate(facet_label = factor(facet_label, levels = unique(coords$facet_label)))
      
      # 🔹 Ustal kolejność zgodną z cluster_summary_figure
      cluster_order <- names(cluster_summary_figure)
      if (is.null(cluster_order)) cluster_order <- cluster_summary_figure
      selected_clusters$cluster_label <- factor(selected_clusters$cluster_label,
                                                levels = cluster_summary_figure)
      
      text_blocks <- selected_clusters |>
        dplyr::arrange(factor(cluster_label, levels = cluster_summary_figure)) |>
        dplyr::group_by(facet_label, x, y) |>
        dplyr::summarise(
          text_label = paste(text_label, collapse = "\n"),
          .groups = "drop"
        )
      
      p <- p + ggplot2::geom_text(
        data = text_blocks,
        ggplot2::aes(x = x, y = y, label = text_label),
        hjust = 0, vjust = 0, color = "black", size = label_size,
        inherit.aes = FALSE
      )
    }
  }
  
  return(p)
}



# seurat_featurePlot_binary_clusterSummary(
#   seurat_obj = integrated_pilotDexCortCtrl_filtered,
#   gene_pattern = "PEG3",
#   split.by = "treatment",
#   split_levels = c("CTRL", "CORT", "DEX"),
#   cluster_col = "res0.3_RNA_noScaled_brain",
#   label_size = 12,
#   point_size = 0.4,
#   cluster_summary_figure = c("Astrocytes" = "Astrocytes", "Oligodendrocytes" = "Oligodendrocytes",
#                              "Oligodendrocyte precursor cells" = "OPCs", "Radial glial cells" = "RGCs")
# )

