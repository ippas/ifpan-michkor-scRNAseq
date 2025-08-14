plot_umap_cluster_faceted <- function(
    seurat_data,
    facet_var = "sample_label",
    plot_title = "UMAP faceted plot",
    facet_order = NULL,
    facet_nrow = NULL,
    facet_ncol = NULL,
    label_clusters = TRUE,
    text_size_axis_title = 10,
    text_size_axis_text = 9,
    text_size_legend_title = 10,
    text_size_legend_text = 9,
    text_size_facet = 10,
    text_size_plot_title = 12
) {
  stopifnot("umap" %in% names(seurat_data@reductions))
  
  # 1. UMAP + metadane
  umap_df <- as.data.frame(Embeddings(seurat_data, reduction = "umap")) %>%
    rownames_to_column("cell_id")
  
  metadata_df <- seurat_data@meta.data %>%
    rownames_to_column("cell_id")
  
  # Sprawdź, czy facet_var istnieje
  if (!(facet_var %in% colnames(metadata_df))) {
    stop(paste0("Column '", facet_var, "' not found in Seurat object's metadata."))
  }
  
  # 2. Połączenie i przygotowanie danych
  umap_df <- umap_df %>%
    left_join(metadata_df, by = "cell_id") %>%
    mutate(
      cluster = factor(seurat_clusters),
      facet_group = .data[[facet_var]]
    )
  
  # 3. Kolejność facetów
  if (!is.null(facet_order)) {
    umap_df$facet_group <- factor(umap_df$facet_group, levels = facet_order)
  } else {
    umap_df$facet_group <- factor(umap_df$facet_group)
  }
  
  # 4. Pozycje etykiet klastrów
  if (label_clusters) {
    label_df <- umap_df %>%
      group_by(cluster, facet_group) %>%
      summarise(
        x = median(umap_1),
        y = median(umap_2),
        .groups = "drop"
      )
  }
  
  # 5. Wykres
  p <- ggplot(umap_df, aes(x = umap_1, y = umap_2, color = cluster)) +
    geom_point(size = 0.5, alpha = 0.7)
  
  if (label_clusters) {
    p <- p + geom_text(
      data = label_df,
      aes(x = x, y = y, label = cluster),
      inherit.aes = FALSE,
      size = 3,
      color = "black"
    )
  }
  
  p <- p +
    labs(
      title = plot_title,
      color = "Cluster"
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(size = text_size_legend_title, face = "bold"),
      legend.text = element_text(size = text_size_legend_text),
      axis.title = element_text(size = text_size_axis_title),
      axis.text = element_text(size = text_size_axis_text),
      strip.text = element_text(size = text_size_facet),
      plot.title = element_text(size = text_size_plot_title, hjust = 0.5, face = "bold")
    ) +
    guides(
      color = guide_legend(nrow = 1, override.aes = list(size = 3))
    )
  
  # 6. Facetowanie
  facet_formula <- as.formula("~ facet_group")
  if (!is.null(facet_nrow) || !is.null(facet_ncol)) {
    p <- p + facet_wrap(facet_formula, nrow = facet_nrow, ncol = facet_ncol)
  } else {
    p <- p + facet_wrap(facet_formula)
  }
  
  return(p)
}



plot_umap_cluster_faceted <- function(
    seurat_data,
    facet_var = "sample_label",
    plot_title = "UMAP faceted plot",
    facet_order = NULL,
    facet_rename = NULL,
    facet_nrow = NULL,
    facet_ncol = NULL,
    label_clusters = TRUE,
    text_size_axis_title = 10,
    text_size_axis_text = 9,
    text_size_legend_title = 10,
    text_size_legend_text = 9,
    text_size_facet = 10,
    text_size_plot_title = 12
) {
  stopifnot("umap" %in% names(seurat_data@reductions))
  
  umap_df <- as.data.frame(Embeddings(seurat_data, reduction = "umap")) %>%
    rownames_to_column("cell_id")
  
  metadata_df <- seurat_data@meta.data %>%
    rownames_to_column("cell_id")
  
  if (!(facet_var %in% colnames(metadata_df))) {
    stop(paste0("Column '", facet_var, "' not found in metadata"))
  }
  
  umap_df <- umap_df %>%
    left_join(metadata_df, by = "cell_id") %>%
    mutate(
      cluster = factor(seurat_clusters),
      facet_group = .data[[facet_var]]
    )
  
  if (!is.null(facet_order)) {
    umap_df$facet_group <- factor(umap_df$facet_group, levels = facet_order)
  } else {
    umap_df$facet_group <- factor(umap_df$facet_group)
  }
  
  # Podmiana nazw facetów (tylko tych z facet_rename)
  if (!is.null(facet_rename)) {
    # zamieniamy tylko te które są w nazwach
    levels(umap_df$facet_group) <- sapply(levels(umap_df$facet_group), function(lv) {
      if (lv %in% names(facet_rename)) facet_rename[[lv]] else lv
    })
  }
  
  if (label_clusters) {
    label_df <- umap_df %>%
      group_by(cluster, facet_group) %>%
      summarise(
        x = median(umap_1),
        y = median(umap_2),
        .groups = "drop"
      )
  }
  
  p <- ggplot(umap_df, aes(x = umap_1, y = umap_2, color = cluster)) +
    geom_point(size = 0.5, alpha = 0.7)
  
  if (label_clusters) {
    p <- p + geom_text(
      data = label_df,
      aes(x = x, y = y, label = cluster),
      inherit.aes = FALSE,
      size = 3,
      color = "black"
    )
  }
  
  p <- p +
    labs(
      title = plot_title,
      color = "Cluster"
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(size = text_size_legend_title, face = "bold", color = "black"),
      legend.text = element_text(size = text_size_legend_text, color = "black"),
      axis.title = element_text(size = text_size_axis_title, color = "black"),
      axis.text = element_text(size = text_size_axis_text, color = "black"),
      strip.text = element_text(size = text_size_facet, color = "black"),
      plot.title = element_text(size = text_size_plot_title, hjust = 0.5, face = "bold", color = "black")
    ) +
    guides(
      color = guide_legend(nrow = 1, override.aes = list(size = 3))
    )
  
  facet_formula <- as.formula("~ facet_group")
  if (!is.null(facet_nrow) || !is.null(facet_ncol)) {
    p <- p + facet_wrap(facet_formula, nrow = facet_nrow, ncol = facet_ncol)
  } else {
    p <- p + facet_wrap(facet_formula)
  }
  
  return(p)
}

plot_umap_cluster_faceted <- function(
    seurat_data,
    facet_var = "sample_label",
    plot_title = "UMAP faceted plot",
    facet_order = NULL,
    facet_rename = NULL,
    facet_nrow = NULL,
    facet_ncol = NULL,
    label_clusters = TRUE,
    text_size_label = 3,
    text_size_axis_title = 10,
    text_size_axis_text = 9,
    text_size_legend_title = 10,
    text_size_legend_text = 9,
    text_size_facet = 10,
    text_size_plot_title = 12
) {
  stopifnot("umap" %in% names(seurat_data@reductions))
  
  umap_df <- as.data.frame(Embeddings(seurat_data, reduction = "umap")) %>%
    rownames_to_column("cell_id")
  
  metadata_df <- seurat_data@meta.data %>%
    rownames_to_column("cell_id")
  
  if (!(facet_var %in% colnames(metadata_df))) {
    stop(paste0("Column '", facet_var, "' not found in metadata"))
  }
  
  umap_df <- umap_df %>%
    left_join(metadata_df, by = "cell_id") %>%
    mutate(
      cluster = factor(seurat_clusters),
      facet_group = .data[[facet_var]]
    )
  
  if (!is.null(facet_order)) {
    umap_df$facet_group <- factor(umap_df$facet_group, levels = facet_order)
  } else {
    umap_df$facet_group <- factor(umap_df$facet_group)
  }
  
  if (!is.null(facet_rename)) {
    levels(umap_df$facet_group) <- sapply(levels(umap_df$facet_group), function(lv) {
      if (lv %in% names(facet_rename)) facet_rename[[lv]] else lv
    })
  }
  
  if (label_clusters) {
    label_df <- umap_df %>%
      group_by(cluster, facet_group) %>%
      summarise(
        x = median(umap_1),
        y = median(umap_2),
        .groups = "drop"
      )
  }
  
  p <- ggplot(umap_df, aes(x = umap_1, y = umap_2, color = cluster)) +
    geom_point(size = 0.5, alpha = 0.7)
  
  if (label_clusters) {
    p <- p + geom_text(
      data = label_df,
      aes(x = x, y = y, label = cluster),
      inherit.aes = FALSE,
      size = text_size_label,
      color = "black"
    )
  }
  
  p <- p +
    labs(
      title = plot_title,
      color = "Cluster"
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(size = text_size_legend_title, face = "bold", color = "black"),
      legend.text = element_text(size = text_size_legend_text, color = "black"),
      axis.title = element_text(size = text_size_axis_title, color = "black"),
      axis.text = element_text(size = text_size_axis_text, color = "black"),
      strip.text = element_text(size = text_size_facet, color = "black"),
      plot.title = element_text(size = text_size_plot_title, hjust = 0.5, face = "bold", color = "black")
    ) +
    guides(
      color = guide_legend(nrow = 1, override.aes = list(size = 3))
    )
  
  facet_formula <- as.formula("~ facet_group")
  if (!is.null(facet_nrow) || !is.null(facet_ncol)) {
    p <- p + facet_wrap(facet_formula, nrow = facet_nrow, ncol = facet_ncol)
  } else {
    p <- p + facet_wrap(facet_formula)
  }
  
  return(p)
}


plot_umap_cluster_faceted <- function(
    seurat_data,
    facet_var = "sample_label",
    plot_title = "UMAP faceted plot",
    facet_order = NULL,
    facet_rename = NULL,
    facet_nrow = NULL,
    facet_ncol = NULL,
    label_clusters = TRUE,
    text_size_label = 3,
    text_size_axis_title = 10,
    text_size_axis_text = 9,
    text_size_legend_title = 10,
    text_size_legend_text = 9,
    text_size_facet = 10,
    text_size_plot_title = 12,
    n_row_legend = 1   # nowy argument dla liczby rzędów legendy
) {
  stopifnot("umap" %in% names(seurat_data@reductions))
  
  umap_df <- as.data.frame(Embeddings(seurat_data, reduction = "umap")) %>%
    rownames_to_column("cell_id")
  
  metadata_df <- seurat_data@meta.data %>%
    rownames_to_column("cell_id")
  
  if (!(facet_var %in% colnames(metadata_df))) {
    stop(paste0("Column '", facet_var, "' not found in metadata"))
  }
  
  umap_df <- umap_df %>%
    left_join(metadata_df, by = "cell_id") %>%
    mutate(
      cluster = factor(seurat_clusters),
      facet_group = .data[[facet_var]]
    )
  
  if (!is.null(facet_order)) {
    umap_df$facet_group <- factor(umap_df$facet_group, levels = facet_order)
  } else {
    umap_df$facet_group <- factor(umap_df$facet_group)
  }
  
  if (!is.null(facet_rename)) {
    levels(umap_df$facet_group) <- sapply(levels(umap_df$facet_group), function(lv) {
      if (lv %in% names(facet_rename)) facet_rename[[lv]] else lv
    })
  }
  
  if (label_clusters) {
    label_df <- umap_df %>%
      group_by(cluster, facet_group) %>%
      summarise(
        x = median(umap_1),
        y = median(umap_2),
        .groups = "drop"
      )
  }
  
  p <- ggplot(umap_df, aes(x = umap_1, y = umap_2, color = cluster)) +
    geom_point(size = 0.5, alpha = 0.7)
  
  if (label_clusters) {
    p <- p + geom_text(
      data = label_df,
      aes(x = x, y = y, label = cluster),
      inherit.aes = FALSE,
      size = text_size_label,
      color = "black"
    )
  }
  
  p <- p +
    labs(
      title = plot_title,
      color = "Cluster"
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(size = text_size_legend_title, face = "bold", color = "black"),
      legend.text = element_text(size = text_size_legend_text, color = "black"),
      axis.title = element_text(size = text_size_axis_title, color = "black"),
      axis.text = element_text(size = text_size_axis_text, color = "black"),
      strip.text = element_text(size = text_size_facet, color = "black"),
      plot.title = element_text(size = text_size_plot_title, hjust = 0.5, face = "bold", color = "black")
    ) +
    guides(
      color = guide_legend(nrow = n_row_legend, override.aes = list(size = 3))
    )
  
  facet_formula <- as.formula("~ facet_group")
  if (!is.null(facet_nrow) || !is.null(facet_ncol)) {
    p <- p + facet_wrap(facet_formula, nrow = facet_nrow, ncol = facet_ncol)
  } else {
    p <- p + facet_wrap(facet_formula)
  }
  
  return(p)
}



plot_umap_cluster_faceted <- function(
    seurat_data,
    facet_var = "sample_label",       # nazwa kolumny do facetowania lub NULL dla jednego wykresu
    plot_title = "UMAP faceted plot",
    facet_order = NULL,
    facet_rename = NULL,
    facet_nrow = NULL,
    facet_ncol = NULL,
    label_clusters = TRUE,
    text_size_label = 3,
    text_size_axis_title = 10,
    text_size_axis_text = 9,
    text_size_legend_title = 10,
    text_size_legend_text = 9,
    text_size_facet = 10,
    text_size_plot_title = 12,
    n_row_legend = 1
) {
  stopifnot("umap" %in% names(seurat_data@reductions))
  
  umap_df <- as.data.frame(Embeddings(seurat_data, reduction = "umap")) %>%
    rownames_to_column("cell_id")
  
  metadata_df <- seurat_data@meta.data %>%
    rownames_to_column("cell_id")
  
  facet_group_present <- TRUE
  if (is.null(facet_var)) {
    facet_group_present <- FALSE
  } else if (!(facet_var %in% colnames(metadata_df))) {
    stop(paste0("Column '", facet_var, "' not found in metadata"))
  }
  
  if (facet_group_present) {
    umap_df <- umap_df %>%
      left_join(metadata_df, by = "cell_id") %>%
      mutate(
        cluster = factor(seurat_clusters),
        facet_group = .data[[facet_var]]
      )
    
    if (!is.null(facet_order)) {
      umap_df$facet_group <- factor(umap_df$facet_group, levels = facet_order)
    } else {
      umap_df$facet_group <- factor(umap_df$facet_group)
    }
    
    if (!is.null(facet_rename)) {
      levels(umap_df$facet_group) <- sapply(levels(umap_df$facet_group), function(lv) {
        if (lv %in% names(facet_rename)) facet_rename[[lv]] else lv
      })
    }
    
    if (label_clusters) {
      label_df <- umap_df %>%
        group_by(cluster, facet_group) %>%
        summarise(
          x = median(umap_1),
          y = median(umap_2),
          .groups = "drop"
        )
    }
    
  } else {
    # bez facetingu - cały dataset
    umap_df <- umap_df %>%
      left_join(metadata_df, by = "cell_id") %>%
      mutate(cluster = factor(seurat_clusters))
    
    if (label_clusters) {
      label_df <- umap_df %>%
        group_by(cluster) %>%
        summarise(
          x = median(umap_1),
          y = median(umap_2),
          .groups = "drop"
        )
    }
  }
  
  # Rysowanie wykresu
  p <- ggplot(umap_df, aes(x = umap_1, y = umap_2, color = cluster)) +
    geom_point(size = 0.5, alpha = 0.7)
  
  if (label_clusters) {
    p <- p + geom_text(
      data = label_df,
      aes(x = x, y = y, label = cluster),
      inherit.aes = FALSE,
      size = text_size_label,
      color = "black"
    )
  }
  
  p <- p +
    labs(
      title = plot_title,
      color = "Cluster"
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(size = text_size_legend_title, face = "bold", color = "black"),
      legend.text = element_text(size = text_size_legend_text, color = "black"),
      axis.title = element_text(size = text_size_axis_title, color = "black"),
      axis.text = element_text(size = text_size_axis_text, color = "black"),
      strip.text = element_text(size = text_size_facet, color = "black"),
      plot.title = element_text(size = text_size_plot_title, hjust = 0.5, face = "bold", color = "black")
    ) +
    guides(
      color = guide_legend(nrow = n_row_legend, override.aes = list(size = 3))
    )
  
  if (facet_group_present) {
    facet_formula <- as.formula("~ facet_group")
    if (!is.null(facet_nrow) || !is.null(facet_ncol)) {
      p <- p + facet_wrap(facet_formula, nrow = facet_nrow, ncol = facet_ncol)
    } else {
      p <- p + facet_wrap(facet_formula)
    }
  }
  
  return(p)
}


plot_umap_cluster_faceted <- function(
    seurat_data,
    facet_var = "sample_label",       # nazwa kolumny do facetowania lub NULL dla jednego wykresu
    cluster_var = "seurat_clusters",  # nazwa kolumny z klastrami
    plot_title = "UMAP faceted plot",
    facet_order = NULL,
    facet_rename = NULL,
    facet_nrow = NULL,
    facet_ncol = NULL,
    label_clusters = TRUE,
    text_size_label = 3,
    text_size_axis_title = 10,
    text_size_axis_text = 9,
    text_size_legend_title = 10,
    text_size_legend_text = 9,
    text_size_facet = 10,
    text_size_plot_title = 12,
    n_row_legend = 1
) {
  stopifnot("umap" %in% names(seurat_data@reductions))
  
  umap_df <- as.data.frame(Embeddings(seurat_data, reduction = "umap")) %>%
    rownames_to_column("cell_id")
  
  metadata_df <- seurat_data@meta.data %>%
    rownames_to_column("cell_id")
  
  facet_group_present <- TRUE
  if (is.null(facet_var)) {
    facet_group_present <- FALSE
  } else if (!(facet_var %in% colnames(metadata_df))) {
    stop(paste0("Column '", facet_var, "' not found in metadata"))
  }
  
  if (!(cluster_var %in% colnames(metadata_df))) {
    stop(paste0("Column '", cluster_var, "' not found in metadata"))
  }
  
  if (facet_group_present) {
    umap_df <- umap_df %>%
      left_join(metadata_df, by = "cell_id") %>%
      mutate(
        cluster = factor(.data[[cluster_var]]),
        facet_group = .data[[facet_var]]
      )
    
    if (!is.null(facet_order)) {
      umap_df$facet_group <- factor(umap_df$facet_group, levels = facet_order)
    } else {
      umap_df$facet_group <- factor(umap_df$facet_group)
    }
    
    if (!is.null(facet_rename)) {
      levels(umap_df$facet_group) <- sapply(levels(umap_df$facet_group), function(lv) {
        if (lv %in% names(facet_rename)) facet_rename[[lv]] else lv
      })
    }
    
    if (label_clusters) {
      label_df <- umap_df %>%
        group_by(cluster, facet_group) %>%
        summarise(
          x = median(umap_1),
          y = median(umap_2),
          .groups = "drop"
        )
    }
    
  } else {
    # bez facetingu - cały dataset
    umap_df <- umap_df %>%
      left_join(metadata_df, by = "cell_id") %>%
      mutate(cluster = factor(.data[[cluster_var]]))
    
    if (label_clusters) {
      label_df <- umap_df %>%
        group_by(cluster) %>%
        summarise(
          x = median(umap_1),
          y = median(umap_2),
          .groups = "drop"
        )
    }
  }
  
  # Rysowanie wykresu
  p <- ggplot(umap_df, aes(x = umap_1, y = umap_2, color = cluster)) +
    geom_point(size = 0.5, alpha = 0.7)
  
  if (label_clusters) {
    p <- p + geom_text(
      data = label_df,
      aes(x = x, y = y, label = cluster),
      inherit.aes = FALSE,
      size = text_size_label,
      color = "black"
    )
  }
  
  p <- p +
    labs(
      title = plot_title,
      color = "Cluster"
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(size = text_size_legend_title, face = "bold", color = "black"),
      legend.text = element_text(size = text_size_legend_text, color = "black"),
      axis.title = element_text(size = text_size_axis_title, color = "black"),
      axis.text = element_text(size = text_size_axis_text, color = "black"),
      strip.text = element_text(size = text_size_facet, color = "black"),
      plot.title = element_text(size = text_size_plot_title, hjust = 0.5, face = "bold", color = "black")
    ) +
    guides(
      color = guide_legend(nrow = n_row_legend, override.aes = list(size = 3))
    )
  
  if (facet_group_present) {
    facet_formula <- as.formula("~ facet_group")
    if (!is.null(facet_nrow) || !is.null(facet_ncol)) {
      p <- p + facet_wrap(facet_formula, nrow = facet_nrow, ncol = facet_ncol)
    } else {
      p <- p + facet_wrap(facet_formula)
    }
  }
  
  return(p)
}



plot_umap_cluster_faceted <- 
  function(
    seurat_data,
    facet_var = "sample_label",
    cluster_var = "seurat_clusters",
    plot_title = "UMAP faceted plot",
    facet_order = NULL,
    facet_rename = NULL,
    facet_nrow = NULL,
    facet_ncol = NULL,
    label_clusters = TRUE,
    text_size_label = 3,
    text_size_axis_title = 10,
    text_size_axis_text = 9,
    text_size_legend_title = 10,
    text_size_legend_text = 9,
    text_size_facet = 10,
    text_size_plot_title = 12,
    n_row_legend = 1
  ) {
    stopifnot("umap" %in% names(seurat_data@reductions))
    
    # Dane UMAP
    umap_df <- base::as.data.frame(Seurat::Embeddings(seurat_data, reduction = "umap")) |>
      tibble::rownames_to_column("cell_id")
    
    # Metadata
    metadata_df <- seurat_data@meta.data |>
      tibble::rownames_to_column("cell_id")
    
    facet_group_present <- TRUE
    if (base::is.null(facet_var)) {
      facet_group_present <- FALSE
    } else if (!(facet_var %in% base::colnames(metadata_df))) {
      base::stop(base::paste0("Column '", facet_var, "' not found in metadata"))
    }
    
    if (!(cluster_var %in% base::colnames(metadata_df))) {
      base::stop(base::paste0("Column '", cluster_var, "' not found in metadata"))
    }
    
    # Połączenie UMAP z metadanymi
    if (facet_group_present) {
      umap_df <- umap_df |>
        dplyr::left_join(metadata_df, by = "cell_id") |>
        dplyr::mutate(
          cluster = base::factor(!!rlang::sym(cluster_var)),
          facet_group = !!rlang::sym(facet_var)
        )
      
      # Kolejność facetów
      if (!base::is.null(facet_order)) {
        umap_df[["facet_group"]] <- base::factor(umap_df[["facet_group"]], levels = facet_order)
      } else {
        umap_df[["facet_group"]] <- base::factor(umap_df[["facet_group"]])
      }
      
      # Zmiana nazw facetów
      if (!base::is.null(facet_rename)) {
        base::levels(umap_df[["facet_group"]]) <- base::sapply(
          base::levels(umap_df[["facet_group"]]),
          function(lv) if (lv %in% base::names(facet_rename)) facet_rename[[lv]] else lv
        )
      }
      
      # Pozycje etykiet klastrów
      if (label_clusters) {
        label_df <- umap_df |>
          dplyr::group_by(cluster, facet_group) |>
          dplyr::summarise(
            x = stats::median(umap_1),
            y = stats::median(umap_2),
            .groups = "drop"
          )
      }
      
    } else {
      umap_df <- umap_df |>
        dplyr::left_join(metadata_df, by = "cell_id") |>
        dplyr::mutate(
          cluster = base::factor(!!rlang::sym(cluster_var))
        )
      
      if (label_clusters) {
        label_df <- umap_df |>
          dplyr::group_by(cluster) |>
          dplyr::summarise(
            x = stats::median(umap_1),
            y = stats::median(umap_2),
            .groups = "drop"
          )
      }
    }
    
    # Wykres UMAP
    p <- ggplot2::ggplot(umap_df, ggplot2::aes(x = umap_1, y = umap_2, color = cluster)) +
      ggplot2::geom_point(size = 0.5, alpha = 0.7)
    
    if (label_clusters) {
      p <- p + ggplot2::geom_text(
        data = label_df,
        ggplot2::aes(x = x, y = y, label = cluster),
        inherit.aes = FALSE,
        size = text_size_label,
        color = "black"
      )
    }
    
    p <- p +
      ggplot2::labs(
        title = plot_title,
        color = "Cluster"
      ) +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::theme(
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.title = ggplot2::element_text(size = text_size_legend_title, face = "bold", color = "black"),
        legend.text = ggplot2::element_text(size = text_size_legend_text, color = "black"),
        axis.title  = ggplot2::element_text(size = text_size_axis_title, color = "black"),
        axis.text   = ggplot2::element_text(size = text_size_axis_text, color = "black"),
        strip.text  = ggplot2::element_text(size = text_size_facet, color = "black"),
        plot.title  = ggplot2::element_text(size = text_size_plot_title, hjust = 0.5, face = "bold", color = "black")
      ) +
      ggplot2::guides(
        color = ggplot2::guide_legend(nrow = n_row_legend, override.aes = list(size = 3))
      )
    
    # Faceting
    if (facet_group_present) {
      facet_formula <- stats::as.formula("~ facet_group")
      if (!base::is.null(facet_nrow) || !base::is.null(facet_ncol)) {
        p <- p + ggplot2::facet_wrap(facet_formula, nrow = facet_nrow, ncol = facet_ncol)
      } else {
        p <- p + ggplot2::facet_wrap(facet_formula)
      }
    }
    
    return(p)
  }
