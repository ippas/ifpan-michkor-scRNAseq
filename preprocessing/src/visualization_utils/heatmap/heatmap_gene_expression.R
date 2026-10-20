heatmap_gene_expression <- function(
    mat, genes,
    quantile_normalize = FALSE,
    scale = TRUE,
    limits = c(-2, 2),
    x_text_angle = 90, x_text_size = 8, y_text_size = 6
) {
  # --- (1) opcjonalna normalizacja kwantylowa na całej macierzy ---
  if (isTRUE(quantile_normalize)) {
    rn <- rownames(mat); cn <- colnames(mat)
    mat <- as.matrix(mat)
    mat <- preprocessCore::normalize.quantiles(mat, copy = TRUE)  # wywołanie bez attach pakietu
    rownames(mat) <- rn; colnames(mat) <- cn
  }
  
  # --- (2) subset genów ---
  found <- intersect(genes, rownames(mat))
  if (length(found) == 0) stop("Żaden z podanych genów nie jest w rownames(mat).")
  sel_mat <- as.matrix(mat[found, , drop = FALSE])
  rownames(sel_mat) <- found
  
  # --- (3) skalowanie per gen ---
  if (isTRUE(scale)) {
    sel_mat <- t(scale(t(sel_mat)))
    sel_mat[is.na(sel_mat)] <- 0
  }
  
  # --- (4) przycięcie do limits ---
  sel_mat[sel_mat < limits[1]] <- limits[1]
  sel_mat[sel_mat > limits[2]] <- limits[2]
  
  # --- (5) long + rysunek ---
  df_long <- as.data.frame(sel_mat)
  df_long$gene <- rownames(sel_mat)
  df_long <- pivot_longer(df_long, -gene, names_to = "sample", values_to = "expression")
  df_long$gene   <- factor(df_long$gene, levels = rev(found))
  df_long$sample <- factor(df_long$sample, levels = colnames(sel_mat))
  
  ann_df <- df_long %>%
    dplyr::select(sample) %>%
    dplyr::distinct() %>%
    dplyr::mutate(sample = factor(sample, levels = col_order)) %>%
    dplyr::arrange(sample) %>%
    dplyr::left_join(sm, by = c("sample" = sample_col))
  
  # wymuszenie levels grupy
  ann_df[[group_col]] <- factor(
    ann_df[[group_col]],
    levels = group_levels,   # np. c("CTRL","DEX","CORT")
    ordered = TRUE
  )
  
  
  ggplot(df_long, aes(x = sample, y = gene, fill = expression)) +
    geom_tile() +
    scale_fill_gradient2(name = "expression",
                         low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = limits, oob = scales::squish) +
    theme_minimal(base_size = 11) +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(angle = x_text_angle, hjust = 1, vjust = 1, size = x_text_size),
      axis.text.y = element_text(size = y_text_size)
    )
}


# --- prosta kwantylówka w base R (bez pakietów) ---
.qn_base <- function(m) {
  m <- as.matrix(m)
  n <- nrow(m); p <- ncol(m)
  ord <- apply(m, 2, order, na.last = TRUE)           # indeksy sortowania ↑
  m_sorted <- m
  for (j in seq_len(p)) m_sorted[, j] <- m[ord[, j], j]
  avg <- rowMeans(m_sorted, na.rm = TRUE)             # średnia w każdym rzędzie po sortowaniu
  m_qn <- matrix(NA_real_, n, p, dimnames = dimnames(m))
  for (j in seq_len(p)) {
    inv <- integer(n); inv[ord[, j]] <- seq_len(n)    # odwrotność permutacji
    m_qn[, j] <- avg[inv]                             # wstaw średnie z powrotem do kolejn. oryginalnej
  }
  m_qn
}

heatmap_gene_expression <- function(
    mat, genes,
    quantile_normalize = FALSE,              # TRUE włącza QN
    quantile_method = c("base","preprocessCore"),
    scale = TRUE,
    limits = c(-2, 2),
    x_text_angle = 90, x_text_size = 8, y_text_size = 6
) {
  quantile_method <- match.arg(quantile_method)
  
  # --- (1) opcjonalna normalizacja kwantylowa na CAŁEJ macierzy ---
  if (isTRUE(quantile_normalize)) {
    rn <- rownames(mat); cn <- colnames(mat)
    mat <- as.matrix(mat)
    if (quantile_method == "base") {
      mat <- .qn_base(mat)
    } else {
      # ścieżka alternatywna, gdybyś chciał użyć preprocessCore
      if (!requireNamespace("preprocessCore", quietly = TRUE)) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
        BiocManager::install("preprocessCore")
      }
      mat <- preprocessCore::normalize.quantiles(mat, copy = TRUE)
    }
    rownames(mat) <- rn; colnames(mat) <- cn
  }
  
  # --- (2) subset genów (zachowaj kolejność wejścia) ---
  found <- intersect(genes, rownames(mat))
  if (length(found) == 0) stop("Żaden z podanych genów nie jest w rownames(mat).")
  sel_mat <- as.matrix(mat[found, , drop = FALSE])
  rownames(sel_mat) <- found
  
  # --- (3) skalowanie PER GEN (wiersz) ---
  if (isTRUE(scale)) {
    sel_mat <- t(scale(t(sel_mat)))
    sel_mat[is.na(sel_mat)] <- 0
  }
  
  # --- (4) przycięcie do limits ---
  sel_mat[sel_mat < limits[1]] <- limits[1]
  sel_mat[sel_mat > limits[2]] <- limits[2]
  
  # --- (5) long + rysunek ---
  df_long <- as.data.frame(sel_mat)
  df_long$gene <- rownames(sel_mat)
  df_long <- pivot_longer(df_long, -gene, names_to = "sample", values_to = "expression")
  df_long$gene   <- factor(df_long$gene, levels = rev(found))
  df_long$sample <- factor(df_long$sample, levels = colnames(sel_mat))
  
  ggplot(df_long, aes(x = sample, y = gene, fill = expression)) +
    geom_tile() +
    scale_fill_gradient2(name = "expression",
                         low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = limits, oob = scales::squish) +
    theme_minimal(base_size = 11) +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(angle = x_text_angle, hjust = 1, vjust = 1, size = x_text_size),
      axis.text.y = element_text(size = y_text_size)
    )
}

heatmap_gene_expression <- function(
    mat, genes,
    quantile_normalize = FALSE,
    scale = TRUE,
    limits = c(-2, 2),
    cluster_row = FALSE,
    cluster_col = FALSE,
    dist_method = "euclidean",
    hclust_method = "complete",
    x_text_angle = 90, x_text_size = 8, y_text_size = 6
) {
  # 1) opcjonalna kwantylówka na całej macierzy
  if (isTRUE(quantile_normalize)) {
    rn <- rownames(mat); cn <- colnames(mat)
    mat <- .qn_base(as.matrix(mat))
    rownames(mat) <- rn; colnames(mat) <- cn
  }
  
  # 2) subset genów (kolejność jak w wejściu)
  found <- intersect(genes, rownames(mat))
  if (length(found) == 0) stop("Żaden z podanych genów nie jest w rownames(mat).")
  sel_mat <- as.matrix(mat[found, , drop = FALSE])
  rownames(sel_mat) <- found
  
  # 3) skalowanie per gen (wiersz)
  if (isTRUE(scale)) {
    sel_mat <- t(scale(t(sel_mat)))
    sel_mat[is.na(sel_mat)] <- 0
  }
  
  # 4) klastrowanie (po przekształceniach)
  row_order <- found
  col_order <- colnames(sel_mat)
  if (isTRUE(cluster_row) && nrow(sel_mat) > 1) {
    row_order <- rownames(sel_mat)[hclust(dist(sel_mat, method = dist_method), method = hclust_method)$order]
  }
  if (isTRUE(cluster_col) && ncol(sel_mat) > 1) {
    col_order <- colnames(sel_mat)[hclust(dist(t(sel_mat), method = dist_method), method = hclust_method)$order]
  }
  
  # 5) przycięcie do limits
  sel_mat[sel_mat < limits[1]] <- limits[1]
  sel_mat[sel_mat > limits[2]] <- limits[2]
  
  # 6) long + rysunek
  df_long <- as.data.frame(sel_mat)
  df_long$gene <- rownames(sel_mat)
  df_long <- pivot_longer(df_long, -gene, names_to = "sample", values_to = "expression")
  
  df_long$gene   <- factor(df_long$gene,   levels = rev(row_order))  # od góry
  df_long$sample <- factor(df_long$sample, levels = col_order)
  
  ggplot(df_long, aes(x = sample, y = gene, fill = expression)) +
    geom_tile() +
    scale_fill_gradient2(name = "expression",
                         low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = limits, oob = scales::squish) +
    theme_minimal(base_size = 11) +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(angle = x_text_angle, hjust = 1, vjust = 1, size = x_text_size),
      axis.text.y = element_text(size = y_text_size)
    )
}

heatmap_gene_expression <- function(
    mat, genes,
    sample_meta,
    quantile_normalize = FALSE,
    scale = TRUE,
    show_dendrogram = TRUE,                # alias dla cluster_row
    limits = c(-2, 2),
    group_col = "treatment",
    group_levels = NULL,
    group_palette = NULL,
    legend_position = "top",
    annotation_height = 0.08,
    col_low = "blue",                      # nowe
    col_mid = "white",
    col_high = "red"
) {
  stopifnot(all(c("sample_id", group_col, "sample_order") %in% colnames(sample_meta)))
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(RColorBrewer)
  
  .qn_base <- function(x) {
    ox <- apply(x, 2, order, na.last = TRUE)
    sx <- apply(x, 2, sort,  na.last = TRUE)
    row_means <- rowMeans(sx, na.rm = TRUE)
    res <- x
    for (j in seq_len(ncol(x))) res[ox[, j], j] <- row_means
    res
  }
  
  genes_found <- intersect(genes, rownames(mat))
  if (length(genes_found) == 0) stop("Żaden z podanych genów nie występuje w mat.")
  
  sel_mat <- as.matrix(mat[genes_found, , drop = FALSE])
  rownames(sel_mat) <- genes_found
  
  if (quantile_normalize) {
    rn <- rownames(sel_mat); cn <- colnames(sel_mat)
    sel_mat <- .qn_base(sel_mat)
    rownames(sel_mat) <- rn; colnames(sel_mat) <- cn
  }
  
  if (scale) {
    sel_mat <- t(scale(t(sel_mat)))
    sel_mat[is.na(sel_mat)] <- 0
  }
  
  sel_mat[sel_mat < limits[1]] <- limits[1]
  sel_mat[sel_mat > limits[2]] <- limits[2]
  
  row_order <- genes_found
  if (show_dendrogram && nrow(sel_mat) > 1) {
    row_order <- rownames(sel_mat)[hclust(dist(sel_mat))$order]
  }
  
  sample_meta <- sample_meta %>%
    filter(sample_id %in% colnames(sel_mat)) %>%
    arrange(sample_order)
  
  col_order <- sample_meta$sample_id
  
  df_long <- as.data.frame(sel_mat) %>%
    tibble::rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "sample", values_to = "expression") %>%
    mutate(
      gene = factor(gene, levels = rev(row_order)),
      sample = factor(sample, levels = col_order)
    )
  
  ann_df <- sample_meta %>%
    select(sample = sample_id, group = !!sym(group_col)) %>%
    mutate(sample = factor(sample, levels = col_order))
  
  if (!is.null(group_levels)) {
    ann_df$group <- factor(ann_df$group, levels = group_levels, ordered = TRUE)
  } else {
    ann_df$group <- factor(ann_df$group, ordered = TRUE)
    group_levels <- levels(ann_df$group)
  }
  
  if (is.null(group_palette)) {
    base_cols <- brewer.pal(max(3, length(group_levels)), "Set1")[seq_along(group_levels)]
    group_palette <- setNames(base_cols, group_levels)
  } else {
    group_palette <- group_palette[group_levels]
  }
  
  # pasek treatmentów
  p_ann <- ggplot(ann_df, aes(x = sample, y = 1, fill = group)) +
    geom_tile() +
    scale_fill_manual(values = group_palette, name = group_col, drop = FALSE) +
    theme_void() +
    theme(legend.position = legend_position)
  
  # heatmapa
  p_heat <- ggplot(df_long, aes(x = sample, y = gene, fill = expression)) +
    geom_tile() +
    scale_fill_gradient2(
      name = "expression",
      low = col_low,
      mid = col_mid,
      high = col_high,
      midpoint = 0,
      limits = limits,
      oob = scales::squish
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 8),
      axis.text.y = element_text(size = 6),
      axis.title = element_blank(),
      legend.position = legend_position
    )
  
  # złożenie
  (p_ann / p_heat) +
    plot_layout(heights = c(annotation_height, 1), guides = "collect") &
    theme(legend.position = legend_position)
}


heatmap_gene_expression <- function(
    mat, genes,
    sample_meta,
    quantile_normalize = FALSE,
    scale = TRUE,
    cluster_row = TRUE,
    limits = c(-2, 2),
    group_col = "treatment",
    group_levels = NULL,
    group_palette = NULL,
    legend_position = "top",
    annotation_height = 0.08,
    col_low = "blue",
    col_mid = "white",
    col_high = "red"
) {
  stopifnot(all(c("sample_id", group_col, "sample_order") %in% colnames(sample_meta)))
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(RColorBrewer)
  
  .qn_base <- function(x) {
    ox <- apply(x, 2, order, na.last = TRUE)
    sx <- apply(x, 2, sort,  na.last = TRUE)
    row_means <- rowMeans(sx, na.rm = TRUE)
    res <- x
    for (j in seq_len(ncol(x))) res[ox[, j], j] <- row_means
    res
  }
  
  genes_found <- intersect(genes, rownames(mat))
  if (length(genes_found) == 0) stop("Żaden z podanych genów nie występuje w mat.")
  
  sel_mat <- as.matrix(mat[genes_found, , drop = FALSE])
  rownames(sel_mat) <- genes_found
  
  if (quantile_normalize) {
    rn <- rownames(sel_mat); cn <- colnames(sel_mat)
    sel_mat <- .qn_base(sel_mat)
    rownames(sel_mat) <- rn; colnames(sel_mat) <- cn
  }
  
  if (scale) {
    sel_mat <- t(scale(t(sel_mat)))
    sel_mat[is.na(sel_mat)] <- 0
  }
  
  sel_mat[sel_mat < limits[1]] <- limits[1]
  sel_mat[sel_mat > limits[2]] <- limits[2]
  
  row_order <- genes_found
  if (cluster_row && nrow(sel_mat) > 1) {
    row_order <- rownames(sel_mat)[hclust(dist(sel_mat))$order]
  }
  
  sample_meta <- sample_meta %>%
    filter(sample_id %in% colnames(sel_mat)) %>%
    arrange(sample_order)
  
  col_order <- sample_meta$sample_id
  
  df_long <- as.data.frame(sel_mat) %>%
    tibble::rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "sample", values_to = "expression") %>%
    mutate(
      gene = factor(gene, levels = rev(row_order)),
      sample = factor(sample, levels = col_order)
    )
  
  ann_df <- sample_meta %>%
    select(sample = sample_id, group = !!sym(group_col)) %>%
    mutate(sample = factor(sample, levels = col_order))
  
  if (!is.null(group_levels)) {
    ann_df$group <- factor(ann_df$group, levels = group_levels, ordered = TRUE)
  } else {
    ann_df$group <- factor(ann_df$group, ordered = TRUE)
    group_levels <- levels(ann_df$group)
  }
  
  if (is.null(group_palette)) {
    base_cols <- brewer.pal(max(3, length(group_levels)), "Set1")[seq_along(group_levels)]
    group_palette <- setNames(base_cols, group_levels)
  } else {
    group_palette <- group_palette[group_levels]
  }
  
  # pasek treatmentów
  p_ann <- ggplot(ann_df, aes(x = sample, y = 1, fill = group)) +
    geom_tile() +
    scale_fill_manual(values = group_palette, name = group_col, drop = FALSE) +
    theme_void() +
    theme(legend.position = legend_position)
  
  # heatmapa
  p_heat <- ggplot(df_long, aes(x = sample, y = gene, fill = expression)) +
    geom_tile() +
    scale_fill_gradient2(
      name = "expression",
      low = col_low,
      mid = col_mid,
      high = col_high,
      midpoint = 0,
      limits = limits,
      oob = scales::squish
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 8),
      axis.text.y = element_text(size = 6),
      axis.title = element_blank(),
      legend.position = legend_position
    )
  
  # złożenie
  (p_ann / p_heat) +
    plot_layout(heights = c(annotation_height, 1), guides = "collect") &
    theme(legend.position = legend_position)
}
