dotplot_log2FC_correlation_v2(
  res$cluster_2$combined %>% dplyr::rename(FDR = "g_FDR"),
  col_x = "CORT_CTRL_logFC",
  col_y = "DEX_CTRL_logFC",
  FDR_threshold   = 0.05,
  thresholds      = c(-1, 1),      # linie przerywane
  threshold_x_text = 2.5,      # etykiety, jeśli |log2FC.x| >= 1.5
  threshold_y_text = 5,         # etykiety, jeśli |log2FC.y| >= 2
  axis_text_size = 12,
  axis_title_size = 12,
  title_size = 16,
  legend_text_size = 12,
  color_sig_points = "firebrick",
  # color_sig_line   = "firebrick",
  # color_sig_labels = "black",
  label_size = 12,
  # x_limits = c(-7, 10.5),
  # y_limits = c(-6, 10.5),
  title_prefix = "cluster 2"
) 
res$cluster_2$combined %>% head


res$cluster_2$combined %>% head

res$cluster_2$combined %>% dplyr::rename(FDR = "g_FDR") %>% 
  mutate(label = ifelse(FDR < 0.05 & (CORT_CTRL_logFC > 1 | CORT_CTRL_logFC < -1 |DEX_CTRL_logFC > 1 | DEX_CTRL_logFC < -1), "signif", "ns")) %>% 
  mutate(label = ifelse( CORT_CTRL_FDR < 0.05 & DEX_CTRL_FDR > 0.05 & label == "signif", "onlyCORT", label)) 

res$cluster_2$combined %>%
  # filtr bazowy – globalny test
  filter(g_FDR < 0.05 & (abs(CORT_CTRL_logFC) > 1 | abs(DEX_CTRL_logFC) > 1)) %>%
  mutate(
    sig_CORT = CORT_CTRL_FDR < 0.05 & abs(CORT_CTRL_logFC) > 1,
    sig_DEX  = DEX_CTRL_FDR  < 0.05 & abs(DEX_CTRL_logFC)  > 1,
    label = case_when(
      sig_CORT & sig_DEX ~ "CORT_and_DEX",
      sig_CORT & !sig_DEX ~ "onlyCORT",
      !sig_CORT & sig_DEX ~ "onlyDEX",
      TRUE ~ "ns"
    )
  ) %>%
  count(label)


dotplot_log2FC_correlation_v3 <- function(
    df,
    col_x = "logFC.treatmentDEX",
    col_y = "logFC.treatmentCORT",
    method = "pearson",
    fdr_col = "FDR",
    FDR_threshold = NULL,
    thresholds = c(-1, 1),
    threshold_x_text = NULL,
    threshold_y_text = NULL,
    point_alpha = 0.6,
    point_size  = 1.5,
    label_genes = TRUE,
    # wizual
    axis_text_size = 12,
    axis_title_size = 14,
    title_size = 16,
    legend_text_size = 12,
    label_size = 12,
    x_limits = NULL,
    y_limits = NULL,
    # kolory (tylko do punktów)
    color_col = NULL,                 # np. "label"
    color_values = c(                 # mapowanie kolorów jeśli color_col != NULL
      "ns"       = "grey70",
      "signif"   = "red",
      "onlyCORT" = "dodgerblue",
      "onlyDEX"  = "darkorange",
      "both"     = "purple"
    ),
    color_labels = NULL,              # opcjonalne etykiety legendy
    legend_title = "",
    # fallback gdy brak color_col
    color_sig_points = "red",
    # tytuł
    title_prefix = ""
) {
  if (label_genes && !requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required for labeling.")
  }
  stopifnot(all(c(col_x, col_y) %in% colnames(df)))
  
  label_size_conv <- label_size / 3.88
  
  # --- Istotność (używana TYLKO do licznika n_sig i do etykiet tekstowych) ---
  df$.__signif__. <- "Not significant"
  sig_idx <- rep(FALSE, nrow(df))
  if (!is.null(FDR_threshold) && fdr_col %in% colnames(df)) {
    abs_thr <- if (length(thresholds)) max(thresholds) else 0
    sig_idx <- df[[fdr_col]] <= FDR_threshold &
      (abs(df[[col_x]]) >= abs_thr | abs(df[[col_y]]) >= abs_thr)
    df$.__signif__.[sig_idx] <- "Significant"
  }
  n_sig <- sum(sig_idx, na.rm = TRUE)
  
  # --- Korelacja globalna ---
  cor_all <- stats::cor(df[[col_x]], df[[col_y]], method = method, use = "pairwise.complete.obs")
  
  # --- Mappings do koloru PUNKTÓW ---
  use_color_col <- !is.null(color_col) && color_col %in% colnames(df)
  
  if (use_color_col) {
    # Dopasuj poziomy do zdefiniowanych kolorów (jeśli istnieją)
    present_lvls <- intersect(names(color_values), unique(as.character(df[[color_col]])))
    if (length(present_lvls) == 0) {
      # jeśli nie ma zdefiniowanych w color_values – użyj unikalnych z danych (na szaro)
      present_lvls <- unique(as.character(df[[color_col]]))
      color_values <- setNames(rep("grey50", length(present_lvls)), present_lvls)
    }
    df[[color_col]] <- factor(df[[color_col]], levels = present_lvls)
    aes_color <- ggplot2::aes(color = .data[[color_col]])
  } else {
    # fallback: koloruj po istotności
    aes_color <- ggplot2::aes(color = .data$.__signif__.)
  }
  
  # --- Wykres ---
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[col_x]], y = .data[[col_y]])) +
    ggplot2::geom_point(aes_color, alpha = point_alpha, size = point_size) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, color = "black", fullrange = TRUE) +
    ggplot2::labs(
      x = paste0("log2FC (", gsub("logFC\\.", "", col_x), ")"),
      y = paste0("log2FC (", gsub("logFC\\.", "", col_y), ")"),
      title = paste0(
        if (nzchar(title_prefix)) paste0(title_prefix, " | ") else "",
        sprintf("Correlation (%s): r = %.3f | n significant = %d", method, cor_all, n_sig)
      ),
      color = legend_title
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "bottom",
      axis.text  = ggplot2::element_text(size = axis_text_size,  color = "black"),
      axis.title = ggplot2::element_text(size = axis_title_size, color = "black"),
      plot.title = ggplot2::element_text(size = title_size, hjust = 0.5, color = "black"),
      legend.text = ggplot2::element_text(size = legend_text_size, color = "black")
    )
  
  # Skale kolorów (tylko do koloru punktów)
  if (use_color_col) {
    lvls <- levels(df[[color_col]])
    vals <- color_values[lvls]
    if (is.null(color_labels)) {
      color_labels <- setNames(lvls, lvls)
    } else {
      color_labels <- color_labels[lvls]
    }
    p <- p + ggplot2::scale_color_manual(values = vals, labels = color_labels, drop = FALSE)
  } else {
    p <- p + ggplot2::scale_color_manual(
      values = c("Not significant" = "grey70", "Significant" = color_sig_points),
      labels = c("Not significant", if (!is.null(FDR_threshold))
        sprintf("Significant (FDR < %.2f & |log2FC| ≥ %.1f)", FDR_threshold,
                if (length(thresholds)) max(thresholds) else 0) else "Significant")
    )
  }
  
  # Linie progów
  if (!is.null(thresholds) && length(thresholds) > 0) {
    for (thr in thresholds) {
      p <- p +
        ggplot2::geom_vline(xintercept = c(thr, -thr), linetype = "dashed", color = "grey50") +
        ggplot2::geom_hline(yintercept = c(thr, -thr), linetype = "dashed", color = "grey50")
    }
  }
  
  if (!is.null(x_limits)) p <- p + ggplot2::xlim(x_limits)
  if (!is.null(y_limits)) p <- p + ggplot2::ylim(y_limits)
  
  # --- Etykiety TEKSTOWE (niezależne od color_col) ---
  if (label_genes && !is.null(FDR_threshold) && "gene_symbol" %in% colnames(df)) {
    label_idx <- sig_idx &
      (( !is.null(threshold_x_text) && abs(df[[col_x]]) >= threshold_x_text ) |
         ( !is.null(threshold_y_text) && abs(df[[col_y]]) >= threshold_y_text ))
    df_labels <- subset(df, label_idx & !is.na(gene_symbol))
    if (nrow(df_labels) > 0) {
      p <- p + ggrepel::geom_text_repel(
        data = df_labels,
        ggplot2::aes(label = gene_symbol),
        size = label_size_conv,
        max.overlaps = Inf
      )
      message(sprintf("Labeling %d genes (by FDR/log2FC thresholds).", nrow(df_labels)))
    }
  }
  
  list(cor_all = cor_all, n_sig = n_sig, plot = p)
}

dotplot_log2FC_correlation_v3 <- function(
    df,
    col_x = "logFC.treatmentDEX",
    col_y = "logFC.treatmentCORT",
    method = "pearson",
    fdr_col = "FDR",
    FDR_threshold = NULL,
    thresholds = c(-1, 1),
    threshold_x_text = NULL,
    threshold_y_text = NULL,
    point_alpha = 0.6,
    point_size  = 1.5,
    label_genes = TRUE,
    # wizual
    axis_text_size = 12,
    axis_title_size = 14,
    title_size = 16,
    legend_text_size = 12,
    label_size = 12,
    x_limits = NULL,
    y_limits = NULL,
    # kolory (tylko do punktów)
    color_col = NULL,                 # np. "label"
    color_values = c(                 # mapowanie kolorów jeśli color_col != NULL
      "ns"       = "grey70",
      "signif"   = "red",
      "onlyCORT" = "dodgerblue",
      "onlyDEX"  = "darkorange",
      "both"     = "purple"
    ),
    color_labels = NULL,              # opcjonalne etykiety legendy
    legend_title = "",
    # fallback gdy brak color_col
    color_sig_points = "red",
    # tytuł
    title_prefix = ""
) {
  if (label_genes && !requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required for labeling.")
  }
  stopifnot(all(c(col_x, col_y) %in% colnames(df)))
  
  label_size_conv <- label_size / 3.88
  
  # --- Istotność (tylko do licznika n_sig i etykiet tekstowych) ---
  df$.__signif__. <- "Not significant"
  sig_idx <- rep(FALSE, nrow(df))
  if (!is.null(FDR_threshold) && fdr_col %in% colnames(df)) {
    abs_thr <- if (length(thresholds)) max(thresholds) else 0
    sig_idx <- df[[fdr_col]] <= FDR_threshold &
      (abs(df[[col_x]]) >= abs_thr | abs(df[[col_y]]) >= abs_thr)
    df$.__signif__.[sig_idx] <- "Significant"
  }
  n_sig <- sum(sig_idx, na.rm = TRUE)
  
  # --- Korelacja globalna ---
  cor_all <- stats::cor(df[[col_x]], df[[col_y]], method = method, use = "pairwise.complete.obs")
  
  # --- Mappings do koloru PUNKTÓW ---
  use_color_col <- !is.null(color_col) && color_col %in% colnames(df)
  
  if (use_color_col) {
    present_lvls <- intersect(names(color_values), unique(as.character(df[[color_col]])))
    if (length(present_lvls) == 0) {
      present_lvls <- unique(as.character(df[[color_col]]))
      color_values <- setNames(rep("grey50", length(present_lvls)), present_lvls)
    }
    df[[color_col]] <- factor(df[[color_col]], levels = present_lvls)
    aes_color <- ggplot2::aes(color = .data[[color_col]])
  } else {
    aes_color <- ggplot2::aes(color = .data$.__signif__.)
  }
  
  # --- Wykres ---
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[col_x]], y = .data[[col_y]])) +
    ggplot2::geom_point(aes_color, alpha = point_alpha, size = point_size) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, color = "black", fullrange = TRUE) +
    ggplot2::labs(
      x = paste0("log2FC (", gsub("logFC\\.", "", col_x), ")"),
      y = paste0("log2FC (", gsub("logFC\\.", "", col_y), ")"),
      title = paste0(
        if (nzchar(title_prefix)) paste0(title_prefix, " | ") else "",
        sprintf("Correlation (%s): r = %.3f | n significant = %d", method, cor_all, n_sig)
      ),
      color = legend_title
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "bottom",
      axis.text  = ggplot2::element_text(size = axis_text_size,  color = "black"),
      axis.title = ggplot2::element_text(size = axis_title_size, color = "black"),
      plot.title = ggplot2::element_text(size = title_size, hjust = 0.5, color = "black"),
      legend.text = ggplot2::element_text(size = legend_text_size, color = "black")
    )
  
  # Skale kolorów
  if (use_color_col) {
    lvls <- levels(df[[color_col]])
    vals <- color_values[lvls]
    if (is.null(color_labels)) {
      color_labels <- setNames(lvls, lvls)
    } else {
      color_labels <- color_labels[lvls]
    }
    p <- p + ggplot2::scale_color_manual(values = vals, labels = color_labels, drop = FALSE)
  } else {
    p <- p + ggplot2::scale_color_manual(
      values = c("Not significant" = "grey70", "Significant" = color_sig_points),
      labels = c("Not significant", if (!is.null(FDR_threshold))
        sprintf("Significant (FDR < %.2f & |log2FC| ≥ %.1f)", FDR_threshold,
                if (length(thresholds)) max(thresholds) else 0) else "Significant")
    )
  }
  
  # Linie progów
  if (!is.null(thresholds) && length(thresholds) > 0) {
    for (thr in thresholds) {
      p <- p +
        ggplot2::geom_vline(xintercept = c(thr, -thr), linetype = "dashed", color = "grey50") +
        ggplot2::geom_hline(yintercept = c(thr, -thr), linetype = "dashed", color = "grey50")
    }
  }
  
  if (!is.null(x_limits)) p <- p + ggplot2::xlim(x_limits)
  if (!is.null(y_limits)) p <- p + ggplot2::ylim(y_limits)
  
  # --- Etykiety TEKSTOWE ---
  if (label_genes && !is.null(FDR_threshold) && "gene_symbol" %in% colnames(df)) {
    cond_x <- if (!is.null(threshold_x_text)) abs(df[[col_x]]) >= threshold_x_text else FALSE
    cond_y <- if (!is.null(threshold_y_text)) abs(df[[col_y]]) >= threshold_y_text else FALSE
    label_idx <- sig_idx & (cond_x | cond_y)
    
    df_labels <- subset(df, label_idx & !is.na(gene_symbol))
    if (nrow(df_labels) > 0) {
      p <- p + ggrepel::geom_text_repel(
        data = df_labels,
        ggplot2::aes(label = gene_symbol),
        size = label_size_conv,
        max.overlaps = Inf
      )
      message(sprintf("Labeling %d genes (by FDR/log2FC thresholds).", nrow(df_labels)))
    }
  }
  
  list(cor_all = cor_all, n_sig = n_sig, plot = p)
}



df2 <- res$cluster_2$combined %>%
  dplyr::rename(FDR = "g_FDR") %>%
  dplyr::mutate(
    label = dplyr::if_else(
      FDR < 0.05 & (abs(CORT_CTRL_logFC) > 1 | abs(DEX_CTRL_logFC) > 1),
      "signif", "ns"
    ),
    label = dplyr::if_else(CORT_CTRL_FDR < 0.05 & DEX_CTRL_FDR > 0.05 & label == "signif", "onlyCORT", label)
  )

out <- dotplot_log2FC_correlation_v3(
  df2,
  col_x = "CORT_CTRL_logFC",
  col_y = "DEX_CTRL_logFC",
  fdr_col = "FDR",
  FDR_threshold = 0.05,
  thresholds = 1,
  color_col = "label",              # <<— tylko kolorowanie punktów wg label
  legend_title = "",
  color_values = c(                 # mapowanie kolorów jeśli color_col != NULL
    "ns"       = "grey70",
    "signif"   = "darkgreen",
    "onlyCORT" = "darkred"),
  label_genes = TRUE,          # włączamy etykiety
  threshold_x_text = 2.5,      # etykietuj geny gdy |log2FC.x| ≥ 2.5
  threshold_y_text = 5
)
out$plot

dev.off()
svg("results/tmp/figure2poster_v4.svg", width = 10, height = 7.5)
out$plot
dev.off()

