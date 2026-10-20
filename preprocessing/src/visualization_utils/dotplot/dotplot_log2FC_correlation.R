dotplot_log2FC_correlation <- function(
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
    # visual parameters
    axis_text_size = 12,
    axis_title_size = 14,
    title_size = 16,
    legend_text_size = 12,
    label_size = 12,              # now in same units as axis/title (pt)
    x_limits = NULL,
    y_limits = NULL,
    # color parameters
    color_sig_points = "red",
    color_sig_line   = "red",
    color_sig_labels = "red"
) {
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required for labeling.")
  }
  stopifnot(all(c(col_x, col_y) %in% colnames(df)))
  
  # convert label size from pt to ggplot scale (geom_text_repel uses ~3.88 factor)
  label_size_conv <- label_size / 3.88
  
  # classify significant genes
  df$highlight <- "Not significant"
  if (!is.null(FDR_threshold) && fdr_col %in% colnames(df)) {
    sig_idx <- df[[fdr_col]] <= FDR_threshold &
      (abs(df[[col_x]]) >= max(thresholds) | abs(df[[col_y]]) >= max(thresholds))
    df$highlight[sig_idx] <- "Significant"
  }
  
  # correlations
  cor_all <- cor(df[[col_x]], df[[col_y]], method = method, use = "pairwise.complete.obs")
  df_sig  <- subset(df, highlight == "Significant")
  cor_sig <- if (nrow(df_sig) >= 3) {
    cor(df_sig[[col_x]], df_sig[[col_y]], method = method, use = "pairwise.complete.obs")
  } else NA
  n_sig <- nrow(df_sig)
  
  # legend labels
  legend_labels <- c(
    "Not significant" = "No significant",
    "Significant" = sprintf("Significant (FDR < %.2f & |log2FC| ≥ %.1f)", 
                            FDR_threshold, max(thresholds))
  )
  legend_labels <- gsub("0+$", "", legend_labels)  
  legend_labels <- gsub("\\.$", "", legend_labels)
  
  # scatter
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[col_x]], y = .data[[col_y]])) +
    ggplot2::geom_point(ggplot2::aes(color = highlight), alpha = point_alpha, size = point_size) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, color = "black", fullrange = TRUE) +
    {if (nrow(df_sig) >= 3) ggplot2::geom_smooth(
      data = df_sig,
      ggplot2::aes(x = .data[[col_x]], y = .data[[col_y]]),
      method = "lm", se = FALSE, color = color_sig_line, fullrange = TRUE)} +
    ggplot2::scale_color_manual(
      values = c("Not significant" = "grey70", "Significant" = color_sig_points),
      labels = legend_labels
    ) +
    ggplot2::labs(
      x = paste0("log2FC (", gsub("logFC\\.", "", col_x), ")"),
      y = paste0("log2FC (", gsub("logFC\\.", "", col_y), ")"),
      title = sprintf("Correlation (%s): all = %.3f, significant = %.3f (n = %d)", 
                      method, cor_all, cor_sig, n_sig),
      color = ""
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom",
      axis.text = ggplot2::element_text(size = axis_text_size, color = "black"),
      axis.title = ggplot2::element_text(size = axis_title_size, color = "black"),
      plot.title = ggplot2::element_text(size = title_size, hjust = 0.5, color = "black"),
      legend.text = ggplot2::element_text(size = legend_text_size, color = "black")
    )
  
  # dashed lines
  if (!is.null(thresholds) && length(thresholds) > 0) {
    for (thr in thresholds) {
      p <- p +
        ggplot2::geom_vline(xintercept = c(thr, -thr), linetype = "dashed", color = "grey50") +
        ggplot2::geom_hline(yintercept = c(thr, -thr), linetype = "dashed", color = "grey50")
    }
  }
  
  # limits
  if (!is.null(x_limits)) p <- p + ggplot2::xlim(x_limits)
  if (!is.null(y_limits)) p <- p + ggplot2::ylim(y_limits)
  
  # labels
  if (label_genes && !is.null(FDR_threshold) && "gene_symbol" %in% colnames(df)) {
    label_idx <- df[[fdr_col]] <= FDR_threshold &
      (( !is.null(threshold_x_text) & abs(df[[col_x]]) >= threshold_x_text ) |
         ( !is.null(threshold_y_text) & abs(df[[col_y]]) >= threshold_y_text ))
    
    df_labels <- subset(df, label_idx & !is.na(gene_symbol))
    if (nrow(df_labels) > 0) {
      p <- p + ggrepel::geom_text_repel(
        data = df_labels,
        ggplot2::aes(label = gene_symbol),
        size = label_size_conv,   # converted to geom_text scale
        max.overlaps = Inf,
        color = color_sig_labels
      )
      message(sprintf("Labeling %d genes that passed FDR and text thresholds.", nrow(df_labels)))
    }
  }
  
  return(list(cor_all = cor_all, cor_sig = cor_sig, plot = p))
}
