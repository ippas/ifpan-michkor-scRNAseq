barplot_enrichr_results <- function(
    df,
    min_genes = 2,
    n_top = 10,
    title = "Enrichr results",
    axis_text_size = 12,
    axis_title_size = 14,
    legend_text_size = 10,
    legend_title_size = 12,
    plot_title_size = 16,
    label_text_size = 3.5,
    fill_limits = NULL,       # np. c(0, 10) lub NULL
    squish_to_limits = FALSE, # NEW: TRUE => wartości poza zakresem są „dociśnięte” do granic
    symmetric_limits = FALSE  # opcjonalnie: TRUE => automatyczne symetryczne limity ±max(|log2_or|)
) {
  df2 <- df %>%
    dplyr::filter(n_genes >= min_genes) %>%
    dplyr::arrange(fdr, dplyr::desc(odds_ratio)) %>%
    dplyr::slice_head(n = n_top) %>%
    dplyr::mutate(
      term = factor(term, levels = rev(term)),
      log2_or = log2(odds_ratio)
    )
  
  # Jeśli chcesz automatycznie symetryczne limity (przydatne dla gradient2 z midpoint=0)
  if (isTRUE(symmetric_limits) && is.null(fill_limits)) {
    m <- max(abs(df2$log2_or), na.rm = TRUE)
    fill_limits <- c(-m, m)
  }
  
  oob_fun <- if (isTRUE(squish_to_limits)) scales::squish else scales::censor
  
  ggplot(df2, aes(x = term, y = -log10(fdr), fill = log2_or)) +
    geom_col() +
    geom_text(aes(label = overlap), hjust = -0.2, size = label_text_size, color = "black") +
    coord_flip() +
    scale_fill_gradient2(
      low = scales::muted("blue"),
      mid = "white",
      high = scales::muted("red"),
      midpoint = 0,
      limits = fill_limits,
      oob = oob_fun              # <-- klucz: „squish” zamiast domyślnego „censor”
    ) +
    labs(
      title = title,
      x = NULL,
      y = expression(-log[10](FDR)),
      fill = expression(log[2]~"(Odds Ratio)")
    ) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = legend_title_size, hjust = 0.5, color = "black"),
      legend.text = element_text(size = legend_text_size, color = "black"),
      axis.text = element_text(size = axis_text_size, color = "black"),
      axis.title = element_text(size = axis_title_size, color = "black"),
      plot.title = element_text(size = plot_title_size, hjust = 0.5, color = "black")
    ) +
    guides(fill = guide_colorbar(title.position = "left", title.hjust = 0.5)) +
    expand_limits(y = max(-log10(df$fdr), na.rm = TRUE) * 1.1)
}
