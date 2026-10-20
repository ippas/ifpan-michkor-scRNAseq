barplot_percentDetection_gene <- function(
    summary_list,
    cluster_name,
    gene_symbol,
    bar_color = "black",
    bar_width = 0.7,
    bar_size = 1,
    y_limits = NULL,
    show_p = TRUE,
    custom_title = NULL,      # manual title
    text_size_values = 12,    # pt
    text_size_p = 12,         # pt
    text_size_title = 18,     # pt
    text_size_axis = 12,      # pt default axis size
    xaxis_text_size = 12,     # NEW
    yaxis_text_size = 12      # NEW
) {
  
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  
  # ------------------------------------------------------------------
  # Internal helper: intelligent p-value formatting
  # ------------------------------------------------------------------
  format_p <- function(p) {
    if (is.na(p)) return("NA")
    if (p >= 0.0001) {
      return(formatC(p, format = "f", digits = 4))
    } else {
      return(formatC(p, format = "e", digits = 2))
    }
  }
  
  # convert pt -> mm
  text_size_values_mm <- text_size_values / 2.83
  text_size_p_mm      <- text_size_p / 2.83
  
  # ------------------------------------------------------------------
  # Extract & reshape data
  # ------------------------------------------------------------------
  df_raw <- summary_list[[cluster_name]] %>%
    dplyr::filter(gene_symbol == !!gene_symbol)
  
  df <- df_raw %>%
    dplyr::select(CTRL_percent, CORT_percent, DEX_percent) %>%
    tidyr::pivot_longer(cols = everything(),
                        names_to = "condition",
                        values_to = "percent") %>%
    dplyr::mutate(
      condition = dplyr::recode(condition,
                                "CTRL_percent" = "CTRL",
                                "CORT_percent" = "CORT",
                                "DEX_percent"  = "DEX"),
      condition = factor(condition, levels = c("CTRL","CORT","DEX"))
    )
  
  ymax <- max(df$percent, na.rm = TRUE)
  pad  <- ymax * 0.20
  
  # ------------------------------------------------------------------
  # Title selection
  # ------------------------------------------------------------------
  if (is.null(custom_title)) {
    plot_title <- paste0(gene_symbol, " detection % in ", cluster_name)
  } else {
    plot_title <- custom_title
  }
  
  # ------------------------------------------------------------------
  # Plot object
  # ------------------------------------------------------------------
  p <- ggplot(df, aes(x = condition, y = percent)) +
    geom_col(color = bar_color, fill = NA, width = bar_width, linewidth = bar_size) +
    geom_text(aes(label = paste0(round(percent,1),"%")),
              vjust = -1.3, size = text_size_values_mm) +
    scale_x_discrete(limits = c("CTRL","CORT","DEX")) +
    theme_classic(base_size = text_size_axis) +
    labs(
      title = plot_title,
      x = "",
      y = "% of cells with gene detected"
    ) +
    theme(
      plot.title = element_text(size = text_size_title, face = "plain", hjust = 0.5),
      axis.text.x = element_text(size = xaxis_text_size, face = "plain"),
      axis.text.y = element_text(size = yaxis_text_size, face = "plain"),
      axis.title  = element_text(size = text_size_axis, face = "plain")
    )
  
  if (!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits)
  }
  
  # ------------------------------------------------------------------
  # p-value brackets
  # ------------------------------------------------------------------
  if (show_p) {
    pval_cort <- suppressWarnings(as.numeric(df_raw$p_chi2_CTRL_vs_CORT))
    pval_dex  <- suppressWarnings(as.numeric(df_raw$p_chi2_CTRL_vs_DEX))
    
    top1 <- ymax + pad
    top2 <- ymax + pad * 1.8
    tick <- pad * 0.15
    
    p <- p +
      # CTRL vs CORT
      geom_segment(aes(x = 1, xend = 2, y = top1, yend = top1), linewidth = 0.8) +
      geom_segment(aes(x = 1, xend = 1, y = top1, yend = top1 - tick), linewidth = 0.8) +
      geom_segment(aes(x = 2, xend = 2, y = top1, yend = top1 - tick), linewidth = 0.8) +
      annotate("text", x = 1.5, y = top1 + pad * 0.20,
               label = paste0("p = ", format_p(pval_cort)),
               size = text_size_p_mm) +
      
      # CTRL vs DEX
      geom_segment(aes(x = 1, xend = 3, y = top2, yend = top2), linewidth = 0.8) +
      geom_segment(aes(x = 1, xend = 1, y = top2, yend = top2 - tick), linewidth = 0.8) +
      geom_segment(aes(x = 3, xend = 3, y = top2, yend = top2 - tick), linewidth = 0.8) +
      annotate("text", x = 2, y = top2 + pad * 0.20,
               label = paste0("p = ", format_p(pval_dex)),
               size = text_size_p_mm)
  }
  
  return(p)
}


barplot_percentDetection_gene <- function(
    summary_list,
    cluster_name,
    gene_symbol,
    bar_color = "black",
    bar_width = 0.7,
    bar_size = 1,
    y_limits = NULL,
    show_p = TRUE,
    custom_title = NULL,      # manual title
    text_size_values = 12,    # pt
    text_size_p = 12,         # pt
    text_size_title = 18,     # pt
    text_size_axis = 12,      # pt default axis size
    xaxis_text_size = 12,     # NEW
    yaxis_text_size = 12,     # NEW
    
    # ---------------- NEW ARGUMENTS ----------------
    stats_bracket_mult = 0.75,   # wysokość linii p-value (brackets)
    stats_text_mult    = 0.28    # wysokość tekstu p-value
) {
  
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  
  # ------------------------------------------------------------------
  # Internal helper: intelligent p-value formatting
  # ------------------------------------------------------------------
  format_p <- function(p) {
    if (is.na(p)) return("NA")
    if (p >= 0.0001) {
      return(formatC(p, format = "f", digits = 4))
    } else {
      return(formatC(p, format = "e", digits = 2))
    }
  }
  
  # convert pt -> mm for ggplot
  text_size_values_mm <- text_size_values / 2.83
  text_size_p_mm      <- text_size_p / 2.83
  
  # ------------------------------------------------------------------
  # Extract & reshape data
  # ------------------------------------------------------------------
  df_raw <- summary_list[[cluster_name]] %>%
    dplyr::filter(gene_symbol == !!gene_symbol)
  
  df <- df_raw %>%
    dplyr::select(CTRL_percent, CORT_percent, DEX_percent) %>%
    tidyr::pivot_longer(cols = everything(),
                        names_to = "condition",
                        values_to = "percent") %>%
    dplyr::mutate(
      condition = dplyr::recode(condition,
                                "CTRL_percent" = "CTRL",
                                "CORT_percent" = "CORT",
                                "DEX_percent"  = "DEX"),
      condition = factor(condition, levels = c("CTRL","CORT","DEX"))
    )
  
  ymax <- max(df$percent, na.rm = TRUE)
  pad  <- ymax * 0.20
  
  # ------------------------------------------------------------------
  # Title selection
  # ------------------------------------------------------------------
  if (is.null(custom_title)) {
    plot_title <- paste0(gene_symbol, " detection % in ", cluster_name)
  } else {
    plot_title <- custom_title
  }
  
  # ------------------------------------------------------------------
  # Plot object
  # ------------------------------------------------------------------
  p <- ggplot(df, aes(x = condition, y = percent)) +
    geom_col(color = bar_color, fill = NA, width = bar_width, linewidth = bar_size) +
    geom_text(
      aes(label = paste0(round(percent,1), "%")),
      vjust = -1.3,
      size  = text_size_values_mm
    ) +
    scale_x_discrete(limits = c("CTRL","CORT","DEX")) +
    theme_classic(base_size = text_size_axis) +
    labs(
      title = plot_title,
      x = "",
      y = "% of cells with gene detected"
    ) +
    theme(
      plot.title = element_text(size = text_size_title, face = "plain", hjust = 0.5),
      axis.text.x = element_text(size = xaxis_text_size, face = "plain"),
      axis.text.y = element_text(size = yaxis_text_size, face = "plain"),
      axis.title  = element_text(size = text_size_axis, face = "plain")
    )
  
  if (!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits)
  }
  
  # ------------------------------------------------------------------
  # p-value brackets with dynamic scaling
  # ------------------------------------------------------------------
  if (show_p) {
    pval_cort <- suppressWarnings(as.numeric(df_raw$p_chi2_CTRL_vs_CORT))
    pval_dex  <- suppressWarnings(as.numeric(df_raw$p_chi2_CTRL_vs_DEX))
    
    top1 <- ymax + pad * stats_bracket_mult
    top2 <- ymax + pad * stats_bracket_mult * 1.8
    tick <- pad * 0.15
    
    p <- p +
      # CTRL vs CORT
      geom_segment(aes(x = 1, xend = 2, y = top1, yend = top1), linewidth = 0.8) +
      geom_segment(aes(x = 1, xend = 1, y = top1, yend = top1 - tick), linewidth = 0.8) +
      geom_segment(aes(x = 2, xend = 2, y = top1, yend = top1 - tick), linewidth = 0.8) +
      annotate("text",
               x = 1.5,
               y = top1 + pad * stats_text_mult,
               label = paste0("p = ", format_p(pval_cort)),
               size = text_size_p_mm) +
      
      # CTRL vs DEX
      geom_segment(aes(x = 1, xend = 3, y = top2, yend = top2), linewidth = 0.8) +
      geom_segment(aes(x = 1, xend = 1, y = top2, yend = top2 - tick), linewidth = 0.8) +
      geom_segment(aes(x = 3, xend = 3, y = top2, yend = top2 - tick), linewidth = 0.8) +
      annotate("text",
               x = 2,
               y = top2 + pad * stats_text_mult,
               label = paste0("p = ", format_p(pval_dex)),
               size = text_size_p_mm)
  }
  
  return(p)
}

