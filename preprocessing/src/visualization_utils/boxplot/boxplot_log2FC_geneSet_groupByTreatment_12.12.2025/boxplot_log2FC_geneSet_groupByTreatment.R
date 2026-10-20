boxplot_log2FC_geneSet_groupByTreatment <- function(
    data_list,
    genes_vector,
    
    # === CELL TYPES ===
    celltypes = c(
      "Radial glial cells",
      "Astrocytes",
      "OPCs",
      "Oligodendrocytes",
      "Immature neurons",
      "Mature neurons"
    ),
    
    # === TEXT SIZES ===
    text_size_mean         = 12,
    text_size_all          = 12,
    text_size_thresh       = 12,
    text_size_axis         = 12,
    text_size_axis_title   = 12,
    text_size_legend       = 12,
    text_size_legend_title = 12,
    text_size_title        = 14,
    text_size_strip        = 12,
    text_size_base         = 12,
    
    # === Y POSITION CONTROL ===
    mean_y_offset   = 0.10,
    y_all_offset    = 1.45,
    y_thresh_offset = 1.85,
    
    # === DISPLAY OPTIONS ===
    plot_title      = NULL,
    x_axis_title    = "Treatment | Cell type",
    y_axis_title    = "log2 fold change vs CTRL",
    y_limits        = NULL,
    color_map       = c("CORT" = "firebrick2", "DEX" = "steelblue3"),
    alpha_jitter    = 0.35,
    alpha_box       = 1,
    
    # === LEGEND & THEME ===
    legend_position = "right",
    theme_name      = "classic",
    
    # === OUTPUT ===
    output_svg             = NULL,
    svg_width              = 12,
    svg_height             = 8,
    force_create_outputdir = FALSE
) {
  
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  
  # === SVG OUTPUT HANDLING (UNCHANGED) ===
  if (!is.null(output_svg)) {
    out_dir <- dirname(output_svg)
    if (!dir.exists(out_dir)) {
      if (force_create_outputdir) {
        dir.create(out_dir, recursive = TRUE)
        message("📁 Created output directory: ", out_dir)
      } else {
        stop(
          "Output directory does not exist: ", out_dir,
          "\nSet force_create_outputdir = TRUE to create it."
        )
      }
    }
    svg(filename = output_svg, width = svg_width, height = svg_height)
    on.exit(dev.off(), add = TRUE)
  }
  
  # === text scaling ===
  text_size_mean_scaled   <- text_size_mean   / 2.83
  text_size_all_scaled    <- text_size_all    / 2.83
  text_size_thresh_scaled <- text_size_thresh / 2.83
  
  # ============================================================================
  # DATA PREP (SAME)
  # ============================================================================
  
  df_plot_fc <- data_list[celltypes] %>%
    map(function(x) {
      x %>%
        filter(gene_symbol %in% genes_vector) %>%
        select(gene_symbol, log2FC_CORT_CPM, log2FC_DEX_CPM) %>%
        pivot_longer(
          cols      = starts_with("log2FC_"),
          names_to  = "treatment",
          values_to = "log2FC"
        ) %>%
        mutate(
          treatment = case_when(
            treatment == "log2FC_CORT_CPM" ~ "CORT",
            treatment == "log2FC_DEX_CPM"  ~ "DEX"
          )
        )
    }) %>%
    bind_rows(.id = "celltype") %>%
    mutate(
      treatment = factor(treatment, levels = c("CORT", "DEX")),
      celltype  = factor(celltype, levels = celltypes)
    ) %>%
    # <<< CHANGED >>>
    mutate(
      x_group = factor(
        paste(treatment, celltype, sep = " | "),
        levels = c(
          paste("CORT", celltypes, sep = " | "),
          paste("DEX",  celltypes, sep = " | ")
        )
      )
    )
  
  # ============================================================================
  # Y RANGE (SAME)
  # ============================================================================
  
  if (is.null(y_limits)) {
    y_min <- min(df_plot_fc$log2FC, na.rm = TRUE)
    y_max <- max(df_plot_fc$log2FC, na.rm = TRUE)
  } else {
    y_min <- y_limits[1]
    y_max <- y_limits[2]
  }
  
  y_mean <- y_max - mean_y_offset * (y_max - y_min)
  
  n_genes_total <- n_distinct(df_plot_fc$gene_symbol)
  
  # ============================================================================
  # COUNTS / LABELS (SAME LOGIC, DIFFERENT GROUPING)
  # ============================================================================
  
  df_counts_thresh <- df_plot_fc %>%
    group_by(x_group, treatment) %>%
    nest() %>%
    mutate(
      n_up   = map_int(data, ~ sum(.x$log2FC >=  1, na.rm = TRUE)),
      n_down = map_int(data, ~ sum(.x$log2FC <= -1, na.rm = TRUE))
    ) %>%
    ungroup()
  
  df_labels_thresh <- bind_rows(
    df_counts_thresh %>% mutate(y =  y_thresh_offset,  label = paste0("n = ", n_up)),
    df_counts_thresh %>% mutate(y = -y_thresh_offset, label = paste0("n = ", n_down))
  )
  
  df_counts_all <- df_plot_fc %>%
    group_by(x_group, treatment) %>%
    nest() %>%
    mutate(
      n_all_up   = map_int(data, ~ sum(.x$log2FC >  0, na.rm = TRUE)),
      n_all_down = map_int(data, ~ sum(.x$log2FC <  0, na.rm = TRUE))
    ) %>%
    ungroup()
  
  df_labels_all <- bind_rows(
    df_counts_all %>% mutate(y =  y_all_offset,  label = paste0("all↑ = ", n_all_up)),
    df_counts_all %>% mutate(y = -y_all_offset, label = paste0("all↓ = ", n_all_down))
  )
  
  df_stats_mean <- df_plot_fc %>%
    group_by(x_group, treatment) %>%
    nest() %>%
    mutate(mean_log2FC = map_dbl(data, ~ mean(.x$log2FC, na.rm = TRUE))) %>%
    ungroup() %>%
    mutate(
      y = y_mean,
      label = paste0("x̄ = ", round(mean_log2FC, 2))
    )
  
  if (is.null(plot_title)) {
    plot_title <- paste0(
      "Distribution of log2FC values grouped by treatment (n = ",
      n_genes_total,
      " genes)"
    )
  }
  
  base_theme <- switch(
    theme_name,
    classic = theme_classic(base_size = text_size_base),
    bw      = theme_bw(base_size = text_size_base),
    minimal = theme_minimal(base_size = text_size_base),
    light   = theme_light(base_size = text_size_base),
    dark    = theme_dark(base_size = text_size_base),
    gray    = theme_gray(base_size = text_size_base),
    void    = theme_void(base_size = text_size_base),
    test    = theme_test(base_size = text_size_base),
    theme_classic(base_size = text_size_base)
  )
  
  ggplot(df_plot_fc, aes(x = x_group, y = log2FC, color = treatment)) +
    
    geom_boxplot(
      width = 0.6,
      linewidth = 0.6,
      alpha = alpha_box,
      outlier.shape = NA
    ) +
    
    geom_jitter(
      width = 0.15,
      alpha = alpha_jitter,
      size = 0.8
    ) +
    
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.5) +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed",
               color = "grey40", linewidth = 0.6) +
    
    geom_text(
      data = df_labels_all,
      aes(x = x_group, y = y, label = label),
      size = text_size_all_scaled,
      color = "black"
    ) +
    
    geom_text(
      data = df_labels_thresh,
      aes(x = x_group, y = y, label = label),
      size = text_size_thresh_scaled,
      color = "black"
    ) +
    
    geom_text(
      data = df_stats_mean,
      aes(x = x_group, y = y, label = label),
      size = text_size_mean_scaled,
      color = "black"
    ) +
    
    scale_color_manual(values = color_map) +
    scale_y_continuous(limits = c(y_min - 0.4, y_max + 1.0)) +
    
    base_theme +
    theme(
      legend.position = legend_position,
      axis.text.x      = element_text(angle = 45, hjust = 1, size = text_size_axis),
      axis.text.y      = element_text(size = text_size_axis),
      axis.title       = element_text(size = text_size_axis_title),
      legend.text      = element_text(size = text_size_legend),
      legend.title     = element_text(size = text_size_legend_title),
      plot.title       = element_text(size = text_size_title,
                                      hjust = 0.5, face = "bold"),
      strip.text       = element_text(size = text_size_strip)
    ) +
    
    labs(
      title = plot_title,
      x = x_axis_title,
      y = y_axis_title,
      color = "Treatment"
    )
}

# 
# boxplot_log2FC_geneSet_groupByTreatment(
#   combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
#   res_heatmap$genes_per_cluster$cluster_4,
#   # output_svg = "figures/log2FC_cluster1_byTreatment.svg",
#   # force_create_outputdir = TRUE
# )
