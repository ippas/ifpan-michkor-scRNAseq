boxplot_log2FC_geneSet <- function(
    data_list,
    genes_vector,
    
    # === TEXT SIZES ===
    text_size_all            = 4,
    text_size_thresh         = 4,
    text_size_mean           = 4,
    text_size_axis           = 11,
    text_size_axis_title     = 12,
    text_size_legend         = 11,
    text_size_legend_title   = 12,
    text_size_title          = 14,
    text_size_strip          = 11,
    text_size_base           = 11
) {
  
  # ----------------------------------------------------------------------------
  # Visualizes log2FC distributions for a predefined gene set (genes_vector).
  # Not genome-wide.
  # ----------------------------------------------------------------------------
  
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  
  # ============================================================================
  # 1. DATA PREP
  # ============================================================================
  
  df_plot_fc <- data_list[c(
    "Radial glial cells",
    "Astrocytes",
    "OPCs",
    "Oligodendrocytes",
    "Immature neurons",
    "Mature neurons"
  )] %>%
    purrr::map(function(x) {
      x %>%
        filter(gene_symbol %in% genes_vector) %>%
        select(
          gene_symbol,
          log2FC_CORT_CPM,
          log2FC_DEX_CPM
        ) %>%
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
      celltype  = factor(
        celltype,
        levels = c(
          "Radial glial cells",
          "Astrocytes",
          "OPCs",
          "Oligodendrocytes",
          "Immature neurons",
          "Mature neurons"
        )
      )
    )
  
  # ============================================================================
  # 2. GLOBAL PARAMS
  # ============================================================================
  
  y_max <- max(df_plot_fc$log2FC, na.rm = TRUE)
  y_min <- min(df_plot_fc$log2FC, na.rm = TRUE)
  
  y_mean <- y_max - 0.10 * (y_max - y_min)
  
  n_genes_total <- df_plot_fc %>%
    distinct(gene_symbol) %>%
    nrow()
  
  # ============================================================================
  # 3. THRESHOLD COUNTS (|log2FC| >= 1)
  # ============================================================================
  
  df_counts_thresh <- df_plot_fc %>%
    group_by(celltype, treatment) %>%
    nest() %>%
    mutate(
      n_up   = map_int(data, ~ sum(.x$log2FC >=  1, na.rm = TRUE)),
      n_down = map_int(data, ~ sum(.x$log2FC <= -1, na.rm = TRUE))
    ) %>%
    select(celltype, treatment, n_up, n_down) %>%
    ungroup()
  
  df_labels_thresh <- bind_rows(
    df_counts_thresh %>% mutate(y =  1.85, label = paste0("n = ", n_up)),
    df_counts_thresh %>% mutate(y = -1.85, label = paste0("n = ", n_down))
  )
  
  # ============================================================================
  # 4. ALL UP / DOWN
  # ============================================================================
  
  df_counts_all <- df_plot_fc %>%
    group_by(celltype, treatment) %>%
    nest() %>%
    mutate(
      n_all_up   = map_int(data, ~ sum(.x$log2FC >  0, na.rm = TRUE)),
      n_all_down = map_int(data, ~ sum(.x$log2FC <  0, na.rm = TRUE))
    ) %>%
    select(celltype, treatment, n_all_up, n_all_down) %>%
    ungroup()
  
  df_labels_all <- bind_rows(
    df_counts_all %>% mutate(y =  1.45, label = paste0("all↑ = ", n_all_up)),
    df_counts_all %>% mutate(y = -1.45, label = paste0("all↓ = ", n_all_down))
  )
  
  # ============================================================================
  # 5. MEAN
  # ============================================================================
  
  df_stats_mean <- df_plot_fc %>%
    group_by(celltype, treatment) %>%
    nest() %>%
    mutate(
      mean_log2FC = map_dbl(data, ~ mean(.x$log2FC, na.rm = TRUE))
    ) %>%
    select(celltype, treatment, mean_log2FC) %>%
    ungroup() %>%
    mutate(
      y = y_mean,
      label = paste0("x̄ = ", round(mean_log2FC, 2))
    )
  
  # ============================================================================
  # 6. PLOT
  # ============================================================================
  
  ggplot(df_plot_fc, aes(x = celltype, y = log2FC, color = treatment)) +
    
    geom_boxplot(
      aes(group = interaction(celltype, treatment)),
      position = position_dodge(width = 0.75),
      width = 0.6,
      linewidth = 0.6,
      outlier.shape = NA
    ) +
    
    geom_jitter(
      position = position_jitterdodge(0.15, 0.75),
      alpha = 0.35,
      size = 0.8
    ) +
    
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.5) +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
    
    geom_text(
      data = df_labels_all,
      aes(x = celltype, y = y, label = label,
          group = interaction(celltype, treatment)),
      position = position_dodge(0.75),
      size = text_size_all,
      color = "black",
      vjust = ifelse(df_labels_all$y > 0, 0, 1)
    ) +
    
    geom_text(
      data = df_labels_thresh,
      aes(x = celltype, y = y, label = label,
          group = interaction(celltype, treatment)),
      position = position_dodge(0.75),
      size = text_size_thresh,
      color = "black",
      vjust = ifelse(df_labels_thresh$y > 0, 0, 1)
    ) +
    
    geom_text(
      data = df_stats_mean,
      aes(x = celltype, y = y, label = label,
          group = interaction(celltype, treatment)),
      position = position_dodge(0.75),
      size = text_size_mean,
      color = "black",
      vjust = 0
    ) +
    
    scale_color_manual(
      values = c("CORT" = "firebrick2", "DEX" = "steelblue3")
    ) +
    
    scale_y_continuous(limits = c(y_min - 0.4, y_max + 1.0)) +
    
    theme_classic(base_size = text_size_base) +
    theme(
      axis.text        = element_text(size = text_size_axis),
      axis.title       = element_text(size = text_size_axis_title),
      legend.text      = element_text(size = text_size_legend),
      legend.title     = element_text(size = text_size_legend_title),
      plot.title       = element_text(size = text_size_title, hjust = 0.5, face = "bold"),
      strip.text       = element_text(size = text_size_strip)
    ) +
    
    labs(
      title = paste0(
        "Distribution of log2FC values for selected gene set (n = ",
        n_genes_total,
        " genes)"
      ),
      x = "Cell type",
      y = "log2 fold change vs CTRL",
      color = "Treatment"
    )
}


# boxplot_log2FC_geneSet(
#   data_list   = combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
#   genes_vector = res_heatmap$genes_per_cluster$cluster_5,
#   text_size_mean = 4,
#   text_size_all  = 4,
#   text_size_title = 16
# )
