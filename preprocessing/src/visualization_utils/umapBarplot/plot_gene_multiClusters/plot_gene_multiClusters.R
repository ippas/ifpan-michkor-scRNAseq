plot_gene_multiClusters <- function(
    summary_list,
    cluster_names,
    gene_symbol,
    gene_title_size = 22,
    log2_transform = TRUE,
    border_size = 1,      # grubość ramki
    n_rows = 1            # <---- liczba wierszy (np. 1,2,3…)
) {
  require(ggplot2)
  require(patchwork)
  require(purrr)
  require(cowplot)
  require(grid)
  
  add_border <- function(plot_obj, size = 1, color = "black") {
    cowplot::ggdraw(plot_obj) +
      cowplot::draw_plot(plot_obj) +
      cowplot::draw_grob(
        grid::rectGrob(
          gp = grid::gpar(col = color, fill = NA, lwd = size)
        )
      )
  }
  
  make_pair <- function(cl_name) {
    
    p_cpm <- barplot_meanCPM_gene(
      summary_list    = summary_list,
      cluster_name    = cl_name,
      gene_symbol     = gene_symbol,
      custom_title    = cl_name,
      bar_color = "black",
      point_size = 5,
      bar_width = 0.7,
      bar_size = 1,
      log2_transform = log2_transform,
      add_points = TRUE,
      se_width = 0.1,
      show_p = TRUE,
      text_size_values = 0,
      text_size_p = 12,
      text_size_title = 12,
      text_size_axis = 12,
      stats_bracket_mult = 0.75,
      stats_text_mult    = 0.35,
      xaxis_text_size = 12,
      yaxis_text_size = 12
    )
    
    p_pct <- barplot_percentDetection_gene(
      summary_list    = summary_list,
      cluster_name    = cl_name,
      gene_symbol     = gene_symbol,
      custom_title    = cl_name,
      show_p          = TRUE,
      bar_color = "black",
      bar_width = 0.7,
      text_size_values = 10,
      text_size_p = 12,
      text_size_title = 12,
      text_size_axis = 12,
      yaxis_text_size = 12,
      xaxis_text_size = 12,
      stats_bracket_mult = 0.75,
      stats_text_mult    = 0.35
    )
    
    pair_plot <- (p_cpm | p_pct)
    
    return(add_border(pair_plot, size = border_size, color = "black"))
  }
  
  p3 <- purrr::map(cluster_names, make_pair) %>%
    wrap_plots(nrow = n_rows) +     # <---- użycie argumentu n_rows
    plot_annotation(
      title = gene_symbol,
      theme = theme(
        plot.title = element_text(size = gene_title_size,
                                  face = "bold",
                                  hjust = 0.5)
      )
    )
  
  
  seurat_featurePlot_binary_clusterSummary(
    seurat_obj = integrated_pilotDexCortCtrl_filtered,
    gene_pattern = gene_symbol,
    split.by = "treatment",
    split_levels = c("CTRL", "CORT", "DEX"),
    cluster_col = "res0.3_RNA_noScaled_brain",
    label_size = 12,
    point_size = 0.4
  ) -> p1
  
  
  seurat_featurePlot_continuous_v3(
    seurat_obj = integrated_pilotDexCortCtrl_filtered,
    gene_pattern = gene_symbol,
    split.by = "treatment",
    split_levels = c("CTRL", "CORT", "DEX"),
    cluster_col = "res0.3_RNA_noScaled_brain",
    point_size = 0.4,
    winsorize_top = 0.01,
    scale_zeroOne = TRUE
  ) -> p2
  
  ((p1 + p2) / p3) +
    plot_layout(heights = c(1, 4)) -> final_plot
  
  return(final_plot)
}

# example
# svg("results/results_figures_26.11.2025_umapBarplotContinousBinary/umpaBarplots_ZHX3.svg",
#     width = 16,
#     height = 20)
# plot_gene_multiClusters(
#   summary_list  = combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
#   cluster_names = c("Radial glial cells", "Astrocytes",
#                     "Oligodendrocyte precursor cells", "Oligodendrocytes", 
#                     "Immature neurons", "Mature neurons"),
#   gene_symbol   = "ZHX3",
#   gene_title_size = 14,
#   n_rows = 3,
#   log2_transform = FALSE,
#   border_size = 2
# )
# 
# dev.off()
