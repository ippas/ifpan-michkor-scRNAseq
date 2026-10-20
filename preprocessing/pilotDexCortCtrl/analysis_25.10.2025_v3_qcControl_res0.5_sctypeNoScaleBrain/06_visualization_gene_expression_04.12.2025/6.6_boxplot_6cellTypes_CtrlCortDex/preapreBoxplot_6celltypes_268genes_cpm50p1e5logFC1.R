# ============================================================
# PĘTLA: boxploty log2FC dla cluster_1 ... cluster_5
# ============================================================

clusters <- paste0("cluster_", 1:5)

for (cl in clusters) {
  
  genes_vec <- res_heatmap$genes_per_cluster[[cl]]
  n_genes   <- length(genes_vec)
  
  message("📊 Generating boxplot for ", cl,
          " (n genes = ", n_genes, ")")
  
  boxplot_log2FC_geneSet(
    combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
    genes_vec,
    
    # === TEXT SIZES ===
    text_size_mean         = 0,
    text_size_all          = 0,
    text_size_thresh       = 0,
    text_size_axis         = 12,
    text_size_axis_title   = 12,
    text_size_legend       = 12,
    text_size_legend_title = 12,
    text_size_title        = 14,
    text_size_strip        = 12,
    text_size_base         = 12,
    # y_limits =              c(-5, 5),
    
    # === AXIS TITLES ===
    x_axis_title = "Neural cell types",
    y_axis_title = "log2FC (CORT / DEX vs CTRL)",
    
    # === TITLE WITH GENE COUNT ===
    plot_title = paste0(
      "GR-regulated genes — ", cl,
      " (n = ", n_genes, ")"
    ),
    
    # === SVG EXPORT ===
    output_svg = file.path(
      "results/figures/boxplots_log2FC_celltypes_12.12.2025/boxplot_groupByCelltype_15.12.2025",
      paste0(cl, "wthoutDots.svg")
    ),
    svg_width  = 5,
    svg_height = 7,
    
    # === LEGEND ===
    legend_position = "top",
    
    # === DIR CONTROL ===
    force_create_outputdir = TRUE,
    
    add_dots_jitter = F
  )
}




boxplot_log2FC_geneSet_groupByTreatment(
  combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
  res_heatmap$genes_per_cluster$cluster_4,
  output_svg = "results/figures/boxplots_log2FC_celltypes_12.12.2025/test2/test.svg",
  force_create_outputdir = TRUE
)


# ============================================================
# PĘTLA: boxploty log2FC grouped by treatment (cluster_1 ... cluster_5)
# ============================================================

clusters <- paste0("cluster_", 1:5)

for (cl in clusters) {
  
  genes_vec <- res_heatmap$genes_per_cluster[[cl]]
  n_genes   <- length(genes_vec)
  
  message(
    "📊 Generating TREATMENT-grouped boxplot for ", cl,
    " (n genes = ", n_genes, ")"
  )
  
  boxplot_log2FC_geneSet_groupByTreatment(
    combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
    genes_vec,
    
    # === TEXT SIZES ===
    text_size_mean         = 10,
    text_size_all          = 8,
    text_size_thresh       = 8,
    text_size_axis         = 12,
    text_size_axis_title   = 12,
    text_size_legend       = 12,
    text_size_legend_title = 12,
    text_size_title        = 14,
    text_size_strip        = 12,
    text_size_base         = 12,
    
    # === Y RANGE ===
    y_limits = c(-5, 5),
    
    # === AXIS TITLES ===
    x_axis_title = "Treatment | Neural cell type",
    y_axis_title = "log2FC (CORT / DEX vs CTRL)",
    
    # === TITLE WITH GENE COUNT ===
    plot_title = paste0(
      "GR-regulated genes — ", cl,
      " (n = ", n_genes, ")"
    ),
    
    # === SVG EXPORT ===
    output_svg = file.path(
      "results/figures/boxplots_log2FC_celltypes_12.12.2025/boxplot_groupByTreatment/",
      paste0(cl, ".svg")
    ),
    svg_width  = 10,
    svg_height = 7,
    
    # === LEGEND ===
    legend_position = "top",
    
    # === DIR CONTROL ===
    force_create_outputdir = TRUE
  )
}
