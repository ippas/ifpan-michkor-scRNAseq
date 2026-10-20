generate_heatmap_test_clustering <- function(
    data_list,
    celltypes,
    distance_method = c("spearman", "pearson", "euclidean", "manhattan"),
    linkage_method = c("ward.D2", "average", "complete", "single"),
    show_rows = TRUE,
    show_cols = FALSE
) {
  require(dplyr)
  require(tidyr)
  require(tibble)
  require(purrr)
  require(ComplexHeatmap)
  require(circlize)
  require(grid)
  require(RColorBrewer)
  
  message("=== 1) Wybór genów ===")
  genes2heatmap <- data_list[celltypes] %>%
    lapply(function(x){
      x %>%
        filter(meanCPM_CORT > 50 | meanCPM_DEX > 50) %>%
        filter(log2FC_edgerCORT > 1 | log2FC_edgerDEX > 1) %>%
        filter(PValue_edgerCORT < 0.01 | PValue_edgerDEX < 0.01) %>%
        .$gene_symbol
    }) %>% unlist() %>% unique()
  
  message("Liczba genów spełniających kryteria: ", length(genes2heatmap))
  
  message("=== 2) Normalizacja + Z-score ===")
  df_list <- lapply(celltypes, function(cl) {
    
    df_raw <- data_list[[cl]] %>%
      filter(gene_symbol %in% genes2heatmap) %>%
      select(gene_symbol, paste0("QN_A", 1:12))
    
    ctrl_mean <- df_raw %>% select(QN_A1:QN_A4) %>% rowMeans()
    df_norm <- df_raw
    df_norm[, -1] <- sweep(df_norm[, -1], 1, ctrl_mean, "/")
    
    ctrl_sd <- df_norm %>% select(QN_A1:QN_A4) %>% apply(1, sd)
    
    zscores <- as.matrix(df_norm[, -1])
    zscores <- sweep(zscores, 1, 1, "-")
    zscores <- sweep(zscores, 1, ctrl_sd, "/")
    zscores[!is.finite(zscores)] <- 0
    
    df_z <- data.frame(
      gene_symbol = df_norm$gene_symbol,
      as.data.frame(zscores)
    )
    
    colnames(df_z)[-1] <- paste0(gsub(" ", "_", cl), "_", colnames(df_z)[-1])
    df_z
  })
  
  names(df_list) <- celltypes
  
  message("=== 3) Full merge ===")
  merged <- purrr::reduce(df_list, full_join, by = "gene_symbol") %>%
    column_to_rownames("gene_symbol")
  
  mat <- as.matrix(merged)
  mat[!is.finite(mat)] <- 0
  
  message("=== 4) Kolejność CTRL → CORT → DEX ===")
  
  ctrl_cols <- grep("A([1-4])$", colnames(mat), value = TRUE)
  cort_cols <- grep("A(5|6|7|8)$", colnames(mat), value = TRUE)
  dex_cols  <- grep("A(9|1[0-2])$", colnames(mat), value = TRUE)
  
  mat <- mat[, c(ctrl_cols, cort_cols, dex_cols)]
  
  groups <- ifelse(colnames(mat) %in% ctrl_cols, "CTRL",
                   ifelse(colnames(mat) %in% cort_cols, "CORT", "DEX"))
  
  celltype_vec <- sapply(colnames(mat), function(x){
    strsplit(x, "_A")[[1]][1]
  })
  
  celltype_colors <- setNames(
    RColorBrewer::brewer.pal(length(unique(celltype_vec)), "Dark2"),
    unique(celltype_vec)
  )
  
  message("=== 5) Annotation ===")
  
  ha <- HeatmapAnnotation(
    df = data.frame(
      group = groups,
      celltype = celltype_vec
    ),
    col = list(
      group    = c("CTRL"="grey40", "CORT"="firebrick2", "DEX"="steelblue3"),
      celltype = celltype_colors
    ),
    annotation_name_gp = gpar(fontsize = 12, fontface = "bold"),
    annotation_legend_param = list(
      group    = list(title = "group", ncol = 1),
      celltype = list(title = "celltype", ncol = 1)
    )
  )
  
  message("=== 6) Draw heatmap ===")
  
  distance_method <- match.arg(distance_method)
  linkage_method  <- match.arg(linkage_method)
  
  dist_fun <- if(distance_method %in% c("pearson","spearman")){
    function(x) as.dist(1 - cor(t(x), method = distance_method))
  } else { distance_method }
  
  ht <- Heatmap(
    mat,
    name = "Z-score vs CTRL",
    top_annotation = ha,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    clustering_distance_rows = dist_fun,
    clustering_method_rows = linkage_method,
    show_row_names = show_rows,
    show_column_names = FALSE,
    col = colorRamp2(c(-10, 0, 10), c("navy", "white", "firebrick4")),
    row_title = paste0("Genes (n=", nrow(mat), ")"),
    row_title_gp = gpar(fontsize = 16, fontface = "bold"),
    column_title = paste("Clustering:", distance_method, "+", linkage_method),
    column_title_gp = gpar(fontsize = 18, fontface = "bold"),
    heatmap_legend_param = list(
      title = "Z-score vs CTRL",
      ncol = 1,
      title_gp = gpar(fontsize = 14, fontface = "bold")
    )
  )
  
  draw(ht, merge_legend = TRUE)
}


celltypes <- c(
  "Astrocytes",
  "Radial glial cells",
  "Oligodendrocyte precursor cells",
  "Oligodendrocytes",
  "Immature neurons",
  "Mature neurons"
)
generate_heatmap_test_clustering(
  combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
  celltypes,
  distance_method = "spearman",
  linkage_method = "ward.D2"
)



# svg("/home/mateusz/projects/ifpan-michkor-scRNAseq/results/results_heatmap_28.11.2025/heatmap_brain_GR.svg",
#     width = 12, height = 40)
# 
# ht <- generate_heatmap_test_clustering(
#   combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
#   celltypes,
#   distance_method = "spearman",
#   linkage_method = "ward.D2"
# )
# 
# draw(ht, merge_legend = TRUE)
# dev.off()



distance_methods <- c("spearman", "pearson", "euclidean", "manhattan")
linkage_methods  <- c("ward.D2", "average", "complete", "single")

output_dir <- "/home/mateusz/projects/ifpan-michkor-scRNAseq/results/results_heatmap_30.11.2025"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

for (d in distance_methods) {
  for (l in linkage_methods) {
    
    message("\n============================================")
    message(" Generating heatmap for: ", d, " + ", l)
    message("============================================")
    
    svg(
      filename = file.path(
        output_dir,
        paste0("heatmap_", d, "_", l, ".svg")
      ),
      width = 12,
      height = 8
    )
    
    ht <- generate_heatmap_test_clustering(
      combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
      celltypes,
      distance_method = d,
      linkage_method = l
    )
    
    draw(ht, merge_legend = TRUE)
    dev.off()
  }
}

