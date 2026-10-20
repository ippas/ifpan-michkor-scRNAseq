generate_heatmap_test_clustering <- function(
    data_list,
    celltypes,
    distance_method = c("spearman", "pearson", "euclidean", "manhattan"),
    linkage_method  = c("ward.D2", "average", "complete", "single"),
    cpm_threshold   = 50,
    pval_threshold  = 0.01,
    logfc_threshold = 1,
    logfc_mode      = c("up", "down", "both"),
    contrast        = c("CORT", "DEX", "both"),
    zscore_range    = c(-4, 0, 4),
    show_rows       = TRUE,
    show_cols       = FALSE,
    row_dend_width  = unit(40, "mm"),
    cut_k           = NULL,
    cut_height      = NULL,
    output_path     = NULL
) {
  require(dplyr); require(tidyr); require(tibble); require(purrr)
  require(ComplexHeatmap); require(circlize); require(grid)
  require(RColorBrewer); require(dendextend)
  
  logfc_mode      <- match.arg(logfc_mode)
  contrast        <- match.arg(contrast)
  distance_method <- match.arg(distance_method)
  linkage_method  <- match.arg(linkage_method)
  
  # =============================
  # 1) Wybór genów
  # =============================
  message("=== 1) Wybór genów ===")
  
  genes2heatmap <- data_list[celltypes] %>%
    lapply(function(x){
      df <- x %>% filter(meanCPM_CORT > cpm_threshold | meanCPM_DEX > cpm_threshold)
      
      if (logfc_mode == "up")
        df <- df %>% filter(log2FC_edgerCORT > logfc_threshold | log2FC_edgerDEX > logfc_threshold)
      else if (logfc_mode == "down")
        df <- df %>% filter(log2FC_edgerCORT < -logfc_threshold | log2FC_edgerDEX < -logfc_threshold)
      else
        df <- df %>% filter(abs(log2FC_edgerCORT) > logfc_threshold | abs(log2FC_edgerDEX) > logfc_threshold)
      
      if (contrast == "CORT")
        df <- df %>% filter(PValue_edgerCORT < pval_threshold)
      else if (contrast == "DEX")
        df <- df %>% filter(PValue_edgerDEX < pval_threshold)
      else
        df <- df %>% filter(PValue_edgerCORT < pval_threshold | PValue_edgerDEX < pval_threshold)
      
      df$gene_symbol
    }) %>% unlist() %>% unique()
  
  message("Liczba genów spełniających kryteria: ", length(genes2heatmap))
  
  # =============================
  # 2) MERGE surowych danych
  # =============================
  message("=== 2) Full merge raw ===")
  
  df_list <- lapply(celltypes, function(cl) {
    data_list[[cl]] %>%
      filter(gene_symbol %in% genes2heatmap) %>%
      select(gene_symbol, paste0("QN_A", 1:12)) %>%
      rename_with(~ paste0(gsub(" ", "_", cl), "_", .x), -gene_symbol)
  })
  
  merged <- purrr::reduce(df_list, full_join, by = "gene_symbol") %>%
    column_to_rownames("gene_symbol")
  
  mat_raw <- as.matrix(merged)
  mat_raw[!is.finite(mat_raw)] <- 0
  
  # =============================
  # 3) Globalna normalizacja vs CTRL
  # =============================
  message("=== 3) GLOBAL Z-score vs CTRL (ALL CELLS) ===")
  
  ctrl_cols <- grep("A([1-4])$", colnames(mat_raw), value = TRUE)
  
  global_ctrl_mean <- rowMeans(mat_raw[, ctrl_cols], na.rm = TRUE)
  global_ctrl_sd   <- apply(mat_raw[, ctrl_cols], 1, sd, na.rm = TRUE)
  
  mat_norm <- sweep(mat_raw, 1, global_ctrl_mean, "-")
  mat_norm <- sweep(mat_norm, 1, global_ctrl_sd, "/")
  mat_norm[!is.finite(mat_norm)] <- 0
  
  # =============================
  # 4) Sortowanie kolumn CTRL->CORT->DEX
  # =============================
  cort_cols <- grep("A(5|6|7|8)$", colnames(mat_norm), value = TRUE)
  dex_cols  <- grep("A(9|1[0-2])$", colnames(mat_norm), value = TRUE)
  
  mat <- mat_norm[, c(ctrl_cols, cort_cols, dex_cols)]
  
  groups <- ifelse(colnames(mat) %in% ctrl_cols, "CTRL",
                   ifelse(colnames(mat) %in% cort_cols, "CORT", "DEX"))
  
  celltype_vec <- sapply(colnames(mat), function(x) strsplit(x, "_A")[[1]][1])
  
  celltype_colors <- setNames(RColorBrewer::brewer.pal(length(unique(celltype_vec)), "Dark2"),
                              unique(celltype_vec))
  
  ha <- HeatmapAnnotation(
    df = data.frame(group = groups, celltype = celltype_vec),
    col = list(group = c("CTRL"="grey40","CORT"="firebrick2","DEX"="steelblue3"),
               celltype = celltype_colors)
  )
  
  # =============================
  # 5) Dendrogram
  # =============================
  message("=== 5) Clustering dendrogram ===")
  
  dist_fun <- function(x) {
    as.dist(1 - cor(t(x), method = distance_method, use = "pairwise.complete.obs"))
  }
  
  dist_obj <- dist_fun(mat)
  hc <- hclust(dist_obj, method = linkage_method)
  row_dend <- as.dendrogram(hc)
  
  if (!is.null(cut_k)) {
    row_dend <- dendextend::color_branches(row_dend, k = cut_k)
  } else if (!is.null(cut_height)) {
    row_dend <- dendextend::color_branches(row_dend, h = cut_height)
  }
  
  # =============================
  # 6) HEATMAP
  # =============================
  message("=== 6) Draw heatmap ===")
  
  ht <- Heatmap(
    mat,
    name = "Z-score vs CTRL",
    cluster_rows = row_dend,
    cluster_columns = FALSE,
    top_annotation = ha,
    row_dend_width = row_dend_width,
    col = colorRamp2(zscore_range, c("navy", "white", "firebrick4")),
    heatmap_legend_param = list(
      title = "Z-score vs CTRL",
      at = zscore_range,
      labels = as.character(zscore_range)
    ),
    show_row_names = show_rows,
    show_column_names = show_cols
  )
  
  if (!is.null(output_path)) {
    svg(output_path, width = 11, height = 16)
    draw(ht, merge_legend = TRUE)
    dev.off()
  } else {
    draw(ht, merge_legend = TRUE)
  }
  
  return(mat)
}

generate_heatmap_test_clustering(
  combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
  celltypes,
  distance_method = "spearman",
  linkage_method = "ward.D2",
  logfc_mode = "both",
  # cut_k = 6,
  row_dend_width = unit(80,"mm"),
  zscore_range = c(-2,0,2)
)

