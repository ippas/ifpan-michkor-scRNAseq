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
    zscore_range    = c(-4, 0, 4),   # <---- RANGE Z-SCORE
    show_rows       = TRUE,
    show_cols       = FALSE,
    output_path     = NULL
) {
  # =============================
  # ---- LIBRARIES --------------
  # =============================
  require(dplyr); require(tidyr); require(tibble); require(purrr)
  require(ComplexHeatmap); require(circlize); require(grid); require(RColorBrewer)
  
  logfc_mode  <- match.arg(logfc_mode)
  contrast    <- match.arg(contrast)
  distance_method <- match.arg(distance_method)
  linkage_method  <- match.arg(linkage_method)
  
  # =============================
  # ---- 1) WYBÓR GENÓW ---------
  # =============================
  message("=== 1) Wybór genów ===")
  
  genes2heatmap <- data_list[celltypes] %>%
    lapply(function(x){
      df <- x %>%
        filter(meanCPM_CORT > cpm_threshold | meanCPM_DEX > cpm_threshold)
      
      if (logfc_mode == "up") {
        df <- df %>% filter(log2FC_edgerCORT > logfc_threshold | log2FC_edgerDEX > logfc_threshold)
      } else if (logfc_mode == "down") {
        df <- df %>% filter(log2FC_edgerCORT < -logfc_threshold | log2FC_edgerDEX < -logfc_threshold)
      } else {
        df <- df %>% filter(abs(log2FC_edgerCORT) > logfc_threshold | abs(log2FC_edgerDEX) > logfc_threshold)
      }
      
      if (contrast == "CORT") {
        df <- df %>% filter(PValue_edgerCORT < pval_threshold)
      } else if (contrast == "DEX") {
        df <- df %>% filter(PValue_edgerDEX < pval_threshold)
      } else {
        df <- df %>% filter(PValue_edgerCORT < pval_threshold | PValue_edgerDEX < pval_threshold)
      }
      
      df$gene_symbol
    }) %>% unlist() %>% unique()
  
  message("Liczba genów spełniających kryteria: ", length(genes2heatmap))
  
  # =============================
  # ---- 2) NORMALIZACJA --------
  # =============================
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
    
    out <- data.frame(gene_symbol = df_raw$gene_symbol, as.data.frame(zscores))
    colnames(out)[-1] <- paste0(gsub(" ", "_", cl), "_", colnames(out)[-1])
    out
  })
  
  names(df_list) <- celltypes
  
  # =============================
  # ---- 3) FULL MERGE ----------
  # =============================
  message("=== 3) Full merge ===")
  
  merged <- purrr::reduce(df_list, full_join, by = "gene_symbol") %>%
    column_to_rownames("gene_symbol")
  
  mat <- as.matrix(merged)
  mat[!is.finite(mat)] <- 0
  
  # =============================
  # ---- 4) ORDER ---------------
  # =============================
  message("=== 4) Kolejność CTRL → CORT → DEX ===")
  
  ctrl_cols <- grep("A([1-4])$", colnames(mat), value = TRUE)
  cort_cols <- grep("A(5|6|7|8)$", colnames(mat), value = TRUE)
  dex_cols  <- grep("A(9|1[0-2])$", colnames(mat), value = TRUE)
  
  mat <- mat[, c(ctrl_cols, cort_cols, dex_cols)]
  
  groups <- ifelse(colnames(mat) %in% ctrl_cols, "CTRL",
                   ifelse(colnames(mat) %in% cort_cols, "CORT", "DEX"))
  
  celltype_vec <- sapply(colnames(mat), function(x) strsplit(x, "_A")[[1]][1])
  
  celltype_colors <- setNames(
    RColorBrewer::brewer.pal(length(unique(celltype_vec)), "Dark2"),
    unique(celltype_vec)
  )
  
  # =============================
  # ---- 5) ANNOTATION ----------
  # =============================
  message("=== 5) Annotation ===")
  
  ha <- HeatmapAnnotation(
    df = data.frame(
      group = groups,
      celltype = celltype_vec
    ),
    col = list(
      group = c("CTRL"="grey40", "CORT"="firebrick2", "DEX"="steelblue3"),
      celltype = celltype_colors
    )
  )
  
  # =============================
  # ---- 6) HEATMAP -------------
  # =============================
  message("=== 6) Draw heatmap ===")
  
  n_genes <- nrow(mat)
  
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
    show_column_names = show_cols,
    col = colorRamp2(zscore_range, c("navy", "white", "firebrick4")),  # <--- KOLOR
    row_title = paste0("Genes (n=", n_genes, ")"),
    row_title_gp = gpar(fontsize = 16, fontface = "bold"),
    column_title = paste0(
      "Heatmap | genes=", n_genes,
      " | dist=", distance_method,
      " | link=", linkage_method,
      " | CPM≥", cpm_threshold,
      " | p<", pval_threshold,
      " | logFC≥", logfc_threshold,
      " | zrange=", paste(zscore_range, collapse=",")
    ),
    column_title_gp = gpar(fontsize = 18, fontface = "bold")
  )
  
  # =============================
  # ---- 7) SAVE OR SHOW --------
  # =============================
  if (!is.null(output_path)) {
    message("=== Zapis do pliku SVG ===")
    svg(output_path, width = 11, height = 14)
    draw(ht, merge_legend = TRUE)
    dev.off()
    message(paste("Zapisano do pliku:", output_path))
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
  cpm_threshold   = 50,
  pval_threshold  = 0.01,
  logfc_threshold = 1,
  zscore_range = c(-2.5, 0, 2.5)
)
