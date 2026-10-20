generate_heatmap_zscore_global <- function(
    data_list,
    celltypes,
    distance_method = c("spearman", "pearson", "euclidean", "manhattan"),
    linkage_method  = c("ward.D2", "average", "complete", "single"),
    cpm_threshold   = 50,
    pval_threshold  = 0.01,
    logfc_threshold = 1,
    logfc_mode      = c("up", "down", "both"),
    contrast        = c("CORT", "DEX", "both"),
    zscore_range    = c(-4, 0, 4),    # <===== NOWY ARGUMENT
    show_rows       = TRUE,
    show_cols       = FALSE,
    row_dend_width  = unit(40, "mm"),
    cluster_rows_k  = NULL
) {
  
  library(dplyr); library(tidyr); library(ComplexHeatmap)
  library(circlize); library(grid); library(RColorBrewer)
  
  logfc_mode      <- match.arg(logfc_mode)
  contrast        <- match.arg(contrast)
  distance_method <- match.arg(distance_method)
  linkage_method  <- match.arg(linkage_method)
  
  # =============================
  # 1) WYBÓR GENÓW
  # =============================
  message("=== 1) Wybór genów ===")
  
  genes2heatmap <- data_list[celltypes] %>%
    lapply(function(x){
      df <- x %>% filter(meanCPM_CORT > cpm_threshold | meanCPM_DEX > cpm_threshold)
      
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
  # 2) MERGE RAW QN
  # =============================
  message("=== 2) Merge QN ===")
  
  df_list <- lapply(celltypes, function(cl) {
    df_raw <- data_list[[cl]] %>%
      filter(gene_symbol %in% genes2heatmap) %>%
      select(gene_symbol, paste0("QN_A", 1:12))
    
    data.frame(
      gene_symbol = df_raw$gene_symbol,
      df_raw[, -1]
    ) %>% setNames(c("gene_symbol", paste0(gsub(" ", "_", cl), "_", colnames(df_raw)[-1])))
  })
  
  merged <- Reduce(function(x, y) full_join(x, y, by = "gene_symbol"), df_list)
  rownames(merged) <- merged$gene_symbol
  merged$gene_symbol <- NULL
  
  mat <- as.matrix(merged)
  mat[!is.finite(mat)] <- 0
  
  # =============================
  # 3) GLOBAL Z-SCORE
  # =============================
  message("=== 3) Global Z-score ===")
  
  mat <- t(scale(t(mat)))
  mat[!is.finite(mat)] <- 0
  
  # =============================
  # 4) ORDER COLUMNS
  # =============================
  message("=== 4) Order CTRL → CORT → DEX ===")
  
  ctrl_cols <- grep("A([1-4])$", colnames(mat), value = TRUE)
  cort_cols <- grep("A(5|6|7|8)$", colnames(mat), value = TRUE)
  dex_cols  <- grep("A(9|1[0-2])$", colnames(mat), value = TRUE)
  
  mat <- mat[, c(ctrl_cols, cort_cols, dex_cols)]
  
  groups <- ifelse(colnames(mat) %in% ctrl_cols, "CTRL",
                   ifelse(colnames(mat) %in% cort_cols, "CORT", "DEX"))
  
  celltype_vec <- sapply(colnames(mat), function(x) strsplit(x, "_A")[[1]][1])
  
  ha <- HeatmapAnnotation(
    df = data.frame(group = groups, celltype = celltype_vec),
    col = list(
      group = c("CTRL"="grey40", "CORT"="firebrick2", "DEX"="steelblue3"),
      celltype = setNames(brewer.pal(length(unique(celltype_vec)), "Dark2"),
                          unique(celltype_vec))
    )
  )
  
  # =============================
  # 5) HEATMAP
  # =============================
  message("=== 5) Draw heatmap ===")
  
  dist_fun <- function(x) as.dist(1 - cor(t(x), method = distance_method))
  
  ht <- Heatmap(
    mat,
    name = "Global Z-score",
    top_annotation = ha,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    clustering_distance_rows = dist_fun,
    clustering_method_rows = linkage_method,
    row_km = cluster_rows_k,
    row_dend_width = row_dend_width,
    show_row_names = show_rows,
    show_column_names = show_cols,
    col = colorRamp2(zscore_range, c("navy", "white", "firebrick4")),  # <=== Użycie parametru
    column_title = paste0("Global Z-score / genes=", nrow(mat))
  )
  
  draw(ht, merge_legend = TRUE)
  return(mat)
}


generate_heatmap_zscore_global(
  combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
  celltypes,
  distance_method = "spearman",
  linkage_method = "complete",
  logfc_mode = "both",
  # cluster_rows_k = 5,
  row_dend_width = unit(70, "mm"),
  zscore_range = c(-2, 0, 2)
)

