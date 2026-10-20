generate_heatmap_ctrlNormalized_zscore_per_block <- function(
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
    cluster_rows_k  = NULL
) {
  
  # Libraries
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(RColorBrewer)
  
  # match args
  logfc_mode      <- match.arg(logfc_mode)
  contrast        <- match.arg(contrast)
  distance_method <- match.arg(distance_method)
  linkage_method  <- match.arg(linkage_method)
  
  # ============================================
  # 1) GENE SELECTION
  # ============================================
  message("=== 1) Gene selection ===")
  
  genes2heatmap <- data_list[celltypes] %>%
    lapply(function(x){
      
      df <- x %>% filter(meanCPM_CORT > cpm_threshold | meanCPM_DEX > cpm_threshold)
      
      # LOGFC rules
      if (logfc_mode == "up") {
        df <- df %>% filter(log2FC_edgerCORT > logfc_threshold | log2FC_edgerDEX > logfc_threshold)
      } else if (logfc_mode == "down") {
        df <- df %>% filter(log2FC_edgerCORT < -logfc_threshold | log2FC_edgerDEX < -logfc_threshold)
      } else {
        df <- df %>% filter(abs(log2FC_edgerCORT) > logfc_threshold | abs(log2FC_edgerDEX) > logfc_threshold)
      }
      
      # P-values
      if (contrast == "CORT") {
        df <- df %>% filter(PValue_edgerCORT < pval_threshold)
      } else if (contrast == "DEX") {
        df <- df %>% filter(PValue_edgerDEX < pval_threshold)
      } else {
        df <- df %>% filter(PValue_edgerCORT < pval_threshold | PValue_edgerDEX < pval_threshold)
      }
      
      df$gene_symbol
    }) %>% unlist() %>% unique()
  
  message("Liczba wybranych genów: ", length(genes2heatmap))
  if (length(genes2heatmap) == 0) stop("❌ Brak genów do heatmapy.")
  
  
  # ============================================
  # 2) MERGE RAW QUANTILE-NORMALIZED VALUES
  # ============================================
  message("=== 2) Merge raw QN ===")
  
  df_list <- lapply(celltypes, function(cl) {
    df_raw <- data_list[[cl]] %>%
      filter(gene_symbol %in% genes2heatmap) %>%
      select(gene_symbol, paste0("QN_A", 1:12))
    
    df_out <- data.frame(df_raw)
    new_names <- c("gene_symbol",
                   paste0(gsub(" ", "_", cl), "_QN_A", 1:12))
    colnames(df_out) <- new_names
    df_out
  })
  
  merged <- Reduce(function(x, y) full_join(x, y, by = "gene_symbol"), df_list)
  rownames(merged) <- merged$gene_symbol
  merged$gene_symbol <- NULL
  
  mat <- as.matrix(merged)
  mat[!is.finite(mat)] <- 0
  
  
  # ============================================
  # 3) NORMALIZE TO CTRL MEAN
  # ============================================
  message("=== 3) Normalize to CTRL global mean ===")
  
  ctrl_cols <- grep("A([1-4])$", colnames(mat), value = TRUE)
  ctrl_mean <- rowMeans(mat[, ctrl_cols, drop = FALSE])
  
  mat_norm <- sweep(mat, 1, ctrl_mean, "/")
  mat_norm[!is.finite(mat_norm)] <- 0
  
  
  # ============================================
  # 4) BLOCK-WISE Z-SCORE
  # ============================================
  message("=== 4) Z-score per block ===")
  
  block_list <- lapply(celltypes, function(cl) {
    cols <- grep(paste0("^", gsub(" ", "_", cl), "_QN_A"), colnames(mat_norm), value = TRUE)
    block <- mat_norm[, cols, drop = FALSE]
    scaled <- t(scale(t(block)))
    scaled[!is.finite(scaled)] <- 0
    scaled
  })
  
  mat_z <- do.call(cbind, block_list)
  
  
  # ============================================
  # 5) FIX COLUMN ORDER (NO CLUSTERING)
  # ============================================
  message("=== 5) Ordering columns: CTRL → CORT → DEX ===")
  
  ctrl_cols <- grep("A([1-4])$",  colnames(mat_z), value = TRUE)
  dex_cols  <- grep("A(5|6|7|8)$", colnames(mat_z), value = TRUE)
  cort_cols <- grep("A(9|10|11|12)$", colnames(mat_z), value = TRUE)
  
  ordered_cols <- c(ctrl_cols, cort_cols, dex_cols)
  mat_z <- mat_z[, ordered_cols]
  
  groups <- ifelse(colnames(mat_z) %in% ctrl_cols, "CTRL",
                   ifelse(colnames(mat_z) %in% cort_cols, "CORT", "DEX"))
  
  celltype_vec <- sub("_QN_A.*", "", colnames(mat_z))
  
  ha <- HeatmapAnnotation(
    df = data.frame(
      group = groups,
      celltype = celltype_vec
    ),
    col = list(
      group = c("CTRL"="grey40", "CORT"="firebrick2", "DEX"="steelblue3"),
      celltype = setNames(brewer.pal(length(unique(celltype_vec)), "Dark2"),
                          unique(celltype_vec))
    )
  )
  
  
  # ============================================
  # 6) HEATMAP (HARD DISABLE COLUMN CLUSTERING)
  # ============================================
  message("=== 6) Heatmap ===")
  
  dist_fun <- function(x) as.dist(1 - cor(t(x), method = distance_method))
  
  # BLOKADA KLASTROWANIA KOLUMN
  zero_distance <- function(x) as.dist(matrix(0, length(x), length(x)))
  
  ht <- Heatmap(
    mat_z,
    name = "CTRL-normalized Z-score",
    top_annotation = ha,
    cluster_rows = TRUE,
    clustering_distance_rows = dist_fun,
    clustering_method_rows = linkage_method,
    
    # HARD-OFF
    cluster_columns = FALSE,
    clustering_distance_columns = zero_distance,
    clustering_method_columns = "single",
    
    row_km = cluster_rows_k,
    row_dend_width = row_dend_width,
    show_row_names = show_rows,
    show_column_names = show_cols,
    col = colorRamp2(zscore_range, c("navy", "white", "firebrick4")),
    
    column_title = paste0("CTRL-normalized + Z-score per block (n=", nrow(mat_z), ")")
  )
  
  draw(ht, merge_legend = TRUE)
  
  return(mat_z)
}


generate_heatmap_ctrlNormalized_zscore_per_block(
  combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
  celltypes,
  distance_method = "spearman",
  linkage_method = "complete",
  logfc_mode = "both",
  zscore_range = c(-2,0,2),
  cluster_rows_k = 5,
  pval_threshold = 0.00001,
  show_cols = TRUE,
  row_dend_width = unit(70, "mm")
)

generate_heatmap_ctrlNormalized_zscore_per_block(
  combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
  celltypes,
  distance_method = "spearman",
  linkage_method = "complete",
  logfc_mode = "both",
  zscore_range = c(-2,0,2),
  cluster_rows_k = 6,
  pval_threshold = 0.00001,
  show_cols = TRUE,
  row_dend_width = unit(20, "mm")
)


