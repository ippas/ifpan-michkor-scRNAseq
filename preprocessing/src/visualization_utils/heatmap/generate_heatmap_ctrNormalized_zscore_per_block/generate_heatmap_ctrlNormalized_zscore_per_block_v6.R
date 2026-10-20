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
    cpm_filter_mode = c("mean", "max"),
    statistics_mode = c("pvalue", "fdr"),
    zscore_range    = c(-4, 0, 4),
    show_rows       = TRUE,
    show_cols       = FALSE,
    row_dend_width  = grid::unit(5, "mm"),
    cluster_rows_k  = NULL,
    ## ==== NOWE ARGUMENTY ====
    svg_filename    = NULL,   # <- nowa nazwa argumentu
    width           = 14,
    height          = 10
) {
  
  ## ==== LIBS ====
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(RColorBrewer)
  library(purrr)
  
  ## ==== MATCH ARGS ====
  logfc_mode      <- match.arg(logfc_mode)
  contrast        <- match.arg(contrast)
  cpm_filter_mode <- match.arg(cpm_filter_mode)
  statistics_mode <- match.arg(statistics_mode)
  distance_method <- match.arg(distance_method)
  linkage_method  <- match.arg(linkage_method)
  
  ## ============================================================
  ## 1. SELECT GENES
  ## ============================================================
  message("=== 1) Gene selection ===")
  
  gene_list <- list()
  
  for (cl in celltypes) {
    
    df <- data_list[[cl]]
    
    ## ---- CPM FILTER ----
    if (cpm_filter_mode == "mean") {
      df <- df %>%
        filter(meanCPM_CORT > cpm_threshold |
                 meanCPM_DEX  > cpm_threshold)
      
    } else if (cpm_filter_mode == "max") {
      if (!("max_CPM" %in% colnames(df))) {
        stop("❌ Missing 'max_CPM' for cpm_filter_mode = 'max'.")
      }
      df <- df %>% filter(max_CPM > cpm_threshold)
    }
    
    ## ---- logFC FILTER ----
    if (logfc_mode == "up") {
      df <- df %>% filter(log2FC_edgerCORT > logfc_threshold |
                            log2FC_edgerDEX  > logfc_threshold)
      
    } else if (logfc_mode == "down") {
      df <- df %>% filter(log2FC_edgerCORT < -logfc_threshold |
                            log2FC_edgerDEX  < -logfc_threshold)
      
    } else {
      df <- df %>% filter(abs(log2FC_edgerCORT) > logfc_threshold |
                            abs(log2FC_edgerDEX)  > logfc_threshold)
    }
    
    ## ============================================================
    ## ---- STATISTICS MODE: pvalue / fdr ----
    ## ============================================================
    if (statistics_mode == "pvalue") {
      colC <- "PValue_edgerCORT"
      colD <- "PValue_edgerDEX"
    } else {
      colC <- "FDR_edgerCORT"
      colD <- "FDR_edgerDEX"
    }
    
    if (contrast == "CORT") {
      df <- df %>% filter(.data[[colC]] < pval_threshold)
    } else if (contrast == "DEX") {
      df <- df %>% filter(.data[[colD]] < pval_threshold)
    } else {
      df <- df %>% filter(.data[[colC]] < pval_threshold |
                            .data[[colD]] < pval_threshold)
    }
    
    gene_list[[cl]] <- df$gene_symbol
  }
  
  genes2heatmap <- unique(unlist(gene_list))
  
  if (length(genes2heatmap) == 0)
    stop("❌ No genes passed filtering.")
  
  n_genes_selected <- length(genes2heatmap)
  message("Genes selected: ", n_genes_selected)
  
  ## ============================================================
  ## 2. MERGE QN VALUES
  ## ============================================================
  message("=== 2) Merge raw QN ===")
  
  df_list <- lapply(celltypes, function(cl) {
    df <- data_list[[cl]] %>%
      filter(gene_symbol %in% genes2heatmap) %>%
      select(gene_symbol, paste0("QN_A", 1:12))
    
    colnames(df) <- c("gene_symbol",
                      paste0(gsub(" ", "_", cl), "_QN_A", 1:12))
    df
  })
  
  merged <- Reduce(function(a, b) full_join(a, b, by = "gene_symbol"), df_list)
  rownames(merged) <- merged$gene_symbol
  merged$gene_symbol <- NULL
  
  mat <- as.matrix(merged)
  mat[!is.finite(mat)] <- 0
  
  ## ============================================================
  ## 3. CTRL NORMALIZATION
  ## ============================================================
  message("=== 3) Normalize to CTRL ===")
  
  ctrl_cols <- grep("A([1-4])$", colnames(mat), value = TRUE)
  ctrl_mean <- rowMeans(mat[, ctrl_cols, drop = FALSE])
  
  mat_norm <- sweep(mat, 1, ctrl_mean, "/")
  mat_norm[!is.finite(mat_norm)] <- 0
  
  ## ============================================================
  ## 4. Z-SCORE PER BLOCK
  ## ============================================================
  message("=== 4) Z-score per block ===")
  
  block_list <- lapply(celltypes, function(cl) {
    cols <- grep(paste0("^", gsub(" ", "_", cl), "_QN_A"),
                 colnames(mat_norm), value = TRUE)
    z <- t(scale(t(mat_norm[, cols, drop = FALSE])))
    z[!is.finite(z)] <- 0
    z
  })
  
  mat_z <- do.call(cbind, block_list)
  
  message("Genes in matrix: ", nrow(mat_z))
  
  ## ============================================================
  ## 5. ORDER COLUMNS
  ## ============================================================
  message("=== 5) Ordering columns ===")
  
  ctrl_cols <- grep("A([1-4])$", colnames(mat_z), value = TRUE)
  cort_cols <- grep("A(9|10|11|12)$", colnames(mat_z), value = TRUE)
  dex_cols  <- grep("A(5|6|7|8)$",    colnames(mat_z), value = TRUE)
  
  mat_z <- mat_z[, c(ctrl_cols, cort_cols, dex_cols)]
  
  groups <- ifelse(colnames(mat_z) %in% ctrl_cols, "CTRL",
                   ifelse(colnames(mat_z) %in% cort_cols, "CORT", "DEX"))
  celltype_vec <- sub("_QN_A.*", "", colnames(mat_z))
  
  ## ============================================================
  ## 6. ANNOTATION
  ## ============================================================
  ha <- HeatmapAnnotation(
    df = data.frame(group = groups, celltype = celltype_vec),
    col = list(
      group = c("CTRL"="grey40", "CORT"="firebrick2", "DEX"="steelblue3"),
      celltype = setNames(
        brewer.pal(n = length(unique(celltype_vec)), "Dark2"),
        unique(celltype_vec)
      )
    )
  )
  
  ## ============================================================
  ## 7. HEATMAP OBJECT
  ## ============================================================
  message("=== 6) Drawing heatmap ===")
  
  dist_fun  <- function(x) as.dist(1 - cor(t(x), method = distance_method))
  zero_dist <- function(x) as.dist(matrix(0, nrow(x), nrow(x)))
  
  ht <- Heatmap(
    mat_z,
    name = "CTRL-normalized Z-score",
    top_annotation = ha,
    
    cluster_rows = TRUE,
    clustering_distance_rows = dist_fun,
    clustering_method_rows   = linkage_method,
    row_km        = cluster_rows_k,
    row_dend_width = row_dend_width,
    
    cluster_columns = FALSE,
    clustering_distance_columns = zero_dist,
    column_split = factor(groups, levels=c("CTRL","CORT","DEX")),
    
    col = colorRamp2(zscore_range, c("navy","white","firebrick4")),
    
    column_title = paste0(
      "CTRL-normalized Z-score | genes = ", n_genes_selected, "\n",
      "Distance: ", distance_method,
      " | Linkage: ", linkage_method,
      " | ", toupper(statistics_mode), " < ", pval_threshold,
      " | CPM > ", cpm_threshold
    ),
    
    show_row_names    = show_rows,
    show_column_names = show_cols
  )
  
  ## ============================================================
  ## 8. EXPORT TO SVG
  ## ============================================================
  if (!is.null(svg_filename)) {
    message("=== Saving heatmap to SVG: ", svg_filename)
    svg(svg_filename, width = width, height = height)
    draw(ht, merge_legend = TRUE)
    dev.off()
  } else {
    draw(ht, merge_legend = TRUE)
  }
  
  ## ============================================================
  ## 9. RETURN
  ## ============================================================
  return(list(
    heatmap_matrix = mat_z,
    genes_selected = genes2heatmap,
    metadata = list(
      n_genes_selected = n_genes_selected,
      svg_filename     = svg_filename,
      width            = width,
      height           = height,
      statistics_mode  = statistics_mode
    )
  ))
}
