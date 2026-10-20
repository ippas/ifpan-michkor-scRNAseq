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
    cpm_filter_mode = c("mean", "max"),   # <-- NOWY ARGUMENT
    zscore_range    = c(-4, 0, 4),
    show_rows       = TRUE,
    show_cols       = FALSE,
    row_dend_width  = grid::unit(5, "mm"),
    cluster_rows_k  = NULL
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
      # tak jak było wcześniej
      df <- df %>%
        filter(meanCPM_CORT > cpm_threshold |
                 meanCPM_DEX  > cpm_threshold)
      
    } else if (cpm_filter_mode == "max") {
      # filtr TYLKO przez max_CPM
      if (!("max_CPM" %in% colnames(df))) {
        stop("❌ cpm_filter_mode = 'max', ale w data_list[['", cl, "']] brak kolumny 'max_CPM'")
      }
      df <- df %>%
        filter(max_CPM > cpm_threshold)
    }
    
    ## ---- logFC FILTER ----
    if (logfc_mode == "up") {
      df <- df %>%
        filter(log2FC_edgerCORT > logfc_threshold |
                 log2FC_edgerDEX  > logfc_threshold)
      
    } else if (logfc_mode == "down") {
      df <- df %>%
        filter(log2FC_edgerCORT < -logfc_threshold |
                 log2FC_edgerDEX  < -logfc_threshold)
      
    } else { # "both"
      df <- df %>%
        filter(abs(log2FC_edgerCORT) > logfc_threshold |
                 abs(log2FC_edgerDEX)  > logfc_threshold)
    }
    
    ## ---- p-value FILTER ----
    if (contrast == "CORT") {
      df <- df %>% filter(PValue_edgerCORT < pval_threshold)
    } else if (contrast == "DEX") {
      df <- df %>% filter(PValue_edgerDEX < pval_threshold)
    } else { # both
      df <- df %>%
        filter(PValue_edgerCORT < pval_threshold |
                 PValue_edgerDEX  < pval_threshold)
    }
    
    gene_list[[cl]] <- df$gene_symbol
  }
  
  genes2heatmap <- unique(unlist(gene_list))
  n_genes_selected <- length(genes2heatmap)
  
  if (n_genes_selected == 0)
    stop("❌ No genes after filtering (CPM/logFC/p-value)")
  
  message("Genes selected after CPM/logFC/p filtering: ", n_genes_selected)
  
  ## ============================================================
  ## 2. MERGE RAW QN VALUES
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
  ## 3. NORMALIZE TO CTRL
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
  
  message("Genes in final matrix (rows of mat_z): ", nrow(mat_z))
  
  ## ============================================================
  ## 5. ORDER COLUMNS CTRL → CORT → DEX
  ## ============================================================
  message("=== 5) Ordering columns ===")
  
  ctrl_cols <- grep("A([1-4])$",      colnames(mat_z), value = TRUE)
  cort_cols <- grep("A(9|10|11|12)$", colnames(mat_z), value = TRUE)
  dex_cols  <- grep("A(5|6|7|8)$",    colnames(mat_z), value = TRUE)
  
  mat_z <- mat_z[, c(ctrl_cols, cort_cols, dex_cols)]
  
  groups <- ifelse(colnames(mat_z) %in% ctrl_cols, "CTRL",
                   ifelse(colnames(mat_z) %in% cort_cols, "CORT", "DEX"))
  
  celltype_vec <- sub("_QN_A.*", "", colnames(mat_z))
  
  ## ============================================================
  ## 6. CHECK LOST GENES (między filtrem a macierzą)
  ## ============================================================
  all_genes_selected <- genes2heatmap
  genes_in_matrix    <- rownames(mat_z)
  lost_genes <- setdiff(all_genes_selected, genes_in_matrix)
  
  if (length(lost_genes) > 0) {
    message("⚠️ WARNING: Some genes were dropped before heatmap matrix!")
    message("   Genes selected after filtering: ", length(all_genes_selected))
    message("   Genes in final matrix:          ", length(genes_in_matrix))
    message("   Genes lost: ", length(lost_genes))
    message("   Lost genes (first 20): ",
            paste(head(lost_genes, 20), collapse = ", "))
    if (length(lost_genes) > 20)
      message("   ... and ", length(lost_genes) - 20, " more.")
  }
  
  ## ============================================================
  ## 7. HEATMAP ANNOTATION
  ## ============================================================
  ha <- HeatmapAnnotation(
    df = data.frame(group = groups, celltype = celltype_vec),
    col = list(
      group = c("CTRL" = "grey40",
                "CORT" = "firebrick2",
                "DEX"  = "steelblue3"),
      celltype = setNames(
        brewer.pal(n = length(unique(celltype_vec)), "Dark2"),
        unique(celltype_vec)
      )
    )
  )
  
  ## ============================================================
  ## 8. DRAW HEATMAP
  ## ============================================================
  message("=== 7) Drawing heatmap ===")
  
  dist_fun  <- function(x) as.dist(1 - cor(t(x), method = distance_method))
  zero_dist <- function(x) as.dist(matrix(0, length(x), length(x)))
  
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
    
    column_split = factor(groups, levels = c("CTRL", "CORT", "DEX")),
    
    col = colorRamp2(zscore_range, c("navy","white","firebrick4")),
    
    column_title = paste0(
      "CTRL-normalized Z-score | genes selected = ", n_genes_selected, "\n",
      "Distance: ", distance_method,
      " | Linkage: ", linkage_method,
      " | p < ", pval_threshold,
      " | CPM > ", cpm_threshold,
      " | CPM filter mode = ", cpm_filter_mode
    ),
    
    show_row_names    = show_rows,
    show_column_names = show_cols
  )
  
  ht <- draw(ht, merge_legend = TRUE)
  
  ## ============================================================
  ## 9. OPTIONAL: CLUSTER INFO (jeśli używasz row_km)
  ## ============================================================
  cluster_lists      <- NULL
  cluster_assignment <- NULL
  clusters_df        <- NULL
  
  if (!is.null(cluster_rows_k)) {
    
    clust <- row_order(ht)  # lista indeksów dla każdego klastra
    
    cluster_lists <- lapply(clust, function(idx) rownames(mat_z)[idx])
    names(cluster_lists) <- paste0("Cluster_", seq_along(cluster_lists))
    
    cluster_assignment <- rep(NA_integer_, nrow(mat_z))
    for (i in seq_along(clust)) {
      cluster_assignment[clust[[i]]] <- i
    }
    
    clusters_df <- data.frame(
      gene    = rownames(mat_z),
      cluster = cluster_assignment
    )
    
    all_meta <- bind_rows(data_list) %>%
      distinct(gene_symbol, .keep_all = TRUE)
    
    clusters_df <- clusters_df %>%
      left_join(all_meta, by = c("gene" = "gene_symbol"))
  }
  
  ## ============================================================
  ## 10. RETURN
  ## ============================================================
  return(list(
    heatmap_matrix     = mat_z,
    genes_selected     = genes2heatmap,
    genes_in_matrix    = genes_in_matrix,
    lost_genes         = lost_genes,
    metadata = list(
      n_genes_selected = n_genes_selected,
      n_genes_matrix   = length(genes_in_matrix),
      distance_method  = distance_method,
      linkage_method   = linkage_method,
      cpm_threshold    = cpm_threshold,
      pval_threshold   = pval_threshold,
      logfc_threshold  = logfc_threshold,
      logfc_mode       = logfc_mode,
      contrast         = contrast,
      cpm_filter_mode  = cpm_filter_mode
    ),
    cluster_assignment = cluster_assignment,
    cluster_genes      = cluster_lists,
    clusters_df        = clusters_df
  ))
}
