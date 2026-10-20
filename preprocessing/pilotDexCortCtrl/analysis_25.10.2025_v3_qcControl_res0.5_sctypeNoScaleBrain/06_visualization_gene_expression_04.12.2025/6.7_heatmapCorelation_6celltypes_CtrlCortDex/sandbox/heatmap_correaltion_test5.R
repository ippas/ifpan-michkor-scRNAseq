filter_genes_edger_by_celltype <- function(
    data_list,
    celltypes = names(data_list),
    genes = NULL,
    cpm_threshold   = 50,
    logfc_threshold = 1,
    logfc_mode      = c("up", "down", "both"),
    pval_threshold  = 0.01,
    statistics_mode = c("pvalue", "fdr"),
    contrast        = c("CORT", "DEX", "both"),
    cpm_filter_mode = c("mean", "max"),
    cor_method      = c("spearman", "pearson", "kendall"),
    cor_use         = c("pairwise.complete.obs", "complete.obs", "everything")
) {
  
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  
  logfc_mode      <- match.arg(logfc_mode)
  statistics_mode <- match.arg(statistics_mode)
  contrast        <- match.arg(contrast)
  cpm_filter_mode <- match.arg(cpm_filter_mode)
  
  if (!all(cor_method %in% c("spearman", "pearson", "kendall")))
    stop("Invalid cor_method")
  if (!all(cor_use %in% c("pairwise.complete.obs", "complete.obs", "everything")))
    stop("Invalid cor_use")
  
  ## ============================================================
  ## FILTER ONE CELLTYPE
  ## ============================================================
  filter_one <- function(df) {
    
    if (!is.null(genes))
      df <- df %>% filter(gene_symbol %in% genes)
    
    if (cpm_filter_mode == "mean") {
      df <- df %>%
        filter(meanCPM_CORT > cpm_threshold |
                 meanCPM_DEX  > cpm_threshold)
    } else {
      if (!"max_CPM" %in% colnames(df))
        stop("Column 'max_CPM' required for cpm_filter_mode = 'max'")
      df <- df %>% filter(max_CPM > cpm_threshold)
    }
    
    if (logfc_mode == "up") {
      df <- df %>%
        filter(log2FC_edgerCORT >  logfc_threshold |
                 log2FC_edgerDEX  >  logfc_threshold)
    } else if (logfc_mode == "down") {
      df <- df %>%
        filter(log2FC_edgerCORT < -logfc_threshold |
                 log2FC_edgerDEX  < -logfc_threshold)
    } else {
      df <- df %>%
        filter(abs(log2FC_edgerCORT) > logfc_threshold |
                 abs(log2FC_edgerDEX)  > logfc_threshold)
    }
    
    colC <- ifelse(statistics_mode == "pvalue",
                   "PValue_edgerCORT", "FDR_edgerCORT")
    colD <- ifelse(statistics_mode == "pvalue",
                   "PValue_edgerDEX",  "FDR_edgerDEX")
    
    if (contrast == "CORT") {
      df <- df %>% filter(.data[[colC]] < pval_threshold)
    } else if (contrast == "DEX") {
      df <- df %>% filter(.data[[colD]] < pval_threshold)
    } else {
      df <- df %>% filter(.data[[colC]] < pval_threshold |
                            .data[[colD]] < pval_threshold)
    }
    
    df
  }
  
  ## ============================================================
  ## 1. FILTER PER CELLTYPE
  ## ============================================================
  per_celltype_filtered <- lapply(celltypes, function(ct) {
    filter_one(data_list[[ct]])
  })
  names(per_celltype_filtered) <- celltypes
  
  ## ============================================================
  ## 2. GLOBAL GENE SET
  ## ============================================================
  global_genes <- sort(unique(unlist(
    lapply(per_celltype_filtered, \(x) x$gene_symbol)
  )))
  
  ## ============================================================
  ## 3. RE-FILTER INPUT DATA
  ## ============================================================
  per_celltype_allGenes <- lapply(celltypes, function(ct) {
    data_list[[ct]] %>% filter(gene_symbol %in% global_genes)
  })
  names(per_celltype_allGenes) <- celltypes
  
  ## ============================================================
  ## 4. PREPARE LONG DATA
  ## ============================================================
  df_long <- per_celltype_allGenes %>%
    lapply(\(x)
           x %>% select(gene_symbol, log2FC_CORT_CPM, log2FC_DEX_CPM)
    ) %>%
    bind_rows(.id = "celltype") %>%
    pivot_longer(
      cols = c(log2FC_CORT_CPM, log2FC_DEX_CPM),
      names_to  = "treatment",
      values_to = "log2FC"
    ) %>%
    mutate(
      treatment = recode(treatment,
                         log2FC_CORT_CPM = "CORT",
                         log2FC_DEX_CPM  = "DEX"),
      celltype  = str_replace_all(celltype, " ", "_"),
      group     = paste0(celltype, "-", treatment)
    )
  
  ## ============================================================
  ## 5. MATRICES
  ## ============================================================
  make_mat <- function(df) {
    df %>%
      select(-celltype, -treatment) %>%
      pivot_wider(names_from = group, values_from = log2FC) %>%
      column_to_rownames("gene_symbol")
  }
  
  mat_DEX  <- df_long %>% filter(treatment == "DEX")  %>% make_mat()
  mat_CORT <- df_long %>% filter(treatment == "CORT") %>% make_mat()
  mat_BOTH <- df_long %>% make_mat()
  
  mats_by_contrast <- list(
    DEX  = mat_DEX,
    CORT = mat_CORT,
    both = mat_BOTH
  )
  
  ## ============================================================
  ## 6. CORRELATIONS + BETWEEN TREATMENT
  ## ============================================================
  correlation_matrices <- list()
  
  for (m in cor_method) {
    correlation_matrices[[m]] <- list()
    
    for (u in cor_use) {
      
      tmp <- lapply(mats_by_contrast, \(mat)
                    cor(mat, method = m, use = u)
      )
      
      cor_betweenTreatment <- tmp$both %>%
        as.data.frame() %>%
        rownames_to_column("rowname") %>%
        pivot_longer(-rowname, names_to = "colname", values_to = "correlation") %>%
        separate(rowname, into = c("celltype_row", "treatment_row"), sep = "-") %>%
        separate(colname, into = c("celltype_col", "treatment_col"), sep = "-") %>%
        filter(celltype_row == celltype_col,
               treatment_row != treatment_col) %>%
        mutate(
          treatment_pair = paste(
            pmin(treatment_row, treatment_col),
            pmax(treatment_row, treatment_col),
            sep = "_"
          )
        ) %>%
        distinct(celltype_row, treatment_pair, .keep_all = TRUE) %>%
        select(
          celltype = celltype_row,
          treatment_row,
          treatment_col,
          correlation
        )
      
      tmp$cor_betweenTreatment <- cor_betweenTreatment
      correlation_matrices[[m]][[u]] <- tmp
    }
  }
  
  ## ============================================================
  ## 7. RETURN
  ## ============================================================
  list(
    per_celltype_filtered = per_celltype_filtered,
    per_celltype_allGenes = per_celltype_allGenes,
    global_genes          = global_genes,
    correlation_matrices  = correlation_matrices
  )
}



filter_genes_edger_by_celltype(
  data    = combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
  genes = res_heatmap$genes_per_cluster$cluster_5,
  celltypes    = c(
    "Radial glial cells",
    "Astrocytes",
    "OPCs",
    "Oligodendrocytes",
    "Immature neurons",
    "Mature neurons"
  ),
  cpm_threshold   = 50,
  logfc_threshold = 1,
  logfc_mode      = "both",
  pval_threshold  = 0.001,
  statistics_mode = "fdr",
  contrast        = "both",
  cpm_filter_mode = "max",
  cor_method      = c("spearman", "pearson", "kendall"),
  cor_use         = c("pairwise.complete.obs", "complete.obs", "everything")
) -> correlationGeneExpression_6cellTypes


correlationGeneExpression_6cellTypes$correlation_matrices$spearman$pairwise.complete.obs$cor_betweenTreatment

plot_correlation_heatmap <- function(
    cor_matrix,
    title        = NULL,
    cluster      = TRUE,
    add_text     = FALSE,
    digits       = 2,
    text_size    = 7,
    limits       = c(-1, 1),
    colors       = c("#2166ac", "white", "#b2182b"),
    legend_name  = "Correlation"
) {
  
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  
  if (!is.matrix(cor_matrix))
    stop("cor_matrix must be a matrix")
  
  # ============================================================
  # 1. CLUSTERING (FULL MATRIX)
  # ============================================================
  if (cluster) {
    hc <- hclust(dist(1 - cor_matrix), method = "average")
    row_dend <- as.dendrogram(hc)
    col_dend <- row_dend
  } else {
    row_dend <- FALSE
    col_dend <- FALSE
  }
  
  # ============================================================
  # 2. COLOR SCALE
  # ============================================================
  col_fun <- colorRamp2(
    c(limits[1], 0, limits[2]),
    colors
  )
  
  # ============================================================
  # 3. HEATMAP OBJECT
  # ============================================================
  ht <- Heatmap(
    cor_matrix,
    name = legend_name,
    col  = col_fun,
    cluster_rows    = row_dend,
    cluster_columns = col_dend,
    show_row_names    = TRUE,
    show_column_names = TRUE,
    row_names_side    = "left",
    column_names_side = "top",
    column_title = title,
    heatmap_legend_param = list(
      direction = "horizontal",
      title_position = "topcenter"
    ),
    cell_fun = if (add_text) {
      function(j, i, x, y, width, height, fill) {
        grid.text(
          formatC(cor_matrix[i, j], digits = digits, format = "f"),
          x, y,
          gp = gpar(fontsize = text_size)
        )
      }
    } else {
      NULL
    }
  )
  
  # ============================================================
  # 4. DRAW WITH LEGEND AT BOTTOM
  # ============================================================
  draw(
    ht,
    heatmap_legend_side = "bottom"
  )
}



plot_correlation_heatmap(
  correlationGeneExpression_6cellTypes$
    correlation_matrices$pearson$complete.obs$DEX,
  title    = "CORT | Pearson",
  add_text = TRUE,
  digits   = 2,
  text_size = 10
)


correlationGeneExpression_6cellTypes$correlation_matrices$spearman$pairwise.complete.obs$both %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  pivot_longer(
    cols = -rowname,
    names_to = "colname",
    values_to = "correlation"
  ) %>% 
  separate(
    rowname,
    into = c("row_celltype", "row_treatment"),
    sep  = "-"
  ) %>% 
  separate(
    colname,
    into = c("col_celltype", "col_treatment"),
    sep  = "-"
  ) %>% 
  filter(row_celltype == col_celltype) %>% 
  filter(row_treatment != col_treatment) %>% 
  filter(row_treatment != col_treatment) %>%
  mutate(
    treatment_pair = paste(
      pmin(row_treatment, col_treatment),
      pmax(row_treatment, col_treatment),
      sep = "_"
    )
  ) %>%
  distinct(row_celltype, treatment_pair, .keep_all = TRUE) %>%
  select(-treatment_pair) %>% 
  select(-col_celltype) %>% 
  dplyr::rename(celltype = row_celltype)


