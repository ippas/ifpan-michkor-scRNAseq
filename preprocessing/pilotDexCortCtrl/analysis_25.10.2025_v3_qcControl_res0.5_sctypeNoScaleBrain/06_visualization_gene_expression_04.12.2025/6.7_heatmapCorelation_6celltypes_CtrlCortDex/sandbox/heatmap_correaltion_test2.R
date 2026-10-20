# ##############################################################################
# ---- function ----
# ##############################################################################
filter_genes_edger_by_celltype <- function(
    data_list,
    celltypes = names(data_list),
    genes = NULL,   # OPTIONAL: restrict analysis to a gene subset
    cpm_threshold   = 50,
    logfc_threshold = 1,
    logfc_mode      = c("up", "down", "both"),
    pval_threshold  = 0.01,
    statistics_mode = c("pvalue", "fdr"),
    contrast        = c("CORT", "DEX", "both"),
    cpm_filter_mode = c("mean", "max")
) {
  
  ## ============================================================
  ## LIBS
  ## ============================================================
  library(dplyr)
  
  ## ============================================================
  ## MATCH ARGS
  ## ============================================================
  logfc_mode      <- match.arg(logfc_mode)
  statistics_mode <- match.arg(statistics_mode)
  contrast        <- match.arg(contrast)
  cpm_filter_mode <- match.arg(cpm_filter_mode)
  
  ## ============================================================
  ## INTERNAL: FULL FILTER FOR ONE CELLTYPE
  ## ============================================================
  filter_one <- function(df) {
    
    ## ---- OPTIONAL GENE SUBSET (ONLY IF PROVIDED) ----
    if (!is.null(genes)) {
      df <- df %>% filter(gene_symbol %in% genes)
    }
    
    ## ---- CPM FILTER ----
    if (cpm_filter_mode == "mean") {
      df <- df %>%
        filter(
          meanCPM_CORT > cpm_threshold |
            meanCPM_DEX  > cpm_threshold
        )
    } else {
      if (!("max_CPM" %in% colnames(df)))
        stop("❌ Column 'max_CPM' required for cpm_filter_mode = 'max'")
      df <- df %>% filter(max_CPM > cpm_threshold)
    }
    
    ## ---- log2FC FILTER ----
    if (logfc_mode == "up") {
      df <- df %>%
        filter(
          log2FC_edgerCORT >  logfc_threshold |
            log2FC_edgerDEX  >  logfc_threshold
        )
    } else if (logfc_mode == "down") {
      df <- df %>%
        filter(
          log2FC_edgerCORT < -logfc_threshold |
            log2FC_edgerDEX  < -logfc_threshold
        )
    } else {
      df <- df %>%
        filter(
          abs(log2FC_edgerCORT) > logfc_threshold |
            abs(log2FC_edgerDEX)  > logfc_threshold
        )
    }
    
    ## ---- STATISTICS FILTER ----
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
      df <- df %>%
        filter(
          .data[[colC]] < pval_threshold |
            .data[[colD]] < pval_threshold
        )
    }
    
    df
  }
  
  ## ============================================================
  ## 1. FULL FILTER PER CELLTYPE
  ## ============================================================
  per_celltype_filtered <- lapply(celltypes, function(ct) {
    filter_one(data_list[[ct]])
  })
  names(per_celltype_filtered) <- celltypes
  
  ## ============================================================
  ## 2. GLOBAL GENE SET (>= 1 CELLTYPE)
  ## ============================================================
  global_genes <- sort(unique(unlist(
    lapply(per_celltype_filtered, function(df) df$gene_symbol)
  )))
  
  ## ============================================================
  ## 3. RE-FILTER INPUT DATA ONLY BY GLOBAL GENES
  ## ============================================================
  per_celltype_allGenes <- lapply(celltypes, function(ct) {
    data_list[[ct]] %>%
      filter(gene_symbol %in% global_genes)
  })
  names(per_celltype_allGenes) <- celltypes
  
  ## ============================================================
  ## 3. PREPARE DATA TO CORRELATION
  ## ============================================================
  per_celltype_allGenes %>% 
    lapply(., function(x){
      x %>% dplyr::select(gene_symbol, log2FC_CORT_CPM, log2FC_DEX_CPM, log2FC_edgerCORT, log2FC_edgerCORT, 
                          PValue_edgerCORT, FDR_edgerCORT, PValue_edgerDEX, FDR_edgerDEX)
    }) %>% 
    bind_rows(., .id = "celltype") %>% 
    select(gene_symbol, log2FC_CORT_CPM, log2FC_DEX_CPM, celltype) %>% 
    pivot_longer(
      cols      = c(log2FC_CORT_CPM, log2FC_DEX_CPM),
      names_to  = "treatment",
      values_to = "log2FC"
    ) %>% 
    mutate(
      treatment = recode(
        treatment,
        log2FC_CORT_CPM = "CORT",
        log2FC_DEX_CPM  = "DEX"
      )
    ) %>%  
    mutate(celltype = str_replace_all(celltype, " ", "_")) %>% 
    mutate(group = paste0(celltype, "-", treatment)) -> df_preprocessing_to_corr
  
  # ---- DEX ----
  df_preprocessing_to_corr %>% 
    filter(treatment == "DEX") %>% 
    select(-c(celltype, treatment)) %>% 
    pivot_wider(
      names_from  = group,
      values_from = log2FC
    ) %>%
    column_to_rownames("gene_symbol") -> df_wide_toCorrelation
    
  corr_DEX <- cor(df_wide_toCorrelation, use = "pairwise.complete.obs", method = "spearman")
  
  # ---- CORT ----
  df_preprocessing_to_corr %>% 
    filter(treatment == "CORT") %>% 
    select(-c(celltype, treatment)) %>% 
    pivot_wider(
      names_from  = group,
      values_from = log2FC
    ) %>%
    column_to_rownames("gene_symbol") -> df_wide_toCorrelation
  
  corr_CORT <- cor(df_wide_toCorrelation, use = "pairwise.complete.obs", method = "spearman")
  
  # ---- BOTH ----
  df_preprocessing_to_corr %>% 
    select(-c(celltype, treatment)) %>% 
    pivot_wider(
      names_from  = group,
      values_from = log2FC
    ) %>%
    column_to_rownames("gene_symbol") -> df_wide_toCorrelation
  
  corr_both <- cor(df_wide_toCorrelation, use = "pairwise.complete.obs", method = "spearman")
  
  ## ============================================================
  ## 4. RETURN
  ## ============================================================
  return(list(
    per_celltype_filtered = per_celltype_filtered,
    per_celltype_allGenes = per_celltype_allGenes,
    global_genes          = global_genes,
    correlation_matrices  = list(
      corr_CORT = corr_CORT,
      corr_DEX  = corr_DEX,
      corr_both = corr_both
    ),
    metadata = list(
      celltypes       = celltypes,
      n_global_genes  = length(global_genes),
      cpm_threshold   = cpm_threshold,
      cpm_filter_mode = cpm_filter_mode,
      logfc_threshold = logfc_threshold,
      logfc_mode      = logfc_mode,
      pval_threshold  = pval_threshold,
      statistics_mode = statistics_mode,
      contrast        = contrast,
      genes_subset    = !is.null(genes)
    )
  ))
}

filter_genes_edger_by_celltype(
  data    = combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
  # genes = res_heatmap$genes_per_cluster$cluster_1,
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
  cpm_filter_mode = "max"
) -> correlationGeneExpression_6cellTypes


correlationGeneExpression_6cellTypes$correlation_matrices$corr_CORT
