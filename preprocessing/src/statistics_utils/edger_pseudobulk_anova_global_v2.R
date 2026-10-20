edger_pseudobulk_anova_global_v2 <- function(
    agg,                      # matrix (genes x samples) OR list of matrices (clusters)
    sample_meta,
    group_col            = "treatment",
    sample_col           = "sample_id",
    group_levels         = c("CTRL","DEX","CORT"),
    robust               = TRUE,
    verbose              = TRUE,
    FDR_global_threshold = 0.05
) {
  if (!requireNamespace("edgeR", quietly = TRUE)) stop("Package 'edgeR' is required.")
  if (!requireNamespace("limma", quietly = TRUE)) stop("Package 'limma' is required.")
  
  .msg <- function(...) if (isTRUE(verbose)) message(...)
  
  # ---- helpers --------------------------------------------------------------
  .add_ids <- function(df) {
    rn <- rownames(df)
    ensembl_id  <- sub("^([^_-]+)[_-].*$", "\\1", rn)
    gene_symbol <- ifelse(grepl("^[^_-]+[_-]", rn),
                          sub("^[^_-]+[_-](.*)$", "\\1", rn),
                          NA_character_)
    cbind(
      data.frame(
        name_id     = rn,
        ensembl_id  = ensembl_id,
        gene_symbol = gene_symbol,
        stringsAsFactors = FALSE
      ),
      df,
      row.names = NULL
    )
  }
  
  .sum_by_rownames <- function(A, B) {
    all_genes <- union(rownames(A), rownames(B))
    common_samples <- intersect(colnames(A), colnames(B))
    if (length(common_samples) == 0L)
      stop("No overlapping sample columns between matrices.")
    A2 <- matrix(0, nrow = length(all_genes), ncol = length(common_samples),
                 dimnames = list(all_genes, common_samples))
    B2 <- A2
    A2[rownames(A), intersect(colnames(A), common_samples)] <- 
      A[, intersect(colnames(A), common_samples), drop = FALSE]
    B2[rownames(B), intersect(colnames(B), common_samples)] <- 
      B[, intersect(colnames(B), common_samples), drop = FALSE]
    A2 + B2
  }
  
  .as_total_counts <- function(agg) {
    if (is.list(agg)) {
      mats <- Filter(function(m) (is.matrix(m) || inherits(m, "Matrix")) &&
                       nrow(m) > 0 && ncol(m) > 0, agg)
      if (!length(mats)) stop("Input list is empty or contains empty matrices.")
      .msg(sprintf("[edgeR] Merging %d cluster matrix/matrices into one global matrix...", length(mats)))
      mats <- lapply(mats, as.matrix)
      Reduce(.sum_by_rownames, mats)
    } else if (is.matrix(agg) || inherits(agg, "Matrix")) {
      agg <- as.matrix(agg)
      if (nrow(agg) == 0L || ncol(agg) == 0L) stop("Input matrix is empty.")
      .msg("[edgeR] Using provided matrix (genes x samples).")
      agg
    } else {
      stop("'agg' must be a matrix, sparse Matrix, or a list of such objects.")
    }
  }
  
  # ---- metadata prep --------------------------------------------------------
  stopifnot(all(c(sample_col, group_col) %in% colnames(sample_meta)))
  sample_meta <- as.data.frame(sample_meta)
  if (any(duplicated(sample_meta[[sample_col]]))) {
    .msg("[edgeR] Deduplicating sample metadata by '", sample_col, "'.")
    sample_meta <- sample_meta[!duplicated(sample_meta[[sample_col]]), , drop = FALSE]
  }
  if (!is.null(group_levels)) {
    sample_meta[[group_col]] <- factor(sample_meta[[group_col]], levels = group_levels)
  } else {
    sample_meta[[group_col]] <- factor(sample_meta[[group_col]])
  }
  
  # ---- counts assembly ------------------------------------------------------
  counts <- .as_total_counts(agg)   # genes x samples
  smp_ids <- colnames(counts)
  smeta   <- sample_meta[match(smp_ids, sample_meta[[sample_col]]), , drop = FALSE]
  if (anyNA(smeta[[group_col]])) {
    miss <- smp_ids[is.na(smeta[[group_col]])]
    stop("Missing group labels for samples: ", paste(miss, collapse = ", "))
  }
  counts <- counts[, smeta[[sample_col]], drop = FALSE]
  
  .msg(sprintf("[edgeR] Global pseudobulk: %d genes × %d samples.", nrow(counts), ncol(counts)))
  grp_tab <- table(smeta[[group_col]])
  .msg(sprintf("[edgeR] Samples per group: %s",
               paste(sprintf("%s=%d", names(grp_tab), as.integer(grp_tab)), collapse = ", ")))
  
  # ---- edgeR core: filtering & normalization --------------------------------
  .msg("[edgeR] Building DGEList, filtering (filterByExpr), and normalizing (TMM)...")
  dge <- edgeR::DGEList(counts = counts, samples = smeta)
  keep <- edgeR::filterByExpr(dge, group = dge$samples[[group_col]])
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  dge  <- edgeR::calcNormFactors(dge)
  .msg(sprintf("[edgeR] Kept %d/%d genes after filtering.", nrow(dge), nrow(counts)))
  
  # ---- global ANOVA-like test (QL F-test) -----------------------------------
  .msg("[edgeR] Fitting global model (~ group) and running QL F-test (ANOVA-like)...")
  design_g <- stats::model.matrix(stats::as.formula(paste0("~ ", group_col)), data = dge$samples)
  dge <- edgeR::estimateDisp(dge, design_g)
  fit_g <- edgeR::glmQLFit(dge, design_g, robust = robust)
  qlf_global <- edgeR::glmQLFTest(fit_g, coef = 2:ncol(design_g))
  tab_global <- edgeR::topTags(qlf_global, n = Inf)$table
  tab_global <- .add_ids(tab_global)
  sig_n <- sum(tab_global$FDR <= FDR_global_threshold, na.rm = TRUE)
  .msg(sprintf("[edgeR] Global significant (FDR ≤ %.3f): %d gene(s).", FDR_global_threshold, sig_n))
  
  # ---- post-hoc pairwise contrasts ------------------------------------------
  .msg("[edgeR] Post-hoc contrasts: DEX-CTRL, CORT-CTRL, CORT-DEX...")
  design0 <- stats::model.matrix(stats::as.formula(paste0("~ 0 + ", group_col)), data = dge$samples)
  colnames(design0) <- levels(dge$samples[[group_col]])
  fit0 <- edgeR::glmQLFit(dge, design0, robust = robust)
  make_contrast <- function(a, b) limma::makeContrasts(contrasts = paste0(a, "-", b), levels = design0)
  
  q_DEX_CTRL  <- edgeR::glmQLFTest(fit0, contrast = make_contrast("DEX","CTRL"))
  q_CORT_CTRL <- edgeR::glmQLFTest(fit0, contrast = make_contrast("CORT","CTRL"))
  q_CORT_DEX  <- edgeR::glmQLFTest(fit0, contrast = make_contrast("CORT","DEX"))
  
  tab_DEX_CTRL  <- .add_ids(edgeR::topTags(q_DEX_CTRL,  n = Inf)$table)
  tab_CORT_CTRL <- .add_ids(edgeR::topTags(q_CORT_CTRL, n = Inf)$table)
  tab_CORT_DEX  <- .add_ids(edgeR::topTags(q_CORT_DEX,  n = Inf)$table)
  .msg("[edgeR] Post-hoc contrasts completed.")
  
  # ---- one combined table ---------------------------------------------------
  .pick <- function(df, cols) df[, intersect(cols, colnames(df)), drop = FALSE]
  id_cols <- c("name_id","ensembl_id","gene_symbol")
  
  g <- .pick(tab_global, c(id_cols, "logCPM","F","PValue","FDR"))
  colnames(g) <- sub("^logCPM$","g_logCPM", colnames(g))
  colnames(g) <- sub("^F$","g_F", colnames(g))
  colnames(g) <- sub("^PValue$","g_PValue", colnames(g))
  colnames(g) <- sub("^FDR$","g_FDR", colnames(g))
  
  .with_prefix <- function(df, prefix){
    x <- .pick(df, c(id_cols, "logFC","logCPM","F","PValue","FDR"))
    nm <- colnames(x)
    nm <- ifelse(nm %in% id_cols, nm, paste0(prefix, nm))
    colnames(x) <- nm
    x
  }
  d  <- .with_prefix(tab_DEX_CTRL,  "DEX_CTRL_")
  c1 <- .with_prefix(tab_CORT_CTRL, "CORT_CTRL_")
  cd <- .with_prefix(tab_CORT_DEX,  "CORT_DEX_")
  
  combined <- Reduce(function(a,b) merge(a, b, by = id_cols, all = TRUE), list(g, d, c1, cd))
  ord <- c(
    id_cols,
    c("g_logCPM","g_F","g_PValue","g_FDR"),
    grep("^DEX_CTRL_",  colnames(combined), value = TRUE),
    grep("^CORT_CTRL_", colnames(combined), value = TRUE),
    grep("^CORT_DEX_",  colnames(combined), value = TRUE)
  )
  combined <- combined[, ord[ord %in% colnames(combined)], drop = FALSE]
  rownames(combined) <- NULL
  
  .msg("[edgeR] Done.")
  list(
    f_test_global = tab_global,
    DEX_vs_CTRL   = tab_DEX_CTRL,
    CORT_vs_CTRL  = tab_CORT_CTRL,
    CORT_vs_DEX   = tab_CORT_DEX,
    combined      = combined
  )
}
