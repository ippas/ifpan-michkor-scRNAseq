edger_pseudobulk_anova_per_cluster_v3 <- 
  function(
    agg_list,
    sample_meta,
    group_col            = "treatment",
    sample_col           = "sample_id",
    group_levels         = c("CTRL","DEX","CORT"),
    clusters             = NULL,
    robust               = TRUE,
    progress             = TRUE,
    verbose              = TRUE,
    gatekeep             = FALSE,
    FDR_global_threshold = 0.05,
    fix_dispersion_if_no_replicates = TRUE,
    fixed_dispersion_value          = 0.1
  ) {
    if (!requireNamespace("edgeR", quietly = TRUE)) stop("Package 'edgeR' is required.")
    if (gatekeep && !requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required when gatekeep = TRUE.")
    
    if (is.null(clusters)) clusters <- names(agg_list)
    clusters <- intersect(clusters, names(agg_list))
    if (!length(clusters)) stop("No clusters to analyze.")
    n_clusters <- length(clusters)
    
    stopifnot(all(c(sample_col, group_col) %in% colnames(sample_meta)))
    sample_meta <- as.data.frame(sample_meta)
    if (any(duplicated(sample_meta[[sample_col]]))) {
      sample_meta <- sample_meta[!duplicated(sample_meta[[sample_col]]), , drop = FALSE]
    }
    if (!is.null(group_levels)) {
      sample_meta[[group_col]] <- factor(sample_meta[[group_col]], levels = group_levels)
    } else {
      sample_meta[[group_col]] <- factor(sample_meta[[group_col]])
    }
    
    .add_ids <- function(df) {
      rn <- rownames(df)
      ensembl_id  <- sub("^([^_-]+)[_-].*$", "\\1", rn)
      gene_symbol <- ifelse(grepl("^[^_-]+[_-]", rn),
                            sub("^[^_-]+[_-](.*)$", "\\1", rn),
                            NA_character_)
      out <- cbind(
        data.frame(
          name_id    = rn,
          ensembl_id = ensembl_id,
          gene_symbol = gene_symbol,
          stringsAsFactors = FALSE
        ),
        df,
        row.names = NULL
      )
      out
    }
    
    if (verbose) {
      message(sprintf(
        "[edgeR] Starting analysis for %d cluster(s); robust=%s; gatekeep=%s; fix_dispersion_if_no_replicates=%s (%.2f)",
        n_clusters, robust, gatekeep, fix_dispersion_if_no_replicates, fixed_dispersion_value
      ))
      message(sprintf("[edgeR] Groups (levels): %s", paste(levels(sample_meta[[group_col]]), collapse = ", ")))
    }
    
    t_total_start <- Sys.time()
    per_cluster_sec <- numeric(n_clusters)
    if (progress) pb <- utils::txtProgressBar(min = 0, max = n_clusters, style = 3)
    out <- setNames(vector("list", n_clusters), clusters)
    
    for (i in seq_along(clusters)) {
      cl <- clusters[i]
      t0 <- Sys.time()
      if (verbose) message(sprintf("\n[edgeR] Cluster %s (%d/%d)...", cl, i, n_clusters))
      
      mat <- agg_list[[cl]]
      if (is.null(mat) || ncol(mat) == 0L) {
        warning("Empty matrix for ", cl, "; skipping.")
        if (progress) utils::setTxtProgressBar(pb, i)
        next
      }
      
      smp_ids <- colnames(mat)
      smeta <- sample_meta[match(smp_ids, sample_meta[[sample_col]]), , drop = FALSE]
      mat <- mat[, smeta[[sample_col]], drop = FALSE]
      
      grp_tab <- table(smeta[[group_col]])
      min_reps <- min(grp_tab)
      no_replicates <- any(grp_tab < 2)
      
      if (verbose) {
        message(sprintf("[edgeR]   Samples per group: %s", 
                        paste(sprintf("%s=%d", names(grp_tab), as.integer(grp_tab)), collapse = ", ")))
        if (no_replicates) message(sprintf("[edgeR]   ⚠️ Detected group(s) with <2 replicates (min=%d)", min_reps))
      }
      
      dge <- edgeR::DGEList(counts = mat, samples = smeta)
      keep <- edgeR::filterByExpr(dge, group = dge$samples[[group_col]])
      dge <- dge[keep, , keep.lib.sizes = FALSE]
      dge <- edgeR::calcNormFactors(dge)
      
      design_g <- stats::model.matrix(stats::as.formula(paste0("~ ", group_col)),
                                      data = dge$samples)
      
      if (fix_dispersion_if_no_replicates && no_replicates) {
        # --- brak replikatów: ustaw dyspersję ręcznie ---
        if (verbose) message(sprintf("[edgeR]   Using fixed dispersion = %.3f (no biological replicates)", fixed_dispersion_value))
        dge$common.dispersion <- fixed_dispersion_value
        fit_g <- edgeR::glmFit(dge, design_g)
        qlf_global <- edgeR::glmLRT(fit_g, coef = 2:ncol(design_g))
      } else {
        # --- klasyczny przypadek ---
        dge <- edgeR::estimateDisp(dge, design_g)
        fit_g <- edgeR::glmQLFit(dge, design_g, robust = robust)
        qlf_global <- edgeR::glmQLFTest(fit_g, coef = 2:ncol(design_g))
      }
      
      tab_global <- edgeR::topTags(qlf_global, n = Inf)$table %>% .add_ids()
      
      # --- Post-hoc kontrasty ---
      design0 <- stats::model.matrix(stats::as.formula(paste0("~ 0 + ", group_col)),
                                     data = dge$samples)
      colnames(design0) <- levels(dge$samples[[group_col]])
      
      if (fix_dispersion_if_no_replicates && no_replicates) {
        fit0 <- edgeR::glmFit(dge, design0)
        make_contrast <- function(a, b) limma::makeContrasts(contrasts = paste0(a, "-", b), levels = design0)
        q_DEX_CTRL  <- edgeR::glmLRT(fit0, contrast = make_contrast("DEX","CTRL"))
        q_CORT_CTRL <- edgeR::glmLRT(fit0, contrast = make_contrast("CORT","CTRL"))
        q_CORT_DEX  <- edgeR::glmLRT(fit0, contrast = make_contrast("CORT","DEX"))
      } else {
        fit0 <- edgeR::glmQLFit(dge, design0, robust = robust)
        make_contrast <- function(a, b) limma::makeContrasts(contrasts = paste0(a, "-", b), levels = design0)
        q_DEX_CTRL  <- edgeR::glmQLFTest(fit0, contrast = make_contrast("DEX","CTRL"))
        q_CORT_CTRL <- edgeR::glmQLFTest(fit0, contrast = make_contrast("CORT","CTRL"))
        q_CORT_DEX  <- edgeR::glmQLFTest(fit0, contrast = make_contrast("CORT","DEX"))
      }
      
      tab_DEX_CTRL  <- .add_ids(edgeR::topTags(q_DEX_CTRL,  n = Inf)$table)
      tab_CORT_CTRL <- .add_ids(edgeR::topTags(q_CORT_CTRL, n = Inf)$table)
      tab_CORT_DEX  <- .add_ids(edgeR::topTags(q_CORT_DEX,  n = Inf)$table)
      
      if (gatekeep) {
        sig_ids <- subset(tab_global, FDR <= FDR_global_threshold)[, c("ensembl_id","gene_symbol")]
        if (nrow(sig_ids)) {
          tab_DEX_CTRL  <- dplyr::semi_join(tab_DEX_CTRL,  sig_ids, by = c("ensembl_id","gene_symbol"))
          tab_CORT_CTRL <- dplyr::semi_join(tab_CORT_CTRL, sig_ids, by = c("ensembl_id","gene_symbol"))
          tab_CORT_DEX  <- dplyr::semi_join(tab_CORT_DEX,  sig_ids, by = c("ensembl_id","gene_symbol"))
        } else {
          tab_DEX_CTRL  <- tab_CORT_CTRL <- tab_CORT_DEX <- tab_global[0, , drop = FALSE]
        }
      }
      
      # --- combined summary ---
      .pick <- function(df, cols) df[, intersect(cols, colnames(df)), drop = FALSE]
      id_cols <- c("name_id","ensembl_id","gene_symbol")
      g  <- .pick(tab_global, c(id_cols, "logCPM","F","PValue","FDR"))
      colnames(g) <- sub("^logCPM$","g_logCPM", colnames(g))
      colnames(g) <- sub("^F$","g_F", colnames(g))
      colnames(g) <- sub("^PValue$","g_PValue", colnames(g))
      colnames(g) <- sub("^FDR$","g_FDR", colnames(g))
      
      .with_prefix <- function(df, prefix){
        x <- .pick(df, c(id_cols, "logFC","logCPM","F","PValue","FDR"))
        colnames(x) <- ifelse(colnames(x) %in% id_cols, colnames(x), paste0(prefix, colnames(x)))
        x
      }
      d  <- .with_prefix(tab_DEX_CTRL,  "DEX_CTRL_")
      c1 <- .with_prefix(tab_CORT_CTRL, "CORT_CTRL_")
      cd <- .with_prefix(tab_CORT_DEX,  "CORT_DEX_")
      comb <- Reduce(function(a,b) merge(a,b,by=id_cols,all=TRUE), list(g,d,c1,cd))
      
      out[[cl]] <- list(
        f_test_global = tab_global,
        DEX_vs_CTRL   = tab_DEX_CTRL,
        CORT_vs_CTRL  = tab_CORT_CTRL,
        CORT_vs_DEX   = tab_CORT_DEX,
        combined      = comb
      )
      
      per_cluster_sec[i] <- as.numeric(difftime(Sys.time(), t0, units="secs"))
      if (progress) utils::setTxtProgressBar(pb, i)
    }
    if (progress) close(pb)
    
    combined_all <- do.call(rbind, lapply(names(out), function(cl){
      if (is.null(out[[cl]]$combined)) return(NULL)
      cbind(cluster = cl, out[[cl]]$combined)
    }))
    
    attr(out, "combined_all") <- combined_all
    invisible(out)
  }
