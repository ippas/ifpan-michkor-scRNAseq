edger_pseudobulk_anova_per_cluster_v2 <- 
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
    FDR_global_threshold = 0.05
  ) {
    if (!requireNamespace("edgeR", quietly = TRUE)) stop("Package 'edgeR' is required.")
    if (gatekeep && !requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required when gatekeep = TRUE.")
    
    # Choose clusters
    if (is.null(clusters)) clusters <- names(agg_list)
    clusters <- intersect(clusters, names(agg_list))
    if (!length(clusters)) stop("No clusters to analyze.")
    n_clusters <- length(clusters)
    
    # Meta prep
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
    
    # Keep original rownames as name_id; add parsed IDs
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
        "[edgeR] Starting analysis for %d cluster(s); gatekeep=%s (FDR_global_threshold=%.3f); robust=%s",
        n_clusters, gatekeep, FDR_global_threshold, robust
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
      
      if (verbose) message(sprintf("\n[edgeR] Cluster %s (%d/%d) — preparing data...", cl, i, n_clusters))
      
      mat <- agg_list[[cl]]  # genes x samples
      if (is.null(mat) || ncol(mat) == 0L) {
        warning("Empty matrix for ", cl, "; skipping.")
        if (progress) utils::setTxtProgressBar(pb, i)
        per_cluster_sec[i] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        next
      }
      
      # Align sample metadata to columns
      smp_ids <- colnames(mat)
      smeta <- sample_meta[match(smp_ids, sample_meta[[sample_col]]), , drop = FALSE]
      if (anyNA(smeta[[group_col]])) {
        miss <- smp_ids[is.na(smeta[[group_col]])]
        stop("Cluster ", cl, ": missing group labels for samples: ", paste(miss, collapse=", "))
      }
      mat <- mat[, smeta[[sample_col]], drop = FALSE]
      
      if (verbose) {
        grp_tab <- table(smeta[[group_col]])
        message(sprintf("[edgeR]   Samples per group: %s",
                        paste(sprintf("%s=%d", names(grp_tab), as.integer(grp_tab)), collapse = ", ")))
        message("[edgeR]   Building DGEList / filtering / normalization...")
      }
      
      # edgeR core
      dge <- edgeR::DGEList(counts = mat, samples = smeta)
      keep <- edgeR::filterByExpr(dge, group = dge$samples[[group_col]])
      dge <- dge[keep, , keep.lib.sizes = FALSE]
      dge <- edgeR::calcNormFactors(dge)
      
      if (verbose) message("[edgeR]   Global ANOVA-like (QL F-test)...")
      design_g <- stats::model.matrix(stats::as.formula(paste0("~ ", group_col)),
                                      data = dge$samples)
      dge <- edgeR::estimateDisp(dge, design_g)
      fit_g <- edgeR::glmQLFit(dge, design_g, robust = robust)
      qlf_global <- edgeR::glmQLFTest(fit_g, coef = 2:ncol(design_g))
      tab_global <- edgeR::topTags(qlf_global, n = Inf)$table
      tab_global <- .add_ids(tab_global)
      
      if (verbose) {
        sig_n <- sum(tab_global$FDR <= FDR_global_threshold, na.rm = TRUE)
        message(sprintf("[edgeR]   Global significant (FDR <= %.3f): %d gene(s)", FDR_global_threshold, sig_n))
        message("[edgeR]   Post-hoc contrasts (DEX-CTRL, CORT-CTRL, CORT-DEX)...")
      }
      
      # Post-hoc contrasts on no-intercept design
      design0 <- stats::model.matrix(stats::as.formula(paste0("~ 0 + ", group_col)),
                                     data = dge$samples)
      colnames(design0) <- levels(dge$samples[[group_col]])
      fit0 <- edgeR::glmQLFit(dge, design0, robust = robust)
      make_contrast <- function(a, b) limma::makeContrasts(contrasts = paste0(a, "-", b), levels = design0)
      
      q_DEX_CTRL  <- edgeR::glmQLFTest(fit0, contrast = make_contrast("DEX","CTRL"))
      q_CORT_CTRL <- edgeR::glmQLFTest(fit0, contrast = make_contrast("CORT","CTRL"))
      q_CORT_DEX  <- edgeR::glmQLFTest(fit0, contrast = make_contrast("CORT","DEX"))
      
      tab_DEX_CTRL  <- .add_ids(edgeR::topTags(q_DEX_CTRL,  n = Inf)$table)
      tab_CORT_CTRL <- .add_ids(edgeR::topTags(q_CORT_CTRL, n = Inf)$table)
      tab_CORT_DEX  <- .add_ids(edgeR::topTags(q_CORT_DEX,  n = Inf)$table)
      
      # Optional gatekeeping
      if (gatekeep) {
        sig_ids <- subset(tab_global, FDR <= FDR_global_threshold)[, c("ensembl_id","gene_symbol")]
        if (nrow(sig_ids)) {
          tab_DEX_CTRL  <- dplyr::semi_join(tab_DEX_CTRL,  sig_ids, by = c("ensembl_id","gene_symbol"))
          tab_CORT_CTRL <- dplyr::semi_join(tab_CORT_CTRL, sig_ids, by = c("ensembl_id","gene_symbol"))
          tab_CORT_DEX  <- dplyr::semi_join(tab_CORT_DEX,  sig_ids, by = c("ensembl_id","gene_symbol"))
        } else {
          tab_DEX_CTRL  <- tab_DEX_CTRL[0, , drop = FALSE]
          tab_CORT_CTRL <- tab_CORT_CTRL[0, , drop = FALSE]
          tab_CORT_DEX  <- tab_CORT_DEX[0, , drop = FALSE]
        }
        if (verbose) message(sprintf("[edgeR]   Gatekeeping applied (FDR <= %.3f): post-hoc restricted to global hits.", FDR_global_threshold))
      }
      
      # ---- NEW: build one combined table --------------------------------------
      .pick <- function(df, cols) df[, intersect(cols, colnames(df)), drop = FALSE]
      id_cols <- c("name_id","ensembl_id","gene_symbol")
      
      # Global (no logFC here)
      g <- .pick(tab_global, c(id_cols, "logCPM","F","PValue","FDR"))
      colnames(g) <- sub("^logCPM$","g_logCPM", colnames(g))
      colnames(g) <- sub("^F$","g_F", colnames(g))
      colnames(g) <- sub("^PValue$","g_PValue", colnames(g))
      colnames(g) <- sub("^FDR$","g_FDR", colnames(g))
      
      # Helper to prefix columns for pairwise tables
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
      
      # full outer joins by IDs
      comb <- Reduce(function(a,b) merge(a, b, by = id_cols, all = TRUE), list(g, d, c1, cd))
      
      # Porządek kolumn: ID -> global -> DEX_CTRL -> CORT_CTRL -> CORT_DEX
      ord <- c(
        id_cols,
        c("g_logCPM","g_F","g_PValue","g_FDR"),
        grep("^DEX_CTRL_",  colnames(comb), value = TRUE),
        grep("^CORT_CTRL_", colnames(comb), value = TRUE),
        grep("^CORT_DEX_",  colnames(comb), value = TRUE)
      )
      comb <- comb[, ord[ord %in% colnames(comb)], drop = FALSE]
      rownames(comb) <- NULL
      # -------------------------------------------------------------------------
      
      out[[cl]] <- list(
        f_test_global = tab_global,
        DEX_vs_CTRL   = tab_DEX_CTRL,
        CORT_vs_CTRL  = tab_CORT_CTRL,
        CORT_vs_DEX   = tab_CORT_DEX,
        combined      = comb
      )
      
      dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      per_cluster_sec[i] <- dt
      if (verbose) {
        remaining <- n_clusters - i
        mean_so_far <- mean(per_cluster_sec[per_cluster_sec > 0], na.rm = TRUE)
        eta_sec <- if (remaining > 0 && is.finite(mean_so_far)) remaining * mean_so_far else 0
        message(sprintf("[edgeR] Cluster %s done in %.1fs. Remaining: %d. ETA ~ %.1fs",
                        cl, dt, remaining, eta_sec))
      }
      if (progress) utils::setTxtProgressBar(pb, i)
    }
    if (progress) close(pb)
    
    t_total <- as.numeric(difftime(Sys.time(), t_total_start, units = "secs"))
    per_ok <- per_cluster_sec[is.finite(per_cluster_sec) & per_cluster_sec > 0]
    if (verbose) {
      if (length(per_ok)) {
        message("\n[edgeR] Timing summary:")
        message(sprintf("  Total: %.1fs for %d cluster(s)", t_total, n_clusters))
        message(sprintf("  Per cluster (s): mean=%.1f, median=%.1f, min=%.1f, max=%.1f",
                        mean(per_ok), stats::median(per_ok), min(per_ok), max(per_ok)))
      } else {
        message(sprintf("\n[edgeR] Total time: %.1fs", t_total))
      }
    }
    
    # --- NEW: one big table across clusters
    combined_all <- do.call(rbind, lapply(names(out), function(cl){
      if (is.null(out[[cl]]$combined)) return(NULL)
      cbind(cluster = cl, out[[cl]]$combined)
    }))
    
    attr(out, "combined_all") <- combined_all
    return(out)
  }


# edger_pseudobulk_anova_per_cluster 
