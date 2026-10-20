run_fgsea_on_score_matrix <- function(
    mat,
    gene_lists,
    nperm   = 1000,
    minSize = 5,
    maxSize = Inf,
    seed    = 1,
    verbose = TRUE,
    do_fisher           = FALSE,
    fisher_universe     = NULL,
    fisher_min_size     = 1,
    fisher_alternative  = "greater",
    fisher_strict_checks= TRUE,
    fisher_smooth_or    = FALSE,
    rename_pathway_to   = "predict_cell_type",   # set NULL to keep 'pathway'
    rename_size_to      = "n_genes_in_set",      # e.g. "n_genes_overlap"
    add_leadingEdge_n   = TRUE,                  # add 'n_leadingEdge' column
    export_xlsx         = FALSE,
    xlsx_path           = NULL,
    write_first_by_NES  = TRUE,
    write_first_by_pval = TRUE,
    write_first_by_padj = TRUE,
    digits_xlsx         = 3,                     # visible decimals for float cols
    xlsx_float_cols     = c("pval","padj","ES","NES"),
    keep_raw_results    = TRUE,                  # keep original fgsea data.frames
    store_rank_stats    = TRUE,                  # keep stats vectors for plotting
    export_raw_in_xlsx  = FALSE                  # also write raw_* sheets (leadingEdge collapsed)
) {
  # --- 0) Dependencies & checks ---
  stopifnot(requireNamespace("fgsea", quietly = TRUE))
  stopifnot(requireNamespace("dplyr", quietly = TRUE))
  stopifnot(requireNamespace("purrr", quietly = TRUE))
  stopifnot(requireNamespace("tibble", quietly = TRUE))
  if (isTRUE(export_xlsx)) stopifnot(requireNamespace("openxlsx", quietly = TRUE))
  
  if (!is.matrix(mat) && !is.data.frame(mat)) stop("`mat` must be a matrix or data.frame (rows = genes, cols = clusters).")
  mat <- as.data.frame(mat, check.names = FALSE)
  if (is.null(rownames(mat)) || any(rownames(mat) == "")) stop("`mat` must have rownames set to gene symbols (HGNC).")
  if (!is.list(gene_lists) || !length(gene_lists)) stop("`gene_lists` must be a non-empty named list of character vectors.")
  set.seed(seed)
  
  # --- 1) Clean matrix ---
  mat[] <- lapply(mat, function(x) as.numeric(x))
  keep_rows <- apply(as.matrix(mat), 1, function(x) any(is.finite(x)))
  mat <- mat[keep_rows, , drop = FALSE]
  if (verbose) message("Kept ", nrow(mat), " genes after removing all-non-finite rows.")
  
  # --- 2) Intersect gene lists ---
  all_symbols <- rownames(mat)
  cleaned_lists <- purrr::imap(gene_lists, function(v, nm) {
    v <- unique(as.character(v)); v <- v[!is.na(v) & nzchar(v)]
    intersect(v, all_symbols)
  })
  n_before <- length(cleaned_lists)
  cleaned_lists <- cleaned_lists[lengths(cleaned_lists) >= minSize]
  if (verbose) message("Pathway sets: ", n_before, " → ", length(cleaned_lists), " after intersection & minSize ≥ ", minSize, ".")
  if (!length(cleaned_lists)) stop("No gene set remains after intersection and minSize filtering.")
  
  # helpers
  .maybe_rename_pathway <- function(df) {
    if (!is.null(rename_pathway_to) && "pathway" %in% names(df)) {
      names(df)[names(df) == "pathway"] <- rename_pathway_to
    }
    df
  }
  path_col <- if (is.null(rename_pathway_to)) "pathway" else rename_pathway_to
  
  # --- 3) FGSEA for each cluster ---
  if (verbose) message("Running fgsea for ", ncol(mat), " clusters …")
  cluster_order <- colnames(mat)
  
  res_list  <- vector("list", length(cluster_order)); names(res_list)  <- cluster_order
  raw_list  <- if (keep_raw_results)  vector("list", length(cluster_order)) else NULL
  stats_map <- if (store_rank_stats)  vector("list", length(cluster_order)) else NULL
  names(raw_list) <- names(stats_map) <- cluster_order
  
  for (cl in cluster_order) {
    v <- as.numeric(mat[[cl]])
    ok <- is.finite(v)
    stats <- v[ok]; names(stats) <- rownames(mat)[ok]
    stats <- sort(stats, decreasing = TRUE)
    
    out_raw <- suppressWarnings(
      fgsea::fgsea(pathways = cleaned_lists, stats = stats,
                   nperm = nperm, minSize = minSize, maxSize = maxSize)
    )
    if (keep_raw_results)  raw_list[[cl]]  <- out_raw
    if (store_rank_stats)  stats_map[[cl]] <- stats
    
    if (!nrow(out_raw)) {
      df <- tibble::tibble(
        pathway = character(), pval = numeric(), padj = numeric(),
        ES = numeric(), NES = numeric(), nMoreExtreme = integer(),
        size = integer(), cluster = character(),
        genes = character()
      )
    } else {
      out_raw$cluster <- cl
      if (add_leadingEdge_n && "leadingEdge" %in% names(out_raw)) {
        out_raw$n_leadingEdge <- lengths(out_raw$leadingEdge)
      }
      df <- out_raw %>%
        dplyr::mutate(genes = if ("leadingEdge" %in% names(.)) purrr::map_chr(leadingEdge, ~ paste(.x, collapse = "|")) else "") %>%
        dplyr::select(-tidyselect::any_of("leadingEdge")) %>%
        dplyr::arrange(dplyr::desc(NES))
    }
    # rename columns for processed result
    if ("size" %in% names(df) && !is.null(rename_size_to) && rename_size_to != "size") {
      names(df)[names(df) == "size"] <- rename_size_to
    }
    df <- .maybe_rename_pathway(df)
    res_list[[cl]] <- df
  }
  
  # --- 4) Summaries (ties allowed), ordered by cluster ---
  first_by_NES <- res_list %>% purrr::imap(~ if (nrow(.x)) dplyr::slice_max(.x, .x$NES, n = 1, with_ties = TRUE) else NULL) %>%
    purrr::discard(is.null) %>% dplyr::bind_rows()
  if (nrow(first_by_NES)) first_by_NES <- first_by_NES %>%
    dplyr::mutate(cluster = factor(cluster, levels = cluster_order)) %>%
    dplyr::arrange(cluster, dplyr::desc(NES)) %>%
    dplyr::mutate(cluster = as.character(cluster))
  
  first_by_pval <- res_list %>% purrr::imap(~ if (nrow(.x)) dplyr::slice_min(.x, .x$pval, n = 1, with_ties = TRUE) else NULL) %>%
    purrr::discard(is.null) %>% dplyr::bind_rows()
  if (nrow(first_by_pval)) first_by_pval <- first_by_pval %>%
    dplyr::mutate(cluster = factor(cluster, levels = cluster_order)) %>%
    dplyr::arrange(cluster, pval, dplyr::desc(NES)) %>%
    dplyr::mutate(cluster = as.character(cluster))
  
  first_by_padj <- res_list %>% purrr::imap(~ if (nrow(.x)) dplyr::slice_min(.x, .x$padj, n = 1, with_ties = TRUE) else NULL) %>%
    purrr::discard(is.null) %>% dplyr::bind_rows()
  if (nrow(first_by_padj)) first_by_padj <- first_by_padj %>%
    dplyr::mutate(cluster = factor(cluster, levels = cluster_order)) %>%
    dplyr::arrange(cluster, padj, dplyr::desc(NES)) %>%
    dplyr::mutate(cluster = as.character(cluster))
  
  # --- 5) Optional: Fisher tests ---
  fisher_out <- NULL
  if (isTRUE(do_fisher)) {
    if (!exists("compute_geneSet_fisher")) stop("`compute_geneSet_fisher()` not found.")
    if (is.null(fisher_universe)) fisher_universe <- rownames(mat)
    if (verbose) message("Running Fisher tests on marker sets …")
    fisher_out <- compute_geneSet_fisher(
      gene_lists = cleaned_lists, universe = fisher_universe,
      min_size = fisher_min_size, alternative = fisher_alternative,
      strict_checks = fisher_strict_checks, smooth_or = fisher_smooth_or
    )
  }
  
  # --- 6) Optional XLSX export (openxlsx) ---
  if (isTRUE(export_xlsx)) {
    sanitize_sheet_name <- function(nm) {
      nm <- gsub("[\\*/\\?\\[\\]:]", "_", nm); if (nchar(nm) > 31) substr(nm, 1, 31) else nm
    }
    if (is.null(xlsx_path)) xlsx_path <- file.path(getwd(), paste0("fgsea_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx"))
    
    wb <- openxlsx::createWorkbook()
    
    add_sheet_df <- function(wb, sheet_name, df, digits = 3, float_cols = xlsx_float_cols) {
      sheet_name <- sanitize_sheet_name(sheet_name)
      openxlsx::addWorksheet(wb, sheet_name)
      openxlsx::writeData(wb, sheet = sheet_name, x = as.data.frame(df), keepNA = TRUE)
      openxlsx::freezePane(wb, sheet = sheet_name, firstRow = TRUE)
      openxlsx::addFilter(wb, sheet = sheet_name, rows = 1, cols = 1:ncol(df))
      openxlsx::setColWidths(wb, sheet = sheet_name, cols = 1:ncol(df), widths = "auto")
      
      # format ONLY selected float columns (pval, padj, ES, NES by default)
      present <- intersect(float_cols, colnames(df))
      if (length(present) > 0 && nrow(df) > 0) {
        num_fmt <- if (digits <= 0) "0" else paste0("0.", paste(rep("0", digits), collapse = ""))
        num_style <- openxlsx::createStyle(numFmt = num_fmt)
        float_idx <- match(present, colnames(df))
        openxlsx::addStyle(wb, sheet = sheet_name, style = num_style,
                           cols = float_idx, rows = 2:(nrow(df) + 1), gridExpand = TRUE)
      }
    }
    
    # summaries first
    if (write_first_by_NES  && nrow(first_by_NES))  add_sheet_df(wb, "first_by_NES",  first_by_NES,  digits = digits_xlsx)
    if (write_first_by_pval && nrow(first_by_pval)) add_sheet_df(wb, "first_by_pval", first_by_pval, digits = digits_xlsx)
    if (write_first_by_padj && nrow(first_by_padj)) add_sheet_df(wb, "first_by_padj", first_by_padj, digits = digits_xlsx)
    
    # per-cluster processed sheets
    for (cl in cluster_order) {
      add_sheet_df(wb, cl, res_list[[cl]], digits = digits_xlsx)
    }
    
    # optional: raw_* sheets with leadingEdge collapsed (Excel can't store list-cols)
    if (export_raw_in_xlsx && keep_raw_results) {
      for (cl in cluster_order) {
        raw_df <- raw_list[[cl]]
        if (!is.null(raw_df) && nrow(raw_df)) {
          raw_df$leadingEdge_collapsed <- vapply(raw_df$leadingEdge, function(le) paste(le, collapse = "|"), "")
          raw_df <- raw_df %>% dplyr::select(-leadingEdge)
          add_sheet_df(wb, paste0("raw_", cl), raw_df, digits = digits_xlsx)
        } else {
          add_sheet_df(wb, paste0("raw_", cl), data.frame(), digits = digits_xlsx)
        }
      }
    }
    
    if (verbose) message("Writing XLSX to: ", xlsx_path)
    openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)
  }
  
  # --- 7) Return ---
  out <- list(
    results           = res_list,        # processed per-cluster (renamed cols, genes collapsed)
    raw_results       = if (keep_raw_results) raw_list else NULL,   # original fgsea per cluster
    rank_stats        = if (store_rank_stats) stats_map else NULL,  # named vectors for plotting
    first_by_NES      = first_by_NES,
    first_by_pval     = first_by_pval,
    first_by_padj     = first_by_padj,
    gene_lists_used   = cleaned_lists
  )
  if (!is.null(fisher_out)) out$fisher <- fisher_out
  if (isTRUE(export_xlsx))  out$xlsx_path <- xlsx_path
  out
}

run_fgsea_on_score_matrix <- function(
    mat,
    gene_lists,
    nperm   = 1000,
    minSize = 5,
    maxSize = Inf,
    seed    = 1,
    verbose = TRUE,
    do_fisher           = FALSE,
    fisher_universe     = NULL,
    fisher_min_size     = 1,
    fisher_alternative  = "greater",
    fisher_strict_checks= TRUE,
    fisher_smooth_or    = FALSE,
    rename_pathway_to   = "predict_cell_type",   # set NULL to keep 'pathway'
    rename_size_to      = "n_genes_in_set",
    add_leadingEdge_n   = TRUE,
    export_xlsx         = FALSE,
    xlsx_path           = NULL,
    write_first_by_NES  = TRUE,
    write_first_by_pval = TRUE,
    write_first_by_padj = TRUE,
    digits_xlsx         = 3,
    xlsx_float_cols     = c("pval","padj","ES","NES"),
    keep_raw_results    = TRUE,
    store_rank_stats    = TRUE,
    export_raw_in_xlsx  = FALSE
) {
  # --- 0) Dependencies ---
  stopifnot(requireNamespace("fgsea", quietly = TRUE))
  stopifnot(requireNamespace("dplyr", quietly = TRUE))
  stopifnot(requireNamespace("purrr", quietly = TRUE))
  stopifnot(requireNamespace("tibble", quietly = TRUE))
  if (isTRUE(export_xlsx)) stopifnot(requireNamespace("openxlsx", quietly = TRUE))
  
  if (!is.matrix(mat) && !is.data.frame(mat))
    stop("`mat` must be a matrix or data.frame (rows = genes, cols = clusters).")
  mat <- as.data.frame(mat, check.names = FALSE)
  if (is.null(rownames(mat)) || any(rownames(mat) == ""))
    stop("`mat` must have rownames set to gene symbols (HGNC).")
  if (!is.list(gene_lists) || !length(gene_lists))
    stop("`gene_lists` must be a non-empty named list of character vectors.")
  
  set.seed(seed)
  
  # --- 1) Clean matrix ---
  mat[] <- lapply(mat, function(x) as.numeric(x))
  keep_rows <- apply(as.matrix(mat), 1, function(x) any(is.finite(x)))
  mat <- mat[keep_rows, , drop = FALSE]
  if (verbose) message("Kept ", nrow(mat), " genes after removing all-non-finite rows.")
  
  # --- 2) Clean gene lists (no intersect here!) ---
  cleaned_lists <- purrr::imap(gene_lists, function(v, nm) {
    v <- unique(as.character(v))
    v[!is.na(v) & nzchar(v)]
  })
  n_before <- length(cleaned_lists)
  cleaned_lists <- cleaned_lists[lengths(cleaned_lists) >= minSize]
  if (verbose) message("Pathway sets: ", n_before, " → ", length(cleaned_lists), " after minSize ≥ ", minSize, ".")
  if (!length(cleaned_lists)) stop("No gene set remains after minSize filtering.")
  
  # helper rename
  .maybe_rename_pathway <- function(df) {
    if (!is.null(rename_pathway_to) && "pathway" %in% names(df)) {
      names(df)[names(df) == "pathway"] <- rename_pathway_to
    }
    df
  }
  path_col <- if (is.null(rename_pathway_to)) "pathway" else rename_pathway_to
  
  # --- 3) FGSEA per cluster ---
  if (verbose) message("Running fgsea for ", ncol(mat), " clusters …")
  cluster_order <- colnames(mat)
  
  res_list  <- vector("list", length(cluster_order)); names(res_list)  <- cluster_order
  raw_list  <- if (keep_raw_results)  vector("list", length(cluster_order)) else NULL
  stats_map <- if (store_rank_stats)  vector("list", length(cluster_order)) else NULL
  names(raw_list) <- names(stats_map) <- cluster_order
  
  # gene_counts collector
  gene_counts <- list()
  
  for (cl in cluster_order) {
    v <- as.numeric(mat[[cl]])
    ok <- is.finite(v)
    stats <- v[ok]; names(stats) <- rownames(mat)[ok]
    stats <- sort(stats, decreasing = TRUE)
    
    out_raw <- suppressWarnings(
      fgsea::fgsea(pathways = cleaned_lists, stats = stats,
                   nperm = nperm, minSize = minSize, maxSize = maxSize)
    )
    if (keep_raw_results)  raw_list[[cl]]  <- out_raw
    if (store_rank_stats)  stats_map[[cl]] <- stats
    
    if (!nrow(out_raw)) {
      df <- tibble::tibble(
        pathway = character(), pval = numeric(), padj = numeric(),
        ES = numeric(), NES = numeric(), nMoreExtreme = integer(),
        size = integer(), cluster = character(),
        genes = character()
      )
    } else {
      out_raw$cluster <- cl
      if (add_leadingEdge_n && "leadingEdge" %in% names(out_raw)) {
        out_raw$n_leadingEdge <- lengths(out_raw$leadingEdge)
      }
      df <- out_raw %>%
        dplyr::mutate(genes = if ("leadingEdge" %in% names(.)) purrr::map_chr(leadingEdge, ~ paste(.x, collapse = "|")) else "") %>%
        dplyr::select(-tidyselect::any_of("leadingEdge")) %>%
        dplyr::arrange(dplyr::desc(NES))
    }
    # rename cols
    if ("size" %in% names(df) && !is.null(rename_size_to) && rename_size_to != "size") {
      names(df)[names(df) == "size"] <- rename_size_to
    }
    df <- .maybe_rename_pathway(df)
    res_list[[cl]] <- df
    
    # collect gene_counts per cluster
    gene_counts[[cl]] <- purrr::imap_dfr(gene_lists, function(v, nm) {
      v_clean <- unique(as.character(v))
      v_clean <- v_clean[!is.na(v_clean) & nzchar(v_clean)]
      tibble::tibble(
        cluster  = cl,
        gene_set = nm,
        n_input  = length(v),
        n_clean  = length(v_clean),
        n_used   = sum(v_clean %in% names(stats)) # ilu z listy było w stats
      )
    })
  }
  gene_counts <- dplyr::bind_rows(gene_counts)
  
  # --- 4) Summaries ---
  first_by_NES <- res_list %>% purrr::imap(~ if (nrow(.x)) dplyr::slice_max(.x, .x$NES, n = 1, with_ties = TRUE) else NULL) %>%
    purrr::discard(is.null) %>% dplyr::bind_rows()
  first_by_pval <- res_list %>% purrr::imap(~ if (nrow(.x)) dplyr::slice_min(.x, .x$pval, n = 1, with_ties = TRUE) else NULL) %>%
    purrr::discard(is.null) %>% dplyr::bind_rows()
  first_by_padj <- res_list %>% purrr::imap(~ if (nrow(.x)) dplyr::slice_min(.x, .x$padj, n = 1, with_ties = TRUE) else NULL) %>%
    purrr::discard(is.null) %>% dplyr::bind_rows()
  
  # --- 5) Fisher (optional) ---
  fisher_out <- NULL
  if (isTRUE(do_fisher)) {
    if (!exists("compute_geneSet_fisher")) stop("`compute_geneSet_fisher()` not found.")
    if (is.null(fisher_universe)) fisher_universe <- rownames(mat)
    if (verbose) message("Running Fisher tests on marker sets …")
    fisher_out <- compute_geneSet_fisher(
      gene_lists = cleaned_lists, universe = fisher_universe,
      min_size = fisher_min_size, alternative = fisher_alternative,
      strict_checks = fisher_strict_checks, smooth_or = fisher_smooth_or
    )
  }
  
  # --- 6) XLSX export ---
  if (isTRUE(export_xlsx)) {
    sanitize_sheet_name <- function(nm) {
      nm <- gsub("[\\*/\\?\\[\\]:]", "_", nm); if (nchar(nm) > 31) substr(nm, 1, 31) else nm
    }
    if (is.null(xlsx_path)) xlsx_path <- file.path(getwd(), paste0("fgsea_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx"))
    
    wb <- openxlsx::createWorkbook()
    
    add_sheet_df <- function(wb, sheet_name, df, digits = 3, float_cols = xlsx_float_cols) {
      sheet_name <- sanitize_sheet_name(sheet_name)
      openxlsx::addWorksheet(wb, sheet_name)
      openxlsx::writeData(wb, sheet = sheet_name, x = as.data.frame(df), keepNA = TRUE)
      openxlsx::freezePane(wb, sheet = sheet_name, firstRow = TRUE)
      openxlsx::addFilter(wb, sheet = sheet_name, rows = 1, cols = 1:ncol(df))
      openxlsx::setColWidths(wb, sheet = sheet_name, cols = 1:ncol(df), widths = "auto")
      present <- intersect(float_cols, colnames(df))
      if (length(present) > 0 && nrow(df) > 0) {
        num_fmt <- if (digits <= 0) "0" else paste0("0.", paste(rep("0", digits), collapse = ""))
        num_style <- openxlsx::createStyle(numFmt = num_fmt)
        float_idx <- match(present, colnames(df))
        openxlsx::addStyle(wb, sheet = sheet_name, style = num_style,
                           cols = float_idx, rows = 2:(nrow(df) + 1), gridExpand = TRUE)
      }
    }
    
    if (write_first_by_NES  && nrow(first_by_NES))  add_sheet_df(wb, "first_by_NES",  first_by_NES,  digits = digits_xlsx)
    if (write_first_by_pval && nrow(first_by_pval)) add_sheet_df(wb, "first_by_pval", first_by_pval, digits = digits_xlsx)
    if (write_first_by_padj && nrow(first_by_padj)) add_sheet_df(wb, "first_by_padj", first_by_padj, digits = digits_xlsx)
    for (cl in cluster_order) {
      add_sheet_df(wb, cl, res_list[[cl]], digits = digits_xlsx)
    }
    add_sheet_df(wb, "gene_counts", gene_counts, digits = 0, float_cols = character())
    
    if (verbose) message("Writing XLSX to: ", xlsx_path)
    openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)
  }
  
  # --- 7) Return ---
  out <- list(
    results           = res_list,
    raw_results       = if (keep_raw_results) raw_list else NULL,
    rank_stats        = if (store_rank_stats) stats_map else NULL,
    first_by_NES      = first_by_NES,
    first_by_pval     = first_by_pval,
    first_by_padj     = first_by_padj,
    gene_lists_used   = cleaned_lists,
    gene_counts       = gene_counts
  )
  if (!is.null(fisher_out)) out$fisher <- fisher_out
  if (isTRUE(export_xlsx))  out$xlsx_path <- xlsx_path
  out
}

run_celltype_fgseaEnrichemnt_multClusters <- function(
    seurat_obj,
    cluster_cols,                                   # wektor kolumn klastrów (np. różne rozdzielczości)
    gene_lists,
    sample_col      = "sample_id",
    assay           = "RNA",
    layer           = "counts",
    biomart_version = "114",
    nperm   = 1000,
    minSize = 5,
    maxSize = Inf,
    seed    = 1,
    verbose = TRUE,
    do_fisher            = FALSE,
    fisher_universe      = NULL,
    fisher_min_size      = 1,
    fisher_alternative   = "greater",
    fisher_strict_checks = TRUE,
    fisher_smooth_or     = FALSE,
    rename_pathway_to    = "predict_cell_type",
    rename_size_to       = "n_genes_in_set",
    add_leadingEdge_n    = TRUE,
    export_xlsx          = FALSE,
    xlsx_path            = NULL,                    # może być: ścieżka z {cluster}, katalog, albo plik
    write_first_by_NES   = TRUE,
    write_first_by_pval  = TRUE,
    write_first_by_padj  = TRUE,
    digits_xlsx          = 3,
    xlsx_float_cols      = c("pval","padj","ES","NES"),
    keep_raw_results     = TRUE,
    store_rank_stats     = TRUE,
    export_raw_in_xlsx   = FALSE
) {
  # --- sanity checks ---
  if (!exists("prepare_proteincoding_avg_expr")) stop("Brak funkcji prepare_proteincoding_avg_expr().")
  if (!exists("calculate_gene_specificity"))     stop("Brak funkcji calculate_gene_specificity().")
  if (!exists("run_fgsea_on_score_matrix"))      stop("Brak funkcji run_fgsea_on_score_matrix().")
  if (!exists("download_biomart_genes"))         stop("Brak funkcji download_biomart_genes().")
  
  # --- 0) Pobierz geny z BioMarta ---
  bm <- download_biomart_genes(version = biomart_version)
  
  # wybierz właściwą kolumnę symboli
  symbol_col <- if ("hgnc_symbol" %in% colnames(bm)) {
    "hgnc_symbol"
  } else if ("external_gene_name" %in% colnames(bm)) {
    "external_gene_name"
  } else {
    stop("Brak kolumny z nazwą genu (oczekiwano 'hgnc_symbol' albo 'external_gene_name').")
  }
  
  biomart_genes <- bm %>%
    dplyr::filter(
      chromosome_name %in% c(as.character(1:22), "X", "Y"),
      gene_biotype == "protein_coding"
    ) %>%
    dplyr::pull(!!sym(symbol_col)) %>%
    unique()
  
  if (verbose) message("BioMart reference contains ", length(biomart_genes), " protein-coding genes (autosomes + X/Y).")
  
  # --- 1) Przefiltruj gene_lists przez biomart ---
  gene_lists_filtered <- purrr::imap(gene_lists, function(v, nm) {
    v_clean <- unique(as.character(v))
    v_clean <- v_clean[!is.na(v_clean) & nzchar(v_clean)]
    intersect(v_clean, biomart_genes)
  })
  
  # --- 2) Pętla po kolumnach klastrów ---
  out_list <- setNames(vector("list", length(cluster_cols)), cluster_cols)
  
  for (cl_col in cluster_cols) {
    if (verbose) message("\n=== Processing cluster_col: ", cl_col, " ===")
    
    # 2.1 Średnie ekspresje (spójne z biomart)
    avg_mat <- prepare_proteincoding_avg_expr(
      seurat_obj       = seurat_obj,
      cluster_col      = cl_col,
      sample_col       = sample_col,
      assay            = assay,
      layer            = layer,
      biomart_version  = biomart_version,
      verbose          = verbose
    )
    
    # 2.2 Specyficzność genów
    gene_specificity <- calculate_gene_specificity(avg_mat, na_omit = TRUE)
    
    # 2.3 Helper do ścieżki XLSX
    .build_xlsx_path <- function(base, cluster_col) {
      if (is.null(base)) return(NULL)
      if (dir.exists(base)) {
        fname <- paste0("fgsea_", cluster_col, ".xlsx")
        return(file.path(base, fname))
      }
      if (grepl("\\{cluster\\}", base)) {
        out <- sub("\\{cluster\\}", cluster_col, base)
        dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
        return(out)
      }
      dir.create(dirname(base), recursive = TRUE, showWarnings = FALSE)
      bn  <- basename(base)
      ext <- if (grepl("\\.", bn)) sub(".*\\.", "", bn) else "xlsx"
      nm  <- if (grepl("\\.", bn)) sub("\\.[^.]*$", "", bn) else bn
      file.path(dirname(base), paste0(nm, "__", cluster_col, ".", ext))
    }
    xlsx_path_i <- if (isTRUE(export_xlsx)) .build_xlsx_path(xlsx_path, cl_col) else NULL
    
    # 2.4 fgsea
    fgsea_out <- run_fgsea_on_score_matrix(
      mat                   = gene_specificity$combined,
      gene_lists            = gene_lists_filtered,
      nperm                 = nperm,
      minSize               = minSize,
      maxSize               = maxSize,
      seed                  = seed,
      verbose               = verbose,
      do_fisher             = do_fisher,
      fisher_universe       = fisher_universe,
      fisher_min_size       = fisher_min_size,
      fisher_alternative    = fisher_alternative,
      fisher_strict_checks  = fisher_strict_checks,
      fisher_smooth_or      = fisher_smooth_or,
      rename_pathway_to     = rename_pathway_to,
      rename_size_to        = rename_size_to,
      add_leadingEdge_n     = add_leadingEdge_n,
      export_xlsx           = export_xlsx,
      xlsx_path             = xlsx_path_i,
      write_first_by_NES    = write_first_by_NES,
      write_first_by_pval   = write_first_by_pval,
      write_first_by_padj   = write_first_by_padj,
      digits_xlsx           = digits_xlsx,
      xlsx_float_cols       = xlsx_float_cols,
      keep_raw_results      = keep_raw_results,
      store_rank_stats      = store_rank_stats,
      export_raw_in_xlsx    = export_raw_in_xlsx
    )
    
    # 2.5 gene_counts raport
    gene_counts <- purrr::imap_dfr(gene_lists, function(v, nm) {
      v_clean <- unique(as.character(v))
      v_clean <- v_clean[!is.na(v_clean) & nzchar(v_clean)]
      tibble::tibble(
        cluster      = cl_col,
        gene_set     = nm,
        n_input      = length(v),
        n_in_biomart = sum(v_clean %in% biomart_genes),
        n_used       = sum(v_clean %in% rownames(avg_mat))
      )
    })
    
    # 2.6 Zapisz do out_list
    out_list[[cl_col]] <- list(
      cluster_col      = cl_col,
      avg_mat          = avg_mat,
      gene_specificity = gene_specificity,
      fgsea            = fgsea_out,
      gene_counts      = gene_counts,
      xlsx_path        = if (!is.null(fgsea_out$xlsx_path)) fgsea_out$xlsx_path else xlsx_path_i
    )
  }
  
  out_list
}

