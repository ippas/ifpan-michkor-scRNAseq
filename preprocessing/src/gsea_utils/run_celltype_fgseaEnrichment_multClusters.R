run_celltype_fgseaEnrichemnt_multClusters <- function(
    seurat_obj,
    cluster_cols,                                   # wektor kolumn klastrów (np. różne rozdzielczości)
    gene_lists,
    # --- prepare_proteincoding_avg_expr ---
    sample_col      = "sample_id",
    assay           = "RNA",
    layer           = "counts",
    biomart_version = "114",
    # --- fgsea/run_fgsea_on_score_matrix ---
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
  # sprawdzenie dostępności funkcji używanych w pipeline
  if (!exists("prepare_proteincoding_avg_expr")) stop("Brak funkcji prepare_proteincoding_avg_expr().")
  if (!exists("calculate_gene_specificity"))     stop("Brak funkcji calculate_gene_specificity().")
  if (!exists("run_fgsea_on_score_matrix"))      stop("Brak funkcji run_fgsea_on_score_matrix().")
  
  # helper budujący ścieżkę XLSX per cluster_col
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
  
  out_list <- setNames(vector("list", length(cluster_cols)), cluster_cols)
  
  for (cl_col in cluster_cols) {
    if (verbose) message("\n=== Processing cluster_col: ", cl_col, " ===")
    
    # 1) średnie ekspresje dla protein-coding wg danej kolumny klastrów
    avg_mat <- prepare_proteincoding_avg_expr(
      seurat_obj       = seurat_obj,
      cluster_col      = cl_col,
      sample_col       = sample_col,
      assay            = assay,
      layer            = layer,
      biomart_version  = biomart_version,
      verbose          = verbose,
      select_samples   = select_samples,  
    )
    
    # 2) specyficzność genów (max-norm × tau)
    gene_specificity <- calculate_gene_specificity(avg_mat, na_omit = TRUE)
    
    # 3) fgsea + (opcjonalnie) XLSX
    xlsx_path_i <- if (isTRUE(export_xlsx)) .build_xlsx_path(xlsx_path, cl_col) else NULL
    
    fgsea_out <- run_fgsea_on_score_matrix(
      mat                   = gene_specificity$combined,
      gene_lists            = gene_lists,
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
    
    out_list[[cl_col]] <- list(
      cluster_col      = cl_col,
      avg_mat          = avg_mat,
      gene_specificity = gene_specificity,
      fgsea            = fgsea_out,
      xlsx_path        = if (!is.null(fgsea_out$xlsx_path)) fgsea_out$xlsx_path else xlsx_path_i
    )
  }
  
  out_list
}


run_celltype_fgseaEnrichemnt_multClusters <- function(
    seurat_obj,
    cluster_cols,                                   # wektor kolumn klastrów (np. różne rozdzielczości)
    gene_lists,
    sample_col      = "sample_id",
    assay           = "RNA",
    layer           = "counts",
    biomart_version = "114",
    select_samples  = NULL,                         # <--- NOWE
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
    
    # 2.1 Średnie ekspresje (spójne z biomart + opcjonalne select_samples)
    avg_mat <- prepare_proteincoding_avg_expr(
      seurat_obj       = seurat_obj,
      cluster_col      = cl_col,
      sample_col       = sample_col,
      assay            = assay,
      layer            = layer,
      biomart_version  = biomart_version,
      verbose          = verbose,
      select_samples   = select_samples   # <--- PRZYWRÓCONE
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

