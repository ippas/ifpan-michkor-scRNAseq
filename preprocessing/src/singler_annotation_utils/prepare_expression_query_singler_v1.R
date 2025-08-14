prepare_expression_query_singler_v1 <- function(
    seurat_obj,                         # Seurat v5 object (your data)
    assay = "RNA",
    data_layer_regex = "^data",         # which layers to use (e.g., "data")
    normalize_if_missing = TRUE,        # run NormalizeData() if no data-layer
    gene_name_from_rownames = TRUE,     # convert "ENSG...-SYMBOL" -> "SYMBOL"
    gene_sep = "-",
    agg_fun_for_duplicates = c("mean","sum"),
    cluster_col,                        # e.g. "seurat_clusters_res0.17"
    save_rds = NULL,
    verbose = TRUE
) {
  start_time <- Sys.time()
  .msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))
  
  stopifnot(inherits(seurat_obj, "Seurat"))
  DefaultAssay(seurat_obj) <- assay
  agg_fun_for_duplicates <- match.arg(agg_fun_for_duplicates)
  
  # 1) Locate data layers or normalize if missing
  .msg("[1/7] Checking data layers in assay '%s'…", assay)
  if (!assay %in% names(seurat_obj@assays)) {
    stop(sprintf("Assay '%s' not found in Seurat object.", assay))
  }
  data_layers <- grep(data_layer_regex, Layers(seurat_obj[[assay]]), value = TRUE)
  if (length(data_layers) == 0) {
    if (isTRUE(normalize_if_missing)) {
      .msg("       No data layers matched '%s' → running NormalizeData()…", data_layer_regex)
      seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
      data_layers <- "data"
    } else {
      stop("No data layers found and normalize_if_missing=FALSE.")
    }
  }
  .msg("       Using layers: %s", paste(data_layers, collapse=", "))
  
  # 2) Collate layers into one log-expression matrix (genes x cells)
  .msg("[2/7] Extracting matrices from layers and binding columns…")
  mat_list <- lapply(data_layers, function(ly) LayerData(seurat_obj[[assay]], layer = ly))
  logmat <- do.call(cbind, mat_list)
  logmat <- methods::as(logmat, "dgCMatrix")
  .msg("       Raw matrix dims: genes=%d, cells=%d", nrow(logmat), ncol(logmat))
  
  # 3) Optionally convert gene IDs to symbols and aggregate duplicates
  .msg("[3/7] Processing gene names%s…",
       ifelse(gene_name_from_rownames, " (converting & aggregating duplicates)", ""))
  rn <- rownames(logmat)
  if (isTRUE(gene_name_from_rownames)) {
    # convert only if separator present; otherwise keep as-is
    gene_symbols <- ifelse(grepl(gene_sep, rn, fixed = TRUE),
                           sub(paste0("^[^", gene_sep, "]+", gene_sep), "", rn),
                           rn)
  } else {
    gene_symbols <- rn
  }
  
  dup_idx <- duplicated(gene_symbols) | duplicated(gene_symbols, fromLast = TRUE)
  mat_nondup <- logmat[!dup_idx, , drop = FALSE]
  rownames(mat_nondup) <- gene_symbols[!dup_idx]
  
  if (any(dup_idx)) {
    .msg("       Found %d duplicated gene symbols → aggregating by %s…",
         sum(dup_idx), match.arg(agg_fun_for_duplicates))
    aggregate_sparse_matrix <- function(mat, groups, fun = c("mean","sum")) {
      fun <- match.arg(fun)
      groups <- factor(groups, levels = unique(groups))
      G <- Matrix::sparse.model.matrix(~ groups - 1)
      summed <- Matrix::crossprod(G, mat)
      if (fun == "mean") {
        counts <- as.numeric(table(groups))
        summed <- summed / counts
      }
      rownames(summed) <- levels(groups)
      methods::as(summed, "dgCMatrix")
    }
    mat_dup <- logmat[dup_idx, , drop = FALSE]
    mat_dup_agg <- aggregate_sparse_matrix(mat_dup, groups = gene_symbols[dup_idx],
                                           fun = agg_fun_for_duplicates)
    logmat <- rbind(mat_nondup, mat_dup_agg)
  } else {
    logmat <- mat_nondup
  }
  .msg("       Final matrix dims: genes=%d, cells=%d", nrow(logmat), ncol(logmat))
  
  # 4) Build SCE with logcounts
  .msg("[4/7] Building SingleCellExperiment with 'logcounts'…")
  sce_test <- SingleCellExperiment(assays = list(logcounts = logmat))
  
  # 5) Attach clusters from Seurat meta.data
  # [5/7] Attach clusters from Seurat meta.data (by cell barcodes)
  .msg("[5/7] Attaching cluster column '%s' from Seurat metadata…", cluster_col)
  
  if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
    stop(sprintf("Cluster column '%s' not found in Seurat meta.data.", cluster_col))
  }
  
  # sanity: czy nazwy komórek z SCE są w meta.data?
  if (!all(colnames(sce_test) %in% rownames(seurat_obj@meta.data))) {
    n_miss <- sum(!colnames(sce_test) %in% rownames(seurat_obj@meta.data))
    stop(sprintf("Cell name mismatch: %d query cells not found in Seurat meta.data.", n_miss))
  }
  
  clusters_vec <- seurat_obj@meta.data[colnames(sce_test), cluster_col, drop = TRUE]
  sce_test$cluster <- clusters_vec
  
  # 6) Basic sanity checks
  .msg("[6/7] Sanity checks…")
  if (!"logcounts" %in% assayNames(sce_test)) stop("'logcounts' assay missing in sce_test.")
  if (!"cluster"   %in% colnames(colData(sce_test))) stop("'cluster' column missing in sce_test colData.")
  .msg("       OK: logcounts present, cluster column present.")
  
  # 7) Save (optional)
  if (!is.null(save_rds)) {
    .msg("[7/7] Saving test SCE to RDS: %s", save_rds)
    saveRDS(sce_test, file = save_rds, compress = "gzip")
  }
  
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "secs")
  .msg("Done. Test SCE dims: genes=%d, cells=%d; colData cols=%d.", 
       nrow(sce_test), ncol(sce_test), ncol(colData(sce_test)))
  .msg("Total execution time: %.2f seconds.", as.numeric(elapsed))
  
  sce_test
}




prepare_expression_query_singler_v1 <- function(
    seurat_obj,                         # Seurat v5 object (Twoje dane)
    assay = "RNA",
    data_layer_regex = "^data",         # które warstwy wykorzystać (np. "data")
    normalize_if_missing = TRUE,        # uruchom NormalizeData() jeśli brak warstwy data
    gene_name_from_rownames = TRUE,     # konwersja "ENSG...-SYMBOL" -> "SYMBOL"
    gene_sep = "-",
    agg_fun_for_duplicates = c("mean","sum"),
    cluster_col,                        # np. "seurat_clusters_res0.17"
    save_rds = NULL,                    # ścieżka do zapisu RDS (opcjonalnie)
    verbose = TRUE
) {
  start_time <- Sys.time()
  .msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))
  agg_fun_for_duplicates <- match.arg(agg_fun_for_duplicates)
  
  stopifnot(inherits(seurat_obj, "Seurat"))
  DefaultAssay(seurat_obj) <- assay
  
  # [1/7] Warstwy danych
  .msg("[1/7] Checking data layers in assay '%s'…", assay)
  if (!assay %in% names(seurat_obj@assays)) {
    stop(sprintf("Assay '%s' not found in Seurat object.", assay))
  }
  data_layers <- grep(data_layer_regex, Layers(seurat_obj[[assay]]), value = TRUE)
  if (length(data_layers) == 0) {
    if (isTRUE(normalize_if_missing)) {
      .msg("       No data layers matched '%s' → running NormalizeData()…", data_layer_regex)
      seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
      data_layers <- "data"
    } else {
      stop("No data layers found and normalize_if_missing=FALSE.")
    }
  }
  .msg("       Using layers: %s", paste(data_layers, collapse = ", "))
  
  # [2/7] Pobranie macierzy i sklejenie
  .msg("[2/7] Extracting matrices from layers and binding columns…")
  mat_list <- lapply(data_layers, function(ly) LayerData(seurat_obj[[assay]], layer = ly))
  logmat <- do.call(cbind, mat_list)
  logmat <- methods::as(logmat, "dgCMatrix")
  .msg("       Raw matrix dims: genes=%d, cells=%d", nrow(logmat), ncol(logmat))
  
  # [3/7] Nazwy genów i duplikaty
  .msg("[3/7] Processing gene names%s…",
       ifelse(gene_name_from_rownames, " (converting & aggregating duplicates)", ""))
  rn <- rownames(logmat)
  if (isTRUE(gene_name_from_rownames)) {
    gene_symbols <- ifelse(grepl(gene_sep, rn, fixed = TRUE),
                           sub(paste0("^[^", gene_sep, "]+", gene_sep), "", rn),
                           rn)
  } else {
    gene_symbols <- rn
  }
  
  dup_idx <- duplicated(gene_symbols) | duplicated(gene_symbols, fromLast = TRUE)
  mat_nondup <- logmat[!dup_idx, , drop = FALSE]
  rownames(mat_nondup) <- gene_symbols[!dup_idx]
  
  if (any(dup_idx)) {
    .msg("       Found %d duplicated gene symbols → aggregating by %s…",
         sum(dup_idx), agg_fun_for_duplicates)
    aggregate_sparse_matrix <- function(mat, groups, fun = c("mean","sum")) {
      fun <- match.arg(fun)
      groups <- factor(groups, levels = unique(groups))
      G <- Matrix::sparse.model.matrix(~ groups - 1)
      summed <- Matrix::crossprod(G, mat)
      if (fun == "mean") {
        counts <- as.numeric(table(groups))
        summed <- summed / counts
      }
      rownames(summed) <- levels(groups)
      methods::as(summed, "dgCMatrix")
    }
    mat_dup <- logmat[dup_idx, , drop = FALSE]
    mat_dup_agg <- aggregate_sparse_matrix(mat_dup, groups = gene_symbols[dup_idx],
                                           fun = agg_fun_for_duplicates)
    logmat <- rbind(mat_nondup, mat_dup_agg)
  } else {
    logmat <- mat_nondup
  }
  .msg("       Final matrix dims: genes=%d, cells=%d", nrow(logmat), ncol(logmat))
  
  # [4/7] Budowa SCE z logcounts
  .msg("[4/7] Building SingleCellExperiment with 'logcounts'…")
  sce_test <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = logmat))
  
  # [5/7] Dołączenie klastrów (po nazwach komórek = rownames(meta))
  .msg("[5/7] Attaching cluster column '%s' from Seurat metadata…", cluster_col)
  if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
    stop(sprintf("Cluster column '%s' not found in Seurat meta.data.", cluster_col))
  }
  if (!all(colnames(sce_test) %in% rownames(seurat_obj@meta.data))) {
    n_miss <- sum(!colnames(sce_test) %in% rownames(seurat_obj@meta.data))
    stop(sprintf("Cell name mismatch: %d query cells not found in Seurat meta.data.", n_miss))
  }
  clusters_vec <- seurat_obj@meta.data[colnames(sce_test), cluster_col, drop = TRUE]
  sce_test$cluster <- clusters_vec
  
  # [6/7] Sanity
  .msg("[6/7] Sanity checks…")
  if (!"logcounts" %in% SummarizedExperiment::assayNames(sce_test))
    stop("'logcounts' assay missing in sce_test.")
  if (!"cluster" %in% colnames(SummarizedExperiment::colData(sce_test)))
    stop("'cluster' column missing in sce_test colData.")
  .msg("       OK: logcounts present, cluster column present.")
  
  # [7/7] Zapis (opcjonalnie) — szybciej: gzip
  if (!is.null(save_rds)) {
    .msg("[7/7] Saving query SCE to RDS (gzip compression): %s", save_rds)
    saveRDS(sce_test, file = save_rds, compress = "gzip")
  }
  
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "secs")
  .msg("Done. Query SCE dims: genes=%d, cells=%d; colData cols=%d.",
       nrow(sce_test), ncol(sce_test), ncol(SummarizedExperiment::colData(sce_test)))
  .msg("Total execution time: %.2f seconds.", as.numeric(elapsed))
  
  sce_test
}


