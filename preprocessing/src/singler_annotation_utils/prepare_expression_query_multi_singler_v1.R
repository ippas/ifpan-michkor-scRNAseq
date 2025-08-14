# =========================
# 1) Query: wiele kolumn klastrów naraz
# =========================
prepare_expression_query_multi_singler_v1 <- function(
    seurat_obj,
    assay = "RNA",
    data_layer_regex = "^data",
    normalize_if_missing = TRUE,
    gene_name_from_rownames = TRUE,
    gene_sep = "-",
    agg_fun_for_duplicates = c("mean","sum"),
    cluster_cols,                       # c("seurat_clusters_res0.1","seurat_clusters_res0.3", ...)
    alias_first_as = "cluster",         # pierwszy wpis z cluster_cols będzie też skopiowany pod tą nazwę
    save_rds = NULL,
    verbose = TRUE
) {
  start_time <- Sys.time()
  .msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))
  agg_fun_for_duplicates <- match.arg(agg_fun_for_duplicates)
  
  stopifnot(inherits(seurat_obj, "Seurat"))
  if (missing(cluster_cols) || length(cluster_cols) < 1)
    stop("Provide at least one clustering column in 'cluster_cols'.")
  
  DefaultAssay(seurat_obj) <- assay
  
  # [1/8] Warstwy danych
  .msg("[1/8] Checking data layers in assay '%s'…", assay)
  if (!assay %in% names(seurat_obj@assays))
    stop(sprintf("Assay '%s' not found in Seurat object.", assay))
  data_layers <- grep(data_layer_regex, Layers(seurat_obj[[assay]]), value = TRUE)
  if (length(data_layers) == 0) {
    if (isTRUE(normalize_if_missing)) {
      .msg("       No data layers matched '%s' → running NormalizeData()…", data_layer_regex)
      seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
      data_layers <- "data"
    } else stop("No data layers found and normalize_if_missing=FALSE.")
  }
  .msg("       Using layers: %s", paste(data_layers, collapse = ", "))
  
  # [2/8] Pobranie macierzy i sklejenie
  .msg("[2/8] Extracting matrices from layers and binding columns…")
  mat_list <- lapply(data_layers, function(ly) LayerData(seurat_obj[[assay]], layer = ly))
  logmat <- do.call(cbind, mat_list)
  logmat <- methods::as(logmat, "dgCMatrix")
  .msg("       Raw matrix dims: genes=%d, cells=%d", nrow(logmat), ncol(logmat))
  
  # [3/8] Nazwy genów + duplikaty
  .msg("[3/8] Processing gene names%s…",
       ifelse(gene_name_from_rownames, " (converting & aggregating duplicates)", ""))
  rn <- rownames(logmat)
  if (isTRUE(gene_name_from_rownames)) {
    gene_symbols <- ifelse(grepl(gene_sep, rn, fixed = TRUE),
                           sub(paste0("^[^", gene_sep, "]+", gene_sep), "", rn),
                           rn)
  } else gene_symbols <- rn
  
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
  
  # [4/8] Budowa SCE
  .msg("[4/8] Building SingleCellExperiment with 'logcounts'…")
  sce_query <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = logmat))
  
  # [5/8] Walidacja kolumn klastrów
  .msg("[5/8] Validating clustering columns…")
  missing_cols <- setdiff(cluster_cols, colnames(seurat_obj@meta.data))
  if (length(missing_cols))
    stop(sprintf("Missing clustering columns in meta.data: %s",
                 paste(missing_cols, collapse = ", ")))
  
  # [6/8] Spójność nazw komórek
  .msg("[6/8] Checking cell name alignment…")
  if (!all(colnames(sce_query) %in% rownames(seurat_obj@meta.data))) {
    n_miss <- sum(!colnames(sce_query) %in% rownames(seurat_obj@meta.data))
    stop(sprintf("Cell name mismatch: %d query cells not found in Seurat meta.data.", n_miss))
  }
  
  # [7/8] Dołącz wszystkie kolumny klastrów; pierwszą skopiuj pod alias
  .msg("[7/8] Attaching %d clustering columns…", length(cluster_cols))
  for (cc in cluster_cols) {
    sce_query[[cc]] <- seurat_obj@meta.data[colnames(sce_query), cc, drop = TRUE]
  }
  if (!is.null(alias_first_as) && nzchar(alias_first_as)) {
    sce_query[[alias_first_as]] <- sce_query[[cluster_cols[1]]]
    .msg("       Aliased '%s' as '%s'.", cluster_cols[1], alias_first_as)
  }
  
  # [8/8] Sanity + zapis
  .msg("[8/8] Sanity checks…")
  if (!"logcounts" %in% SummarizedExperiment::assayNames(sce_query))
    stop("'logcounts' assay missing in query SCE.")
  .msg("       OK: logcounts present; clustering columns attached: %s",
       paste(c(alias_first_as, cluster_cols), collapse = ", "))
  
  if (!is.null(save_rds)) {
    .msg("Saving query SCE (multi-clustering) to RDS (gzip): %s", save_rds)
    saveRDS(sce_query, file = save_rds, compress = "gzip")
  }
  
  end_time <- Sys.time()
  .msg("Done. Query SCE dims: genes=%d, cells=%d; colData cols=%d. Time: %.2f s.",
       nrow(sce_query), ncol(sce_query), ncol(SummarizedExperiment::colData(sce_query)),
       as.numeric(difftime(end_time, start_time, units = "secs")))
  sce_query
}
