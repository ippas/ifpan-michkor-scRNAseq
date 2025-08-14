predict_celltypes_singler_by_cluster_v1 <- function(
    sce_query,                          # SCE z 'logcounts' i kolumną 'cluster'
    sce_ref,                            # SCE referencyjny z kolumną etykiet
    label_col = "Cell.Type.Full",       # kolumna z etykietami w referencji
    cluster_col = "cluster",            # kolumna z klastrami w query
    min_common_genes = 1000,
    assay.type.test = "logcounts",
    assay.type.ref  = "logcounts",
    save_dir = NULL,                    # jeśli podasz: zapisze RDS + CSV
    verbose = TRUE
) {
  start_time <- Sys.time()
  .msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))
  
  # [1/6] Walidacje wejść
  .msg("[1/6] Validating inputs…")
  if (!label_col %in% colnames(SummarizedExperiment::colData(sce_ref))) {
    stop(sprintf("Label column '%s' not found in reference colData.", label_col))
  }
  if (!cluster_col %in% colnames(SummarizedExperiment::colData(sce_query))) {
    stop(sprintf("Cluster column '%s' not found in query colData.", cluster_col))
  }
  if (!assay.type.ref %in% SummarizedExperiment::assayNames(sce_ref)) {
    stop(sprintf("Assay '%s' missing in reference.", assay.type.ref))
  }
  if (!assay.type.test %in% SummarizedExperiment::assayNames(sce_query)) {
    stop(sprintf("Assay '%s' missing in query.", assay.type.test))
  }
  
  # [2/6] Wspólne geny
  .msg("[2/6] Intersecting genes between reference and query…")
  common <- intersect(rownames(sce_ref), rownames(sce_query))
  .msg("       Common genes: %d", length(common))
  if (length(common) < min_common_genes) {
    stop(sprintf("Too few common genes: %d (< %d).", length(common), min_common_genes))
  }
  ref2   <- sce_ref  [common, ]
  query2 <- sce_query[common, ]
  
  # [3/6] Uruchomienie SingleR (per klaster)
  .msg("[3/6] Running SingleR per cluster…")
  pred <- SingleR::SingleR(
    test            = query2,
    ref             = ref2,
    labels          = SummarizedExperiment::colData(ref2)[[label_col]],
    clusters        = SummarizedExperiment::colData(query2)[[cluster_col]],
    assay.type.test = assay.type.test,
    assay.type.ref  = assay.type.ref
  )
  
  # [4/6] Tabela podsumowująca
  .msg("[4/6] Building summary table…")
  annot_tbl <- data.frame(
    cluster       = rownames(pred),
    label         = pred$labels,
    pruned_label  = if (!is.null(pred$pruned.labels)) pred$pruned.labels else pred$labels,
    best_score    = apply(pred$scores, 1, max),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  annot_tbl <- annot_tbl[order(annot_tbl$cluster), , drop = FALSE]
  
  # [5/6] Zapis (opcjonalnie)
  if (!is.null(save_dir)) {
    .msg("[5/6] Saving outputs to: %s", save_dir)
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
    saveRDS(pred, file = file.path(save_dir, "singler_predictions.rds"), compress = "gzip")
    utils::write.csv(annot_tbl, file = file.path(save_dir, "singler_annotation_table.csv"), row.names = FALSE)
    
    # (opcjonalnie) zapis macierzy score’ów per klaster
    scores_path <- file.path(save_dir, "singler_scores_matrix.csv")
    utils::write.csv(pred$scores, file = scores_path, row.names = TRUE)
  }
  
  # [6/6] Podsumowanie
  end_time <- Sys.time()
  elapsed  <- difftime(end_time, start_time, units = "secs")
  .msg("Done. Predicted %d clusters. Total time: %.2f s.",
       nrow(annot_tbl), as.numeric(elapsed))
  
  list(
    pred = pred,                    # obiekt SingleR z pełnymi wynikami
    annotation_table = annot_tbl,   # zgrabna tabelka: cluster/label/pruned_label/best_score
    n_common_genes = length(common)
  )
}


predict_celltypes_singler_multi_v1 <- function(
    sce_query,
    sce_ref,
    label_col = "Cell.Type.Full",       # label column in reference colData
    cluster_cols = NULL,                # NULL -> auto-detect columns starting with "cluster"
    min_common_genes = 1000,
    assay.type.test = "logcounts",
    assay.type.ref  = "logcounts",
    save_dir = NULL,                    # if provided: save per-clustering outputs under this dir
    save_scores = TRUE,                 # also save the per-cluster score matrix
    verbose = TRUE
) {
  start_time <- Sys.time()
  .msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))
  
  # [1/7] Validate inputs
  .msg("[1/7] Validating inputs…")
  if (!label_col %in% colnames(SummarizedExperiment::colData(sce_ref))) {
    stop(sprintf("Label column '%s' not found in reference colData.", label_col))
  }
  if (!assay.type.ref %in% SummarizedExperiment::assayNames(sce_ref)) {
    stop(sprintf("Assay '%s' missing in reference.", assay.type.ref))
  }
  if (!assay.type.test %in% SummarizedExperiment::assayNames(sce_query)) {
    stop(sprintf("Assay '%s' missing in query.", assay.type.test))
  }
  
  # Auto-detect cluster columns if not provided
  if (is.null(cluster_cols)) {
    cluster_cols <- grep("^cluster", colnames(SummarizedExperiment::colData(sce_query)), value = TRUE)
    if (length(cluster_cols) == 0) {
      stop("No cluster columns found (expected columns starting with 'cluster').")
    }
    .msg("       Auto-detected cluster columns: %s", paste(cluster_cols, collapse = ", "))
  } else {
    missing_cc <- setdiff(cluster_cols, colnames(SummarizedExperiment::colData(sce_query)))
    if (length(missing_cc)) {
      stop(sprintf("Missing cluster columns in query SCE: %s", paste(missing_cc, collapse = ", ")))
    }
  }
  
  # [2/7] Check gene overlap and subset once
  .msg("[2/7] Checking gene overlap between reference and query datasets…")
  .msg("       Reference genes: %d", nrow(sce_ref))
  .msg("       Query genes: %d", nrow(sce_query))
  common <- intersect(rownames(sce_ref), rownames(sce_query))
  .msg("       Common genes: %d", length(common))
  if (length(common) < min_common_genes) {
    stop(sprintf("Too few common genes: %d (< %d).", length(common), min_common_genes))
  }
  ref2   <- sce_ref  [common, ]
  query2 <- sce_query[common, ]
  
  # [3/7] Run SingleR for each clustering column
  .msg("[3/7] Running SingleR for %d clustering columns…", length(cluster_cols))
  results_list <- vector("list", length(cluster_cols))
  names(results_list) <- cluster_cols
  
  if (!is.null(save_dir)) {
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  for (cc in cluster_cols) {
    .msg("       → Running SingleR for clustering: %s", cc)
    
    pred <- SingleR::SingleR(
      test            = query2,
      ref             = ref2,
      labels          = SummarizedExperiment::colData(ref2)[[label_col]],
      clusters        = SummarizedExperiment::colData(query2)[[cc]],
      assay.type.test = assay.type.test,
      assay.type.ref  = assay.type.ref
    )
    
    # Build a tidy summary table per clustering
    annot_tbl <- data.frame(
      clustering    = cc,
      cluster       = rownames(pred),
      label         = pred$labels,
      pruned_label  = if (!is.null(pred$pruned.labels)) pred$pruned.labels else pred$labels,
      best_score    = apply(pred$scores, 1, max),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    annot_tbl <- annot_tbl[order(annot_tbl$cluster), , drop = FALSE]
    
    # Optional save: per-clustering folder with RDS + CSV (and scores)
    if (!is.null(save_dir)) {
      subdir <- file.path(save_dir, cc)
      dir.create(subdir, showWarnings = FALSE, recursive = TRUE)
      saveRDS(pred, file = file.path(subdir, "singler_predictions.rds"), compress = "gzip")
      utils::write.csv(annot_tbl, file = file.path(subdir, "singler_annotation_table.csv"), row.names = FALSE)
      if (isTRUE(save_scores)) {
        utils::write.csv(pred$scores, file = file.path(subdir, "singler_scores_matrix.csv"), row.names = TRUE)
      }
    }
    
    results_list[[cc]] <- list(pred = pred, annotation_table = annot_tbl)
  }
  
  # [4/7] Concatenate all annotation tables (long format)
  .msg("[4/7] Concatenating annotation tables across clusterings…")
  all_annot <- do.call(rbind, lapply(results_list, `[[`, "annotation_table"))
  
  # [5/7] Optional global save
  if (!is.null(save_dir)) {
    utils::write.csv(all_annot, file = file.path(save_dir, "ALL_annotations_long.csv"), row.names = FALSE)
  }
  
  # [6/7] Done + timing
  end_time <- Sys.time()
  elapsed  <- as.numeric(difftime(end_time, start_time, units = "secs"))
  .msg("[5/7] Done. Clusterings processed: %d", length(cluster_cols))
  .msg("[6/7] Total time: %.2f seconds.", elapsed)
  
  # [7/7] Return results
  list(
    results_per_clustering = results_list, # named list: clustering -> list(pred, annotation_table)
    annotations_long       = all_annot,   # combined annotation table (long)
    n_common_genes         = length(common),
    n_genes_ref            = nrow(sce_ref),
    n_genes_query          = nrow(sce_query),
    cluster_cols           = cluster_cols
  )
}
