prepare_expression_reference_singler_v1 <- function(
    expression_input,                  # path or Seurat/SCE/matrix object
    meta_input = NULL,                 # optional: path to CSV or data.frame
    label_col = "Cell.Type",
    full_names_map = NULL,             # named char: short -> full
    assay_name = "RNA",
    log_normalize_if_missing = TRUE,
    save_rds = NULL,
    verbose = TRUE
) {
  start_time <- Sys.time()
  .msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))
  
  # 1) Load/build base object with counts
  if (is.character(expression_input) && file.exists(expression_input)) {
    if (grepl("\\.h5$", expression_input, ignore.case = TRUE)) {
      .msg("[1/7] Loading 10x H5 matrix: %s", basename(expression_input))
      counts_mat <- Read10X_h5(expression_input)
      seurat_ref <- CreateSeuratObject(counts = counts_mat, assay = assay_name)
    } else if (dir.exists(expression_input)) {
      .msg("[1/7] Loading 10x folder (MTX): %s", expression_input)
      counts_mat <- Read10X(expression_input)
      seurat_ref <- CreateSeuratObject(counts = counts_mat, assay = assay_name)
    } else {
      stop("Unsupported file path in 'expression_input' (expected .h5 or 10x folder).")
    }
  } else if (inherits(expression_input, "Seurat")) {
    .msg("[1/7] Using provided Seurat object.")
    seurat_ref <- expression_input
  } else if (inherits(expression_input, "SingleCellExperiment")) {
    .msg("[1/7] Using provided SingleCellExperiment object.")
    sce_ref <- expression_input
  } else if (inherits(expression_input, "matrix") || inherits(expression_input, "dgCMatrix")) {
    .msg("[1/7] Using provided matrix; creating Seurat object.")
    seurat_ref <- CreateSeuratObject(counts = expression_input, assay = assay_name)
  } else {
    stop("Unsupported 'expression_input' type. Provide path (.h5 or 10x folder), Seurat, SCE, or matrix.")
  }
  
  # 2) Metadata attach (optional), then convert to SCE
  if (exists("seurat_ref")) {
    if (!is.null(meta_input)) {
      .msg("[2/7] Attaching metadata…")
      if (is.character(meta_input) && file.exists(meta_input)) {
        meta_df <- utils::read.csv(meta_input)
      } else if (is.data.frame(meta_input)) {
        meta_df <- meta_input
      } else {
        stop("Unsupported 'meta_input'. Provide CSV path or data.frame.")
      }
      if (!"Barcode" %in% colnames(meta_df)) {
        stop("Metadata must contain a 'Barcode' column to align with cell barcodes.")
      }
      rownames(meta_df) <- meta_df$Barcode
      meta_df <- meta_df[colnames(seurat_ref), , drop = FALSE]
      seurat_ref <- AddMetaData(seurat_ref, meta_df)
      .msg("       Metadata attached: %d columns.", ncol(meta_df))
    }
    .msg("[3/7] Converting to SingleCellExperiment…")
    sce_ref <- as.SingleCellExperiment(seurat_ref)
  }
  
  # 3) Validate label column
  .msg("[4/7] Validating label column: '%s'…", label_col)
  if (!label_col %in% colnames(colData(sce_ref))) {
    stop(sprintf("Label column '%s' not found in reference colData.", label_col))
  }
  
  # 4) Optional mapping to full names
  if (!is.null(full_names_map)) {
    .msg("[5/7] Mapping labels to full names (creating 'Cell.Type.Full')…")
    sce_ref$Cell.Type.Full <- plyr::revalue(
      x = as.character(sce_ref[[label_col]]),
      replace = full_names_map,
      warn_missing = FALSE
    )
  } else {
    .msg("[5/7] No mapping provided; using '%s' as 'Cell.Type.Full'.", label_col)
    sce_ref$Cell.Type.Full <- as.character(sce_ref[[label_col]])
  }
  
  # 5) Ensure logcounts
  if (log_normalize_if_missing && !"logcounts" %in% assayNames(sce_ref)) {
    .msg("[6/7] 'logcounts' missing → running scuttle::logNormCounts()…")
    sce_ref <- scuttle::logNormCounts(sce_ref)
  } else {
    .msg("[6/7] 'logcounts' already present or normalization disabled.")
  }
  
  # 6) Save (optional)
  if (!is.null(save_rds)) {
    .msg("[7/7] Saving reference to RDS: %s", save_rds)
    saveRDS(sce_ref, file = save_rds, compress = "gzip")
  }
  
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "secs")
  .msg("Done. Reference dims: genes=%d, cells=%d; colData cols=%d.", 
       nrow(sce_ref), ncol(sce_ref), ncol(colData(sce_ref)))
  .msg("Total execution time: %.2f seconds.", as.numeric(elapsed))
  
  sce_ref
}
