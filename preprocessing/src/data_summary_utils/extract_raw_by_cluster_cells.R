#' Extract per-cell expression matrices split by cluster
#'
#' Returns a named list of matrices (or sparse matrices), one per cluster.
#' Columns are single cells; rows are genes. No aggregation is performed.
#'
#' @param seurat_obj  Seurat object.
#' @param cluster_col Column in meta.data with cluster IDs (default "seurat_clusters").
#' @param sample_col  Optional column in meta.data with sample IDs (used only for ordering/filtering).
#' @param assay       Assay to use (default "RNA").
#' @param layer       Seurat v5 layer (e.g., "counts", "data"; default "counts").
#' @param include_samples Optional character vector of samples to keep (if NULL, keep all).
#' @param order_by_sample Logical; if TRUE and sample_col provided, columns inside each
#'                        cluster are ordered by sample, then by cell name (default TRUE).
#' @param progress    Logical; show text progress bar (default TRUE).
#' @return A named list with elements `cluster_<id>`; each is a genes×cells matrix.
#'         If a cluster has 0 cells after filtering, it is skipped.
#'
#' @examples
#' mats <- extract_raw_by_cluster_cells(seurat_obj, cluster_col="seurat_clusters",
#'                                      sample_col="sample_id", assay="RNA", layer="counts")
extract_raw_by_cluster_cells <- function(
    seurat_obj,
    cluster_col      = "seurat_clusters",
    sample_col       = NULL,
    assay            = "RNA",
    layer            = "counts",
    include_samples  = NULL,
    order_by_sample  = TRUE,
    progress         = TRUE
) {
  # --- 1) Pull expression layer (kept as (sparse) matrix) ---
  expr <- Seurat::GetAssayData(seurat_obj, assay = assay, layer = layer)
  if (!is.matrix(expr) && !inherits(expr, "Matrix")) expr <- as.matrix(expr)
  
  # --- 2) Meta alignment & sanity checks ---
  meta <- seurat_obj@meta.data
  if (!(cluster_col %in% colnames(meta))) {
    stop(sprintf("Column '%s' not found in meta.data", cluster_col))
  }
  if (!is.null(sample_col) && !(sample_col %in% colnames(meta))) {
    stop(sprintf("sample_col '%s' not found in meta.data", sample_col))
  }
  
  cell_ids <- colnames(expr)
  # Ensure meta is in the same order as expr columns
  meta <- meta[cell_ids, , drop = FALSE]
  
  # Optional filtering by sample(s)
  if (!is.null(include_samples) && !is.null(sample_col)) {
    keep <- meta[[sample_col]] %in% include_samples
    expr <- expr[, keep, drop = FALSE]
    meta <- meta[keep, , drop = FALSE]
    cell_ids <- colnames(expr)
  }
  
  # Factorize clusters to preserve natural order
  clusters <- as.factor(meta[[cluster_col]])
  cl_levels <- levels(clusters)
  
  # --- 3) Split columns by cluster ---
  if (progress) pb <- utils::txtProgressBar(min = 0, max = length(cl_levels), style = 3)
  
  out_list <- list()
  step <- 0L
  
  for (cl in cl_levels) {
    step <- step + 1L
    idx <- which(clusters == cl)
    if (length(idx) == 0L) {
      if (progress) utils::setTxtProgressBar(pb, step)
      next
    }
    
    cols <- cell_ids[idx]
    
    # Optional deterministic ordering: by sample then by cell name
    if (!is.null(sample_col) && isTRUE(order_by_sample)) {
      ord <- order(meta[cols, sample_col, drop = TRUE], cols)
      cols <- cols[ord]
    }
    
    mat_cl <- expr[, cols, drop = FALSE]
    out_list[[paste0("cluster_", cl)]] <- mat_cl
    
    if (progress) utils::setTxtProgressBar(pb, step)
  }
  if (progress) close(pb)
  
  # Drop empty clusters (if any filtering made them empty)
  out_list <- out_list[vapply(out_list, ncol, integer(1)) > 0L]
  
  return(out_list)
}
