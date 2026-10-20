#' Find cluster-specific marker genes (strict specificity, with progress bar)
#'
#' @param raw_counts_list list of matrices (genes x cells), one per cluster
#' @param min_mean minimal mean expression in the target cluster
#' @param min_log2fc minimal log2 fold-change vs every other cluster
#' @param min_tau minimal specificity index (tau, 0-1)
#' @param verbose logical; show messages and progress bar
#'
#' @return data.frame with specific marker genes per cluster
find_specific_markers <- function(
    raw_counts_list,
    min_mean   = 0.1,
    min_log2fc = 1,
    min_tau    = 0.8,
    verbose    = TRUE
){
  stopifnot(is.list(raw_counts_list), length(raw_counts_list) > 1)
  
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  cluster_names <- names(raw_counts_list) %||% paste0("cluster_", seq_along(raw_counts_list)-1)
  
  # 1. mean expr per cluster
  cluster_means <- sapply(seq_along(raw_counts_list), function(i){
    Matrix::rowMeans(raw_counts_list[[i]])
  })
  colnames(cluster_means) <- cluster_names
  
  # 2. tau index (specificity)
  tau <- apply(cluster_means, 1, function(x){
    if (max(x) == 0) return(0)
    sum(1 - (x / max(x))) / (length(x) - 1)
  })
  
  clusts <- colnames(cluster_means)
  
  # 3. choose apply function (progress bar or not)
  if (verbose && requireNamespace("pbapply", quietly = TRUE)) {
    lapply_fun <- pbapply::pblapply
  } else {
    lapply_fun <- lapply
    if (verbose) message("Package 'pbapply' not found, progress bar disabled.")
  }
  
  # 4. iterate per cluster
  res_list <- lapply_fun(clusts, function(cl){
    if (verbose) message(sprintf("Processing %s ...", cl))
    
    expr_target <- cluster_means[, cl]
    expr_other  <- cluster_means[, setdiff(clusts, cl), drop=FALSE]
    
    # log2FC vs each other cluster
    logfc_mat <- log2((expr_target + 1e-6) / (expr_other + 1e-6))
    min_logfc <- apply(logfc_mat, 1, min)
    
    keep <- which(
      expr_target >= min_mean &
        min_logfc   >= min_log2fc &
        tau         >= min_tau
    )
    
    if (length(keep) > 0) {
      df <- data.frame(
        gene        = rownames(cluster_means)[keep],
        cluster     = cl,
        mean_target = expr_target[keep],
        min_log2fc  = min_logfc[keep],
        tau         = tau[keep],
        stringsAsFactors = FALSE
      )
      df[order(-df$min_log2fc, -df$mean_target), ]
    } else {
      data.frame(
        gene        = character(),
        cluster     = character(),
        mean_target = numeric(),
        min_log2fc  = numeric(),
        tau         = numeric(),
        stringsAsFactors = FALSE
      )
    }
  })
  
  out <- do.call(rbind, res_list)
  rownames(out) <- NULL
  return(out)
}
