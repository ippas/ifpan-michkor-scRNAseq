select_high_abundant_low_discrete <- function(
    raw_counts_list,
    top_n       = 50,
    min_pct     = 0.6,
    max_gini    = 0.55,
    min_mean    = 0,
    return_long = FALSE,
    sample_cells = NULL,
    verbose     = TRUE
){
  stopifnot(is.list(raw_counts_list), length(raw_counts_list) > 0)
  
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  # Gini helper
  .gini <- function(x){
    x <- as.numeric(x)
    if (!any(x > 0)) return(0)
    x <- sort(x)
    n <- length(x); sx <- sum(x)
    if (sx == 0) return(0)
    (2 * sum(x * seq_len(n)) / (n * sx)) - (n + 1) / n
  }
  
  # pbapply wrapper
  if (verbose && requireNamespace("pbapply", quietly = TRUE)) {
    gini_apply <- function(idx, mat) {
      pbapply::pbsapply(idx, function(r) .gini(mat[r, ]))
    }
  } else {
    gini_apply <- function(idx, mat) {
      sapply(idx, function(r) .gini(mat[r, ]))
    }
  }
  
  cluster_names <- names(raw_counts_list) %||% paste0("cluster_", seq_along(raw_counts_list)-1)
  
  out <- lapply(seq_along(raw_counts_list), function(i){
    cl_name <- cluster_names[i]
    message(sprintf("Processing %s ...", cl_name))
    
    mat <- raw_counts_list[[i]]
    
    # only convert to dense if it's not Matrix already
    if (!inherits(mat, "Matrix") && !is.matrix(mat)) {
      mat <- as.matrix(mat)
    }
    if (is.null(rownames(mat))) rownames(mat) <- paste0("gene_", seq_len(nrow(mat)))
    
    # optional subsample
    if (!is.null(sample_cells) && ncol(mat) > sample_cells) {
      set.seed(1L)
      cols <- sample(seq_len(ncol(mat)), sample_cells)
      mat <- mat[, cols, drop = FALSE]
    }
    
    n_cells <- ncol(mat)
    
    # basic stats (mean & pct ok na sparse)
    mean_expr   <- Matrix::rowMeans(mat)
    pct_expr    <- Matrix::rowSums(mat > 0) / n_cells
    
    # median/sd – uwaga: apply() = dense, więc tylko na subsamplowanych kolumnach
    if (n_cells > 2000) {
      set.seed(1L)
      subcols <- sample(seq_len(n_cells), 2000)
      median_expr <- apply(mat[, subcols, drop=FALSE], 1, stats::median)
      sd_expr     <- apply(mat[, subcols, drop=FALSE], 1, stats::sd)
    } else {
      median_expr <- apply(mat, 1, stats::median)
      sd_expr     <- apply(mat, 1, stats::sd)
    }
    
    # compute gini only for candidates
    pre_ok <- pct_expr >= min_pct & mean_expr >= min_mean
    gini_within <- rep(NA_real_, nrow(mat))
    idx <- which(pre_ok)
    if (length(idx) > 0) {
      gini_within[idx] <- gini_apply(idx, mat)
    }
    
    score <- log1p(mean_expr) * pct_expr * (1 - (gini_within %||% 1))
    
    df <- data.frame(
      gene        = rownames(mat),
      mean_expr   = mean_expr,
      median_expr = median_expr,
      sd_expr     = sd_expr,
      pct_expr    = pct_expr,
      gini_within = gini_within,
      score       = score,
      stringsAsFactors = FALSE
    )
    
    df <- subset(df, pct_expr >= min_pct & mean_expr >= min_mean & gini_within <= max_gini)
    
    if (nrow(df) == 0) {
      return(data.frame(
        gene        = character(),
        mean_expr   = numeric(),
        median_expr = numeric(),
        sd_expr     = numeric(),
        pct_expr    = numeric(),
        gini_within = numeric(),
        score       = numeric(),
        cluster     = character(),  # teraz długość 0, nie 1
        stringsAsFactors = FALSE
      ))
    }
    
    df <- df[order(-df$score, -df$mean_expr), , drop = FALSE]
    if (!is.null(top_n) && is.finite(top_n)) df <- utils::head(df, top_n)
    
    df$cluster <- cl_name
    rownames(df) <- NULL
    df
  })
  
  names(out) <- cluster_names
  if (isTRUE(return_long)) do.call(rbind, out) else out
}
