#' Pairwise Fisher's exact tests between marker sets (Enrichr-style), with summary table
#'
#' Builds 2x2 tables for each pair of sets within a given universe and computes:
#' Fisher p-value (one-sided enrichment), odds ratio (OR), log2(OR).
#' Returns matrices and an aggregated, rowname-free pairs table.
compute_geneSet_fisher <- function(
    gene_lists,
    universe,
    min_size      = 1,
    alternative   = "greater",
    strict_checks = TRUE,
    smooth_or     = FALSE
) {
  stopifnot(is.list(gene_lists), length(gene_lists) > 0)
  stopifnot(is.character(universe), length(universe) > 0)
  
  # Clean background and intersect sets with universe
  U <- unique(universe)
  N <- length(U)
  sets <- lapply(gene_lists, function(v) {
    v <- unique(as.character(v[!is.na(v) & nzchar(v)]))
    intersect(v, U)
  })
  
  # Optional size filter (in-universe size)
  if (min_size > 1) sets <- sets[lengths(sets) >= min_size]
  if (!length(sets)) stop("No gene set remains after intersecting with universe and min_size filtering.")
  
  set_names <- names(sets)
  n_sets <- length(sets)
  sizes <- vapply(sets, length, integer(1))
  
  # Output containers
  overlap_mat <- matrix(0L, n_sets, n_sets, dimnames = list(set_names, set_names))
  or_mat      <- matrix(NA_real_, n_sets, n_sets, dimnames = list(set_names, set_names))
  log2or_mat  <- matrix(NA_real_, n_sets, n_sets, dimnames = list(set_names, set_names))
  pval_mat    <- matrix(NA_real_, n_sets, n_sets, dimnames = list(set_names, set_names))
  
  # Helper: compute OR (with optional 0.5 correction)
  compute_or <- function(a, b, c, d, smooth = FALSE) {
    if (smooth) { a <- a + 0.5; b <- b + 0.5; c <- c + 0.5; d <- d + 0.5 }
    OR <- (a * d) / (b * c)
    if (!is.finite(OR)) {
      if (a == 0 || d == 0) OR <- 0
      if (b == 0 || c == 0) OR <- Inf
    }
    OR
  }
  
  # Pairwise loop
  for (i in seq_len(n_sets)) {
    Ai <- sets[[i]]; ai <- sizes[i]
    for (j in i:n_sets) {
      Aj <- sets[[j]]; aj <- sizes[j]
      
      a <- length(intersect(Ai, Aj))         # overlap
      b <- ai - a                            # in A not in B
      c <- aj - a                            # in B not in A
      union_size <- length(union(Ai, Aj))
      d <- N - union_size                    # in neither
      
      # Fill symmetric overlap
      overlap_mat[i, j] <- overlap_mat[j, i] <- a
      
      # Sanity checks on the contingency structure
      if (strict_checks) {
        if (a + b != ai) stop(sprintf("Check failed: a+b!=|A| for %s vs %s", set_names[i], set_names[j]))
        if (a + c != aj) stop(sprintf("Check failed: a+c!=|B| for %s vs %s", set_names[i], set_names[j]))
        if (a + b + c + d != N) stop(sprintf("Check failed: a+b+c+d!=|U| for %s vs %s", set_names[i], set_names[j]))
      }
      if (min(a, b, c, d) < 0) next
      
      # Fisher exact (one-sided enrichment)
      ft <- suppressWarnings(fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = alternative))
      p  <- ft$p.value
      
      # Odds ratio (optionally smoothed)
      OR <- compute_or(a, b, c, d, smooth = smooth_or)
      log2OR <- if (is.finite(OR) && OR > 0) log2(OR) else if (OR == 0) -Inf else Inf
      
      pval_mat[i, j] <- pval_mat[j, i] <- p
      or_mat[i, j]   <- or_mat[j, i]   <- OR
      log2or_mat[i, j] <- log2or_mat[j, i] <- log2OR
    }
  }
  
  # BH adjust p-values on upper triangle and mirror
  upper_idx <- upper.tri(pval_mat, diag = TRUE)
  pvals <- pval_mat[upper_idx]
  qvals <- stats::p.adjust(pvals, method = "BH")
  qval_mat <- pval_mat
  qval_mat[upper_idx] <- qvals
  qval_mat[lower.tri(qval_mat)] <- t(qval_mat)[lower.tri(qval_mat)]
  
  # Build aggregated pair table (no stat filtering)
  pairs <- data.frame(
    set1 = character(0), set2 = character(0),
    overlap = integer(0),
    size_set1 = integer(0), size_set2 = integer(0),
    union_size = integer(0), sum_sizes = integer(0),
    pval = numeric(0), qval = numeric(0),
    or = numeric(0), log2or = numeric(0),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(n_sets)) {
    for (j in i:n_sets) {
      s1 <- set_names[i]; s2 <- set_names[j]
      a  <- overlap_mat[i, j]
      sz1 <- sizes[i]; sz2 <- sizes[j]
      u  <- length(union(sets[[i]], sets[[j]]))
      pairs <- rbind(pairs, data.frame(
        set1 = s1, set2 = s2,
        overlap = a,
        size_set1 = sz1, size_set2 = sz2,
        union_size = u,
        sum_sizes  = sz1 + sz2,
        pval = pval_mat[i, j],
        qval = qval_mat[i, j],
        or   = or_mat[i, j],
        log2or = log2or_mat[i, j],
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Sort & drop rownames
  ord <- order(pairs$qval, -pairs$overlap, -pairs$log2or, na.last = TRUE)
  pairs <- pairs[ord, , drop = FALSE]
  rownames(pairs) <- NULL
  
  list(
    overlap_mat = overlap_mat,
    or_mat      = or_mat,
    log2or_mat  = log2or_mat,
    pval_mat    = pval_mat,
    qval_mat    = qval_mat,
    sizes       = sizes,
    sets_used   = sets,
    pairs       = pairs
  )
}

