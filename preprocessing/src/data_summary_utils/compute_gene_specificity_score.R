#' Compute specificity exactly like your combine_specificity_mult()
#'
#' - Rowwise max-normalization: x / max(x); zeros-row → 0s
#' - Tau: 1 - mean(x / max(x)); zeros-row → NA
#' - Combined: spec_max * tau (broadcast per row)
#'
#' Returns the same structure as your original function.
compute_gene_specificity_scores <- function(expr_mat) {
  stopifnot(is.matrix(expr_mat) || is.data.frame(expr_mat))
  X <- as.matrix(expr_mat)
  mode(X) <- "numeric"
  
  # spec_max
  spec_max <- t(apply(X, 1, function(x) {
    mx <- max(x, na.rm = TRUE)
    if (!is.finite(mx) || mx == 0) return(rep(0, length(x)))
    x / mx
  }))
  
  # tau
  spec_tau <- apply(X, 1, function(x) {
    mx <- max(x, na.rm = TRUE)
    if (!is.finite(mx) || mx == 0) return(NA_real_)
    1 - (sum(x / mx, na.rm = TRUE) / length(x))
  })
  
  # combined
  spec_combined <- spec_max * spec_tau
  
  # same return shape as yours
  list(max = spec_max, tau = spec_tau, combined = spec_combined)
}
