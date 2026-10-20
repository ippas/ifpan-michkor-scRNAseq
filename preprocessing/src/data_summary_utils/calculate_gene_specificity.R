calculate_gene_specificity <- function(expr_mat, na_omit = TRUE) {
  # Max-normalization
  spec_max <- t(apply(expr_mat, 1, function(x) {
    if (max(x) == 0) return(rep(0, length(x)))
    x / max(x)
  }))
  
  # Tau specificity
  spec_tau <- apply(expr_mat, 1, function(x) {
    if (max(x) == 0) return(NA)
    1 - (sum(x / max(x)) / length(x))
  })
  
  # Combined = spec_max × tau
  spec_combined <- spec_max * spec_tau
  
  # --- automatic assignment of column and row names ---
  if (!is.null(colnames(expr_mat))) {
    colnames(spec_max)      <- colnames(expr_mat)
    colnames(spec_combined) <- colnames(expr_mat)
  }
  rownames(spec_max)      <- rownames(expr_mat)
  rownames(spec_combined) <- rownames(expr_mat)
  
  # Optional removal of rows with NA in combined
  if (na_omit) {
    spec_combined <- na.omit(spec_combined)
  }
  
  list(max = spec_max, tau = spec_tau, combined = spec_combined)
}