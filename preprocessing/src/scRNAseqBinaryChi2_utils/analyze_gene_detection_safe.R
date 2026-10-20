analyze_gene_detection_safe <- function(df, epsilon = 1e-6) {
  results <- df
  
  # Przygotowanie pustych kolumn na wyniki
  results$p_chi2_CTRL_vs_CORT <- NA_real_
  results$OR_CTRL_vs_CORT <- NA_real_
  results$CramerV_CTRL_vs_CORT <- NA_real_
  
  results$p_chi2_CTRL_vs_DEX <- NA_real_
  results$OR_CTRL_vs_DEX <- NA_real_
  results$CramerV_CTRL_vs_DEX <- NA_real_
  
  for (i in seq_len(nrow(df))) {
    # dane wejściowe z epsilonem (dla bezpieczeństwa testu chi2)
    ctrl_detected <- (df$CTRL_detected[i] %||% 0) + epsilon
    ctrl_total    <- (df$CTRL_total[i] %||% 0) + epsilon
    cort_detected <- (df$CORT_detected[i] %||% 0) + epsilon
    cort_total    <- (df$CORT_total[i] %||% 0) + epsilon
    dex_detected  <- (df$DEX_detected[i] %||% 0) + epsilon
    dex_total     <- (df$DEX_total[i] %||% 0) + epsilon
    
    # pomijamy wiersze bez pełnych danych
    if (any(is.na(c(ctrl_detected, ctrl_total, cort_detected, cort_total, dex_detected, dex_total)))) next
    
    # brak ekspresji
    ctrl_absent <- ctrl_total - ctrl_detected
    cort_absent <- cort_total - cort_detected
    dex_absent  <- dex_total - dex_detected
    
    # --- CTRL vs CORT ---
    tbl_cort <- matrix(c(ctrl_detected, ctrl_absent,
                         cort_detected, cort_absent),
                       nrow = 2, byrow = TRUE)
    colnames(tbl_cort) <- c("Detected", "Absent")
    rownames(tbl_cort) <- c("CTRL", "CORT")
    
    chi_cort <- suppressWarnings(chisq.test(tbl_cort, correct = FALSE))
    OR_cort  <- (cort_detected / cort_absent) / (ctrl_detected / ctrl_absent)
    cramers_cort <- sqrt(chi_cort$statistic / (sum(tbl_cort) * (min(dim(tbl_cort)) - 1)))
    
    # --- CTRL vs DEX ---
    tbl_dex <- matrix(c(ctrl_detected, ctrl_absent,
                        dex_detected, dex_absent),
                      nrow = 2, byrow = TRUE)
    colnames(tbl_dex) <- c("Detected", "Absent")
    rownames(tbl_dex) <- c("CTRL", "DEX")
    
    chi_dex <- suppressWarnings(chisq.test(tbl_dex, correct = FALSE))
    OR_dex  <- (dex_detected / dex_absent) / (ctrl_detected / ctrl_absent)
    cramers_dex <- sqrt(chi_dex$statistic / (sum(tbl_dex) * (min(dim(tbl_dex)) - 1)))
    
    # --- zapis wyników ---
    results$p_chi2_CTRL_vs_CORT[i]  <- chi_cort$p.value
    results$OR_CTRL_vs_CORT[i]      <- OR_cort
    results$CramerV_CTRL_vs_CORT[i] <- as.numeric(cramers_cort)
    
    results$p_chi2_CTRL_vs_DEX[i]   <- chi_dex$p.value
    results$OR_CTRL_vs_DEX[i]       <- OR_dex
    results$CramerV_CTRL_vs_DEX[i]  <- as.numeric(cramers_dex)
  }
  
  return(results)
}