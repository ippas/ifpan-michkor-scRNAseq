`%||%` <- function(x, y) if (!is.null(x)) x else y

analyze_gene_detection_signed <- function(df, epsilon = 1e-6, return_long = FALSE) {
  results <- df
  
  # --- Przygotowanie pustych kolumn ---
  new_cols <- c(
    "p_chi2_CTRL_vs_CORT", "OR_CTRL_vs_CORT", "CramerV_CTRL_vs_CORT", "CramerV_signed_CTRL_vs_CORT",
    "p_chi2_CTRL_vs_DEX",  "OR_CTRL_vs_DEX",  "CramerV_CTRL_vs_DEX",  "CramerV_signed_CTRL_vs_DEX"
  )
  for (col in new_cols) results[[col]] <- NA_real_
  
  for (i in seq_len(nrow(df))) {
    # --- Dane wejściowe ---
    ctrl_detected <- (df$CTRL_detected[i] %||% 0) + epsilon
    ctrl_total    <- (df$CTRL_total[i] %||% 0) + epsilon
    cort_detected <- (df$CORT_detected[i] %||% 0) + epsilon
    cort_total    <- (df$CORT_total[i] %||% 0) + epsilon
    dex_detected  <- (df$DEX_detected[i] %||% 0) + epsilon
    dex_total     <- (df$DEX_total[i] %||% 0) + epsilon
    
    # --- Pomijanie niekompletnych wierszy ---
    if (any(is.na(c(ctrl_detected, ctrl_total, cort_detected, cort_total, dex_detected, dex_total)))) next
    
    ctrl_absent <- ctrl_total - ctrl_detected
    cort_absent <- cort_total - cort_detected
    dex_absent  <- dex_total - dex_detected
    
    # --- CTRL vs CORT ---
    tbl_cort <- matrix(c(ctrl_detected, ctrl_absent,
                         cort_detected, cort_absent),
                       nrow = 2, byrow = TRUE)
    chi_cort <- suppressWarnings(chisq.test(tbl_cort, correct = FALSE))
    OR_cort  <- (cort_detected / cort_absent) / (ctrl_detected / ctrl_absent)
    log2OR_cort <- log2(OR_cort)
    cramers_cort <- sqrt(chi_cort$statistic / (sum(tbl_cort) * (min(dim(tbl_cort)) - 1)))
    cramers_signed_cort <- as.numeric(cramers_cort) * ifelse(log2OR_cort > 0, 1, -1)
    
    # --- CTRL vs DEX ---
    tbl_dex <- matrix(c(ctrl_detected, ctrl_absent,
                        dex_detected, dex_absent),
                      nrow = 2, byrow = TRUE)
    chi_dex <- suppressWarnings(chisq.test(tbl_dex, correct = FALSE))
    OR_dex  <- (dex_detected / dex_absent) / (ctrl_detected / ctrl_absent)
    log2OR_dex <- log2(OR_dex)
    cramers_dex <- sqrt(chi_dex$statistic / (sum(tbl_dex) * (min(dim(tbl_dex)) - 1)))
    cramers_signed_dex <- as.numeric(cramers_dex) * ifelse(log2OR_dex > 0, 1, -1)
    
    # --- Zapis wyników ---
    results$p_chi2_CTRL_vs_CORT[i]  <- chi_cort$p.value
    results$OR_CTRL_vs_CORT[i]      <- OR_cort
    results$CramerV_CTRL_vs_CORT[i] <- as.numeric(cramers_cort)
    results$CramerV_signed_CTRL_vs_CORT[i] <- cramers_signed_cort
    
    results$p_chi2_CTRL_vs_DEX[i]   <- chi_dex$p.value
    results$OR_CTRL_vs_DEX[i]       <- OR_dex
    results$CramerV_CTRL_vs_DEX[i]  <- as.numeric(cramers_dex)
    results$CramerV_signed_CTRL_vs_DEX[i] <- cramers_signed_dex
  }
  
  # --- Format long ---
  if (return_long) {
    long_df <- results %>%
      dplyr::select(gene_symbol,
                    p_chi2_CTRL_vs_CORT, OR_CTRL_vs_CORT, CramerV_CTRL_vs_CORT, CramerV_signed_CTRL_vs_CORT,
                    p_chi2_CTRL_vs_DEX,  OR_CTRL_vs_DEX,  CramerV_CTRL_vs_DEX,  CramerV_signed_CTRL_vs_DEX) %>%
      tidyr::pivot_longer(
        cols = -gene_symbol,
        names_to = c(".value", "treatment"),
        names_pattern = "(.*)_CTRL_vs_(CORT|DEX)"
      ) %>%
      dplyr::mutate(
        log2_OR = log2(OR),
        direction = ifelse(log2_OR > 0, "up", "down")
      )
    
    return(long_df)
  }
  
  return(results)
}