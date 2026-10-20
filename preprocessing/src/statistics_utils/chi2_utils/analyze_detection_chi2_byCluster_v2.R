analyze_detection_chi2_byCluster <- function(
    detection_list,
    return_long = FALSE,
    epsilon = 1e-6,
    verbose = TRUE,
    progress = TRUE
) {
  requireNamespace("dplyr")
  requireNamespace("tidyr")
  requireNamespace("purrr")
  
  analyze_df <- function(df, cluster_name) {
    
    ctrl_detected <- df$CTRL_detected + epsilon
    ctrl_absent   <- df$CTRL_total    - df$CTRL_detected + epsilon
    cort_detected <- df$CORT_detected + epsilon
    cort_absent   <- df$CORT_total    - df$CORT_detected + epsilon
    dex_detected  <- df$DEX_detected  + epsilon
    dex_absent    <- df$DEX_total     - df$DEX_detected  + epsilon
    
    tbl_ctrl_cort <- cbind(ctrl_detected, ctrl_absent, cort_detected, cort_absent)
    tbl_ctrl_dex  <- cbind(ctrl_detected, ctrl_absent, dex_detected, dex_absent)
    
    OR_cort  <- (cort_detected / cort_absent) / (ctrl_detected / ctrl_absent)
    OR_dex   <- (dex_detected / dex_absent) / (ctrl_detected / ctrl_absent)
    log2OR_cort <- log2(OR_cort)
    log2OR_dex  <- log2(OR_dex)
    
    chi2_vec <- function(x) suppressWarnings(chisq.test(matrix(x, nrow = 2, byrow = TRUE), correct = FALSE))
    chi2_cort <- suppressWarnings(purrr::map_dbl(seq_len(nrow(df)), ~ chi2_vec(tbl_ctrl_cort[.x, ])$p.value))
    chi2_dex <- suppressWarnings(purrr::map_dbl(seq_len(nrow(df)), ~ chi2_vec(tbl_ctrl_dex[.x, ])$p.value))
    
    cramers <- function(x) {
      chi <- chisq.test(matrix(x, nrow = 2, byrow = TRUE), correct = FALSE)
      sqrt(chi$statistic / (sum(x) * (min(dim(matrix(x,2,2))) - 1)))
    }
    
    CramerV_cort <- purrr::map_dbl(seq_len(nrow(df)), ~ as.numeric(cramers(tbl_ctrl_cort[.x, ])))
    CramerV_dex  <- purrr::map_dbl(seq_len(nrow(df)), ~ as.numeric(cramers(tbl_ctrl_dex[.x, ])))
    
    out <- dplyr::bind_cols(
      df,
      tibble::tibble(
        cluster = cluster_name,
        p_chi2_CTRL_vs_CORT = chi2_cort,
        OR_CTRL_vs_CORT = OR_cort,
        CramerV_CTRL_vs_CORT = CramerV_cort,
        CramerV_signed_CTRL_vs_CORT = CramerV_cort * ifelse(log2OR_cort > 0, 1, -1),
        p_chi2_CTRL_vs_DEX = chi2_dex,
        OR_CTRL_vs_DEX = OR_dex,
        CramerV_CTRL_vs_DEX = CramerV_dex,
        CramerV_signed_CTRL_vs_DEX = CramerV_dex * ifelse(log2OR_dex > 0, 1, -1)
      )
    )
    
    if (return_long) {
      out <- out %>%
        tidyr::pivot_longer(
          cols = c(p_chi2_CTRL_vs_CORT, OR_CTRL_vs_CORT, CramerV_CTRL_vs_CORT,
                   CramerV_signed_CTRL_vs_CORT, p_chi2_CTRL_vs_DEX, OR_CTRL_vs_DEX,
                   CramerV_CTRL_vs_DEX, CramerV_signed_CTRL_vs_DEX),
          names_to = c(".value", "comparison"),
          names_pattern = "(.*)_CTRL_vs_(CORT|DEX)"
        ) %>%
        dplyr::mutate(log2_OR = log2(OR))
    }
    return(out)
  }
  
  if (verbose) message("📦 Running chi² detection analysis for all clusters...")
  
  # LIST output with progress bar
  results_list <- vector("list", length(detection_list))
  names(results_list) <- names(detection_list)
  
  if (progress) pb <- utils::txtProgressBar(min = 0, max = length(detection_list), style = 3)
  
  i <- 0
  for (cl in names(detection_list)) {
    results_list[[cl]] <- analyze_df(detection_list[[cl]], cl)
    i <- i + 1
    if (progress) utils::setTxtProgressBar(pb, i)
  }
  
  if (progress) close(pb)
  
  return(results_list)
}
