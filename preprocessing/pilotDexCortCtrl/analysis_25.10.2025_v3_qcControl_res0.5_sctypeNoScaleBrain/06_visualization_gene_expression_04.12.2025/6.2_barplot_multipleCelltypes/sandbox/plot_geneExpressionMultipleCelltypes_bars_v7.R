plot_geneExpressionMultiCelltypesDexCortEdger_barplot <- function(
    gene_symbol,
    combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
    integrated_pilotDexCortCtrl_filtered,
    celltypes_order = c(
      "Radial glial cells",
      "Astrocytes",
      "Oligodendrocyte precursor cells",
      "Oligodendrocytes",
      "Immature neurons",
      "Mature neurons"
    )
) {
  
  ###############################################################
  ### A) PAKIETY
  ###############################################################
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(tibble)
  
  ###############################################################
  ### C) WYCIĄGNIĘCIE GENU ZE WSZYSTKICH CELL TYPES
  ###############################################################
  df_gene <- combined_summaryEdgerChi2_res0.3RNAnoScaledBrain %>%
    purrr::map(~ dplyr::filter(.x, gene_symbol == !!gene_symbol)) %>%
    .[celltypes_order] %>%           
    bind_rows(.id = "celltype") %>%
    mutate(
      celltype = factor(celltype, levels = celltypes_order)
    )
  
  if (nrow(df_gene) == 0) {
    stop(paste("Brak genu:", gene_symbol))
  }
  
  ###############################################################
  ### D) SAMPLE MAP
  ###############################################################
  sample_map <- integrated_pilotDexCortCtrl_filtered@meta.data %>%
    dplyr::select(bc1_well, treatment) %>%
    dplyr::distinct()
  
  ###############################################################
  ### E) LONG FORMAT (CPM dla A1–A12)
  ###############################################################
  df_long <- df_gene %>%
    dplyr::select(celltype, gene_symbol, starts_with("A")) %>%
    tidyr::pivot_longer(
      cols      = starts_with("A"),
      names_to  = "bc1_well",
      values_to = "CPM"
    ) %>%
    dplyr::left_join(sample_map, by = "bc1_well") %>%
    dplyr::mutate(
      treatment = factor(treatment, levels = c("CTRL", "CORT", "DEX"))
    )
  
  ###############################################################
  ### F) PODSTAWOWE PODSUMOWANIE
  ###############################################################
  df_summary_base <- df_long %>%
    dplyr::group_by(celltype, treatment) %>%
    dplyr::summarise(
      mean_CPM = mean(CPM, na.rm = TRUE),
      sd_CPM   = sd(CPM,   na.rm = TRUE),
      n        = sum(!is.na(CPM)),
      se_CPM   = sd_CPM / sqrt(n),
      .groups  = "drop"
    )
  
  ###############################################################
  ### G) STATYSTYKI edgeR — BEZ LICZENIA
  ###############################################################
  stats_df <- df_gene %>%
    dplyr::select(
      celltype,
      PValue_edgerCORT,
      PValue_edgerDEX,
      log2FC_edgerCORT,
      log2FC_edgerDEX
    ) %>%
    dplyr::distinct()
  
  df_summary <- df_summary_base %>%
    dplyr::left_join(stats_df, by = "celltype") %>%
    dplyr::mutate(
      pval = dplyr::case_when(
        treatment == "CORT" ~ PValue_edgerCORT,
        treatment == "DEX"  ~ PValue_edgerDEX,
        TRUE ~ NA_real_
      ),
      logfc = dplyr::case_when(
        treatment == "CORT" ~ log2FC_edgerCORT,
        treatment == "DEX"  ~ log2FC_edgerDEX,
        TRUE ~ NA_real_
      ),
      p_label = dplyr::case_when(
        is.na(pval) ~ "",
        pval < 1e-4 ~ paste0("p = ", format(pval, scientific = TRUE, digits = 2)),
        TRUE        ~ paste0("p = ", signif(pval, 2))
      ),
      label_fc = ifelse(
        is.na(logfc),
        p_label,
        paste0(p_label, " ; ", round(logfc, 2))
      )
    ) %>%
    dplyr::select(
      celltype, treatment,
      mean_CPM, sd_CPM, n, se_CPM,
      pval, logfc, p_label, label_fc
    )
  
  ###############################################################
  ### H) BRACKETY
  ###############################################################
  df_brackets <- df_summary %>%
    dplyr::group_by(celltype) %>%
    dplyr::summarise(
      y_max = max(mean_CPM + se_CPM, na.rm = TRUE),
      CORT_label = (label_fc[treatment == "CORT"] %>% unique() %>% .[1]),
      DEX_label  = (label_fc[treatment == "DEX"]  %>% unique() %>% .[1]),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      y_CORT = y_max * 1.12,
      y_DEX  = y_max * 1.25,
      x_ctrl = 1,
      x_cort = 2,
      x_dex  = 3
    )
  
  ###############################################################
  ### I) WYKRES
  ###############################################################
  p <- ggplot() +
    geom_jitter(
      data = df_long,
      aes(x = treatment, y = CPM, color = treatment),
      width = 0.15,
      size  = 2.5,
      alpha = 0.85
    ) +
    geom_col(
      data = df_summary,
      aes(x = treatment, y = mean_CPM, fill = treatment),
      width = 0.9,
      alpha = 0.75
    ) +
    geom_errorbar(
      data = df_summary,
      aes(x = treatment,
          ymin = mean_CPM - se_CPM,
          ymax = mean_CPM + se_CPM),
      width = 0.15,
      linewidth = 0.7
    ) +
    geom_segment(
      data = df_brackets,
      aes(x = x_ctrl, xend = x_cort, y = y_CORT, yend = y_CORT),
      linewidth = 0.7
    ) +
    geom_segment(
      data = df_brackets,
      aes(x = x_ctrl, xend = x_ctrl, y = y_CORT*0.98, yend = y_CORT),
      linewidth = 0.7
    ) +
    geom_segment(
      data = df_brackets,
      aes(x = x_cort, xend = x_cort, y = y_CORT*0.98, yend = y_CORT),
      linewidth = 0.7
    ) +
    geom_segment(
      data = df_brackets,
      aes(x = x_ctrl, xend = x_dex, y = y_DEX, yend = y_DEX),
      linewidth = 0.7
    ) +
    geom_segment(
      data = df_brackets,
      aes(x = x_ctrl, xend = x_ctrl, y = y_DEX*0.98, yend = y_DEX),
      linewidth = 0.7
    ) +
    geom_segment(
      data = df_brackets,
      aes(x = x_dex, xend = x_dex, y = y_DEX*0.98, yend = y_DEX),
      linewidth = 0.7
    ) +
    geom_text(
      data = df_brackets,
      aes(x = (x_ctrl + x_cort)/2, y = y_CORT + y_max*0.05, label = CORT_label),
      size = 4.2,
      fontface = "bold"
    ) +
    geom_text(
      data = df_brackets,
      aes(x = (x_ctrl + x_dex)/2, y = y_DEX + y_max*0.05, label = DEX_label),
      size = 4.2,
      fontface = "bold"
    ) +
    facet_wrap(~ celltype, ncol = 6) +
    scale_fill_manual(values = c("CTRL" = "grey40", "CORT" = "indianred2", "DEX" = "steelblue3")) +
    scale_color_manual(values = c("CTRL" = "grey40", "CORT" = "indianred2", "DEX" = "steelblue3")) +
    theme_classic(base_size = 15) +
    labs(
      title = paste0(gene_symbol, " — Expression across 6 cell types"),
      x = "Treatment",
      y = "CPM"
    ) +
    theme(
      legend.position = "none",
      strip.text = element_text(size = 15, face = "bold")
    )
  
  return(p)
  
}

plot_geneExpressionMultiCelltypesDexCortEdger_barplot(
  gene_symbol = "ZBTB16",
  combined_summaryEdgerChi2_res0.3RNAnoScaledBrain,
  integrated_pilotDexCortCtrl_filtered
)
