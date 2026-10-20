plot_logfc_scatter <- function(
    data,
    threshold_x = 1,
    threshold_y = 1,
    color_main = "black",
    color_gray = "gray70",
    line_color = "gray30",
    corr_line_color = "firebrick",
    title = "Porównanie efektów DEX i CORT",
    fixed_axes = TRUE,
    label_min_x = NULL,
    label_max_x = NULL,
    label_min_y = NULL,
    label_max_y = NULL
) {
  library(ggrepel)
  library(dplyr)
  
  # 🔹 Oznaczenie koloru punktów
  data <- data %>%
    mutate(
      region_color = ifelse(
        abs(logFC.treatmentCORT) < threshold_x & abs(logFC.treatmentDEX) < threshold_y,
        color_gray,
        color_main
      )
    )
  
  # 🔹 Korelacja dla całego zbioru
  r_val <- cor(
    data$logFC.treatmentCORT,
    data$logFC.treatmentDEX,
    use = "complete.obs",
    method = "pearson"
  )
  
  # 🔹 Parametry linii korelacji
  mean_x <- mean(data$logFC.treatmentCORT, na.rm = TRUE)
  mean_y <- mean(data$logFC.treatmentDEX, na.rm = TRUE)
  sd_x   <- sd(data$logFC.treatmentCORT, na.rm = TRUE)
  sd_y   <- sd(data$logFC.treatmentDEX, na.rm = TRUE)
  slope  <- r_val * (sd_y / sd_x)
  intercept <- mean_y - slope * mean_x
  
  title_full <- sprintf("%s (r = %.3f)", title, r_val)
  
  # 🔹 Bezpieczne filtrowanie punktów do etykiet
  data_labeled <- data
  if (!is.null(label_min_x))
    data_labeled <- data_labeled %>% filter(logFC.treatmentCORT < label_min_x)
  if (!is.null(label_max_x))
    data_labeled <- data_labeled %>% bind_rows(data %>% filter(logFC.treatmentCORT > label_max_x))
  if (!is.null(label_min_y))
    data_labeled <- data_labeled %>% bind_rows(data %>% filter(logFC.treatmentDEX < label_min_y))
  if (!is.null(label_max_y))
    data_labeled <- data_labeled %>% bind_rows(data %>% filter(logFC.treatmentDEX > label_max_y))
  
  # Usunięcie duplikatów (gdy gen spełnia kilka warunków)
  data_labeled <- data_labeled %>% distinct(name_id, .keep_all = TRUE)
  
  p <- ggplot(data, aes(x = logFC.treatmentCORT, y = logFC.treatmentDEX)) +
    # Szary prostokąt
    annotate(
      "rect",
      xmin = -threshold_x, xmax = threshold_x,
      ymin = -threshold_y, ymax = threshold_y,
      fill = "gray85", alpha = 0.5
    ) +
    # Linie ±threshold
    geom_hline(yintercept = c(-threshold_y, threshold_y),
               color = line_color, linetype = "dashed", linewidth = 0.7) +
    geom_vline(xintercept = c(-threshold_x, threshold_x),
               color = line_color, linetype = "dashed", linewidth = 0.7) +
    # Punkty
    geom_point(aes(color = region_color), size = 1.8, alpha = 0.7, show.legend = FALSE) +
    scale_color_identity() +
    # 🔹 Linia korelacji (dla całego zbioru)
    geom_abline(intercept = intercept, slope = slope,
                color = corr_line_color, linewidth = 1) +
    # 🔹 Etykiety genów
    geom_text_repel(
      data = data_labeled,
      aes(label = gene_symbol),
      size = 3.5,
      color = "black",
      box.padding = 0.4,
      point.padding = 0.2,
      max.overlaps = Inf,
      segment.color = "gray60"
    ) +
    theme_classic(base_size = 14) +
    labs(
      x = "logFC (CORT vs CTRL)",
      y = "logFC (DEX vs CTRL)",
      title = title_full
    )
  
  if (isTRUE(fixed_axes)) p <- p + coord_equal() else p <- p + coord_cartesian()
  return(p)
}



plot_logfc_scatter(
  edgeR_anovaOneFactor_treatment_pseudobulk_res0.5_clusters13_results$cluster_Astrocytes$f_test_global,
  threshold_x = 1,
  threshold_y = 1,
  title = "Astrocytes: DEX vs CORT",
  fixed_axes = F
)

plot_logfc_scatter(
  edgeR_anovaOneFactor_treatment_pseudobulk_res0.5_clusters13_results$`cluster_Glutamatergic neurons`$f_test_global,
  threshold_x = 1,
  threshold_y = 1,
  label_max_x = 5,  # geny z dużym logFC(CORT)
  label_max_y = 2,  # geny z dużym logFC(DEX)
  label_min_x = -4,
  label_min_y = -4,
  title = "Glutaminergic neurons: DEX vs CORT",
  fixed_axes = FALSE
)

