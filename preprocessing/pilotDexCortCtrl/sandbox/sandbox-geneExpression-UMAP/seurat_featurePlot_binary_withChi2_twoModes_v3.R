seurat_featurePlot_binary <- function(
    seurat_obj,
    gene_pattern,
    assay = "RNA",
    slot = "counts",
    split.by = NULL,
    split_levels = NULL,
    cluster_col = NULL,
    reduction = "umap",
    cols = c("gray90", "firebrick"),
    point_size = 0.8,
    legend_position = "none",
    detection_threshold = 0,
    epsilon = 1e-6,
    add_table_on_plot = TRUE,
    add_OR = TRUE,
    return_data = FALSE,
    reference_level = "CTRL",
    mode = c("percent", "cells"),  # 🔹 nowy argument
    # 🔹 rozmiary tekstów (w px)
    font_size_stat_px = 16,
    font_size_table_px = 15,
    font_size_facet_title_px = 18,
    # 🔹 kolory i czcionki
    font_color_stat = "black",
    font_color_table = "black",
    font_color_facet_title = "black",
    font_family_table = "mono",
    font_family_facet_title = "sans"
) {
  library(ggplot2)
  library(dplyr)
  
  mode <- match.arg(mode)
  
  # ✅ przeliczenie px → pt dla geom_text
  px_to_pt <- function(px) px / 2.845
  font_size_stat_pt  <- px_to_pt(font_size_stat_px)
  font_size_table_pt <- px_to_pt(font_size_table_px)
  
  obj_tmp <- seurat_obj
  if (!is.null(assay)) DefaultAssay(obj_tmp) <- assay
  
  all_features <- rownames(GetAssayData(obj_tmp, assay = assay, slot = slot))
  matched <- grep(paste0("-", gene_pattern, "$"), all_features, value = TRUE)
  
  if (length(matched) == 0)
    stop(paste0("❌ No gene found matching pattern: ", gene_pattern))
  if (length(matched) > 1) {
    message("⚠️ Multiple matches found, using first: ", paste(matched, collapse = ", "))
    matched <- matched[1]
  } else {
    message("✅ Found feature: ", matched)
  }
  
  emb <- as.data.frame(Embeddings(obj_tmp, reduction = reduction))
  colnames(emb) <- c("dim1", "dim2")[seq_len(ncol(emb))]
  
  expr_values <- GetAssayData(obj_tmp, assay = assay, slot = slot)[matched, ]
  detected_binary <- as.numeric(expr_values > detection_threshold)
  
  data_df <- emb %>%
    mutate(expression = expr_values, detected = detected_binary)
  
  if (!is.null(split.by)) data_df$treatment <- obj_tmp@meta.data[[split.by]]
  if (!is.null(cluster_col) && cluster_col %in% colnames(obj_tmp@meta.data))
    data_df$cluster <- obj_tmp@meta.data[[cluster_col]]
  
  if (!is.null(split.by)) {
    current_values <- unique(data_df$treatment)
    if (is.null(split_levels)) {
      preferred_order <- c("CTRL", "CORT", "DEX")
      split_levels <- intersect(preferred_order, current_values)
      split_levels <- c(split_levels, setdiff(current_values, split_levels))
    }
    data_df$treatment <- factor(data_df$treatment, levels = split_levels)
  }
  
  # --- summary table ---
  summary_df <- as.data.frame(table(data_df$treatment))
  colnames(summary_df) <- c("treatment", "total_cells")
  det_summary <- as.data.frame(table(data_df$treatment[data_df$detected == 1]))
  colnames(det_summary) <- c("treatment", "detected_count")
  summary_df <- merge(summary_df, det_summary, by = "treatment", all.x = TRUE)
  summary_df$detected_count[is.na(summary_df$detected_count)] <- 0
  summary_df$percent <- round(100 * summary_df$detected_count / summary_df$total_cells, 2)
  
  # --- chi² + OR ---
  table_list <- list()
  if (!is.null(split.by) && reference_level %in% summary_df$treatment) {
    ref_row <- summary_df[summary_df$treatment == reference_level, ]
    
    for (trt in summary_df$treatment) {
      if (trt == reference_level) next
      trt_row <- summary_df[summary_df$treatment == trt, ]
      
      if (mode == "percent") {
        # --- użycie procentów ---
        a <- ref_row$percent; b <- 100 - ref_row$percent
        c <- trt_row$percent; d <- 100 - trt_row$percent
        scale_label <- "%"
      } else {
        # --- użycie liczby komórek ---
        a <- ref_row$detected_count; b <- ref_row$total_cells - ref_row$detected_count
        c <- trt_row$detected_count; d <- trt_row$total_cells - trt_row$detected_count
        scale_label <- ""
      }
      
      # --- tabela kontyngencji ---
      cont <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
      rownames(cont) <- c(reference_level, trt)
      colnames(cont) <- c("Detected+", "Detected-")
      
      # --- test chi² ---
      test <- suppressWarnings(chisq.test(cont, correct = FALSE))
      expected <- round(test$expected, 2)
      
      # --- OR w poprawnym kierunku (trt/ref) ---
      or_value <- ((c + epsilon) * (b + epsilon)) / ((a + epsilon) * (d + epsilon))
      
      # --- drukowanie w terminalu ---
      cat("\n=====================================\n")
      cat(sprintf("Chi² test (%s mode): %s vs %s\n", mode, reference_level, trt))
      cat("=====================================\n\n")
      cat(sprintf("%-10s %-20s %-20s\n", "", "Det.+: obs (exp)", "Det.-: obs (exp)"))
      for (i in 1:2) {
        cat(sprintf("%-10s %-20s %-20s\n",
                    rownames(cont)[i],
                    sprintf("%.2f%s (%.2f%s)", cont[i, 1], scale_label, expected[i, 1], scale_label),
                    sprintf("%.2f%s (%.2f%s)", cont[i, 2], scale_label, expected[i, 2], scale_label)))
      }
      cat("\n")
      cat(sprintf("χ² = %-7.3f | df = %d | OR(trt/CTRL) = %-6.3f | p = %s\n",
                  test$statistic, test$parameter, or_value,
                  formatC(test$p.value, format = "e", digits = 3)))
      
      # --- tekst do wykresu ---
      label_top <- if (add_OR) {
        paste0("χ²=", round(test$statistic, 2),
               "\nOR=", round(or_value, 3),
               "\np=", formatC(test$p.value, format = "e", digits = 2))
      } else {
        paste0("χ²=", round(test$statistic, 2),
               "\np=", formatC(test$p.value, format = "e", digits = 2))
      }
      
      table_text <- paste0(
        sprintf("%-10s %-16s %-16s\n", "", "Det.+: obs (exp)", "Det.-: obs (exp)"),
        paste0(
          sprintf("%-10s %-16s %-16s",
                  rownames(cont)[1],
                  sprintf("%.2f%s (%.2f%s)", cont[1,1], scale_label, expected[1,1], scale_label),
                  sprintf("%.2f%s (%.2f%s)", cont[1,2], scale_label, expected[1,2], scale_label)),
          "\n",
          sprintf("%-10s %-16s %-16s",
                  rownames(cont)[2],
                  sprintf("%.2f%s (%.2f%s)", cont[2,1], scale_label, expected[2,1], scale_label),
                  sprintf("%.2f%s (%.2f%s)", cont[2,2], scale_label, expected[2,2], scale_label))
        )
      )
      
      table_list[[trt]] <- list(
        treatment = trt,
        chi2 = round(test$statistic, 2),
        or = round(or_value, 3),
        p = formatC(test$p.value, format = "e", digits = 2),
        label_top = label_top,
        text = table_text
      )
    }
  }
  
  # --- facet labels ---
  data_df <- merge(data_df, summary_df, by = "treatment", all.x = TRUE)
  data_df$facet_label <- paste0(
    data_df$treatment, "\n",
    data_df$detected_count, "/", data_df$total_cells,
    " (", data_df$percent, "%)"
  )
  
  if (!is.null(split_levels)) {
    data_df$treatment <- factor(data_df$treatment, levels = split_levels)
    facet_labels_ordered <- summary_df %>%
      arrange(match(treatment, split_levels)) %>%
      mutate(facet_label = paste0(
        treatment, "\n",
        detected_count, "/", total_cells, " (", percent, "%)"
      )) %>%
      pull(facet_label)
    data_df$facet_label <- factor(data_df$facet_label, levels = facet_labels_ordered)
  }
  
  # --- plot base ---
  title_text <- paste0(matched, " (binary detection, ", mode, " mode)")
  p <- ggplot(data_df, aes(x = dim1, y = dim2, color = factor(detected))) +
    geom_point(size = point_size) +
    scale_color_manual(values = c("0" = cols[1], "1" = cols[2])) +
    theme_classic() +
    theme(
      legend.position = legend_position,
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      strip.text = element_text(
        size = font_size_facet_title_px,
        color = font_color_facet_title,
        family = font_family_facet_title
      ),
      plot.title = element_text(color = "black")
    ) +
    labs(title = title_text, color = "Detection") +
    facet_wrap(~facet_label, nrow = 1)
  
  # --- annotations ---
  if (length(table_list) > 0) {
    ann_stats <- do.call(rbind, lapply(names(table_list), function(trt) {
      data.frame(
        treatment = trt,
        facet_label = paste0(
          trt, "\n",
          summary_df$detected_count[summary_df$treatment == trt], "/",
          summary_df$total_cells[summary_df$treatment == trt],
          " (", summary_df$percent[summary_df$treatment == trt], "%)"
        ),
        label_top = table_list[[trt]]$label_top,
        label_bottom = table_list[[trt]]$text
      )
    }))
    
    x_max <- max(data_df$dim1, na.rm = TRUE)
    y_max <- max(data_df$dim2, na.rm = TRUE)
    x_min <- min(data_df$dim1, na.rm = TRUE)
    y_min <- min(data_df$dim2, na.rm = TRUE)
    
    p <- p + geom_text(
      data = ann_stats,
      aes(x = x_max, y = y_max, label = label_top),
      inherit.aes = FALSE,
      hjust = 1.1, vjust = 1.1,
      size = font_size_stat_pt,
      color = font_color_stat
    )
    
    if (add_table_on_plot) {
      p <- p + geom_text(
        data = ann_stats,
        aes(x = x_min, y = y_min, label = label_bottom),
        inherit.aes = FALSE,
        hjust = 0, vjust = -0.1,
        size = font_size_table_pt,
        family = font_family_table,
        color = font_color_table,
        lineheight = 1.1
      )
    }
  }
  
  return(p)
}



seurat_featurePlot_binary(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  gene_pattern = "ZBTB16",
  split.by = "treatment",
  split_levels = c("CTRL", "CORT", "DEX"),
  cluster_col = "res0.3_RNA_noScaled_brain",
  detection_threshold = 0,
  legend_position = "none",
  epsilon = 1e-6,
  reference_level = "CTRL",
  mode = "percent",
  add_OR = TRUE,
  add_table_on_plot = TRUE,
  font_size_stat_px = 12,
  font_size_table_px = 12,
  font_size_facet_title_px = 14
)
