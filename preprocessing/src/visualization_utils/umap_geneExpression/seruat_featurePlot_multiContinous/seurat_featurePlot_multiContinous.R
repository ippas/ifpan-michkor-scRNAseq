seurat_featurePlot_multiContinuous <- function(
    seurat_obj,
    gene_patterns,
    assay = "RNA",
    slot = "counts",
    split.by = NULL,
    split_levels = NULL,
    cluster_col = NULL,
    reduction = "umap",
    cols = c("gray90", "firebrick"),
    point_size = 0.8,
    legend_position = "bottom",
    return_data = FALSE,
    title_size = 16,
    facet_title_size = 14,
    aggregation = c("sum", "mean"),
    summary_stat = c("sum", "mean")  # jak liczyć sumy per grupa
) {
  library(ggplot2)
  library(dplyr)
  library(Matrix)
  
  aggregation <- match.arg(aggregation)
  summary_stat <- match.arg(summary_stat)
  
  obj_tmp <- seurat_obj
  if (!is.null(assay)) DefaultAssay(obj_tmp) <- assay
  
  # --- znajdź geny ---
  all_features <- rownames(GetAssayData(obj_tmp[[assay]], slot = slot))
  matched <- unlist(lapply(gene_patterns, function(pat)
    grep(paste0("-", pat, "$"), all_features, value = TRUE)
  ))
  matched <- unique(matched)
  
  if (length(matched) == 0)
    stop("❌ No genes found matching provided patterns.")
  message("✅ Found ", length(matched), " matching genes: ", paste(matched, collapse = ", "))
  
  # --- embeddings ---
  emb <- as.data.frame(Embeddings(obj_tmp, reduction = reduction))
  colnames(emb) <- c("dim1", "dim2")[seq_len(ncol(emb))]
  
  # --- expression matrix ---
  expr_mat <- GetAssayData(obj_tmp[[assay]], slot = slot)[matched, , drop = FALSE]
  
  # --- agregacja ekspresji per komórka ---
  combined_expr <- switch(
    aggregation,
    "sum"  = Matrix::colSums(expr_mat),
    "mean" = Matrix::colMeans(expr_mat)
  )
  
  # --- zbuduj dane ---
  data_df <- emb %>% mutate(expression_combined = combined_expr)
  
  # --- sprawdź split.by ---
  if (!is.null(split.by)) {
    if (!split.by %in% colnames(obj_tmp@meta.data)) {
      stop(paste0("❌ Column '", split.by, "' not found in Seurat meta.data."))
    }
    data_df$treatment <- obj_tmp@meta.data[[split.by]]
  }
  
  # --- ustaw kolejność ---
  if (!is.null(split.by)) {
    current_values <- unique(data_df$treatment)
    if (is.null(split_levels)) {
      preferred_order <- c("CTRL", "CORT", "DEX")
      split_levels <- intersect(preferred_order, current_values)
      split_levels <- c(split_levels, setdiff(current_values, split_levels))
    }
    data_df$treatment <- factor(data_df$treatment, levels = split_levels)
  }
  
  # --- oblicz sumy ekspresji ---
  if (!is.null(split.by)) {
    summary_df <- data_df %>%
      group_by(treatment, .drop = FALSE) %>%
      summarise(
        total_expression = if (summary_stat == "sum") sum(expression_combined) else mean(expression_combined),
        .groups = "drop"
      )
    
    message("\n📊 Summary of combined expression per treatment group:\n")
    print(summary_df)
  } else {
    total_expr <- if (summary_stat == "sum") sum(data_df$expression_combined) else mean(data_df$expression_combined)
    message("\n📊 Total combined expression (all cells): ", signif(total_expr, 4))
    summary_df <- data.frame(total_expression = total_expr)
  }
  
  # --- zwróć dane, jeśli potrzebne ---
  if (return_data) return(list(data = data_df, summary = summary_df))
  
  # --- wykres ---
  p <- ggplot(data_df, aes(x = dim1, y = dim2, color = expression_combined)) +
    geom_point(size = point_size) +
    scale_color_gradientn(colors = cols, name = paste0("Expression (", aggregation, ")")) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = legend_position,
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(size = title_size, color = "black"),
      strip.text = element_text(size = facet_title_size, color = "black"),
      text = element_text(color = "black")
    ) +
    labs(
      title = paste0("Combined expression of ", length(matched), " genes (", aggregation, ")")
    )
  
  if (!is.null(split.by)) {
    p <- p + facet_wrap(~treatment)
  }
  
  return(p)
}



# seurat_featurePlot_multiContinuous(
#   seurat_obj = integrated_pilotDexCortCtrl_filtered,
#   gene_patterns = genes,
#   split.by = "treatment",
#   split_levels = c("CTRL", "CORT", "DEX"),
#   cluster_col = "res0.3_RNA_noScaled_brain",
#   cols = c("gray90", "firebrick"),
#   # detection_threshold = 0,
#   return_data = T,
#   legend_position = "bottom"
# ) -> tmp
# 
# 
# tmp$data %>% 
#   group_by(treatment) %>% 
#   nest() %>% 
#   mutate(sum_expression = map(data, ~ .x %>% .$expression_combined %>% sum)) %>% 
#   unnest(sum_expression)
