get_binaryExpression_multiGenes <- function(
    seurat_obj,
    gene_patterns,
    assay = "RNA",
    slot = "counts",
    split.by = "treatment",
    detection_threshold = 0
) {
  library(dplyr)
  
  DefaultAssay(seurat_obj) <- assay
  all_features <- rownames(GetAssayData(seurat_obj, assay = assay, slot = slot))
  matched <- unique(unlist(lapply(gene_patterns, function(pat)
    grep(paste0("-", pat, "$"), all_features, value = TRUE)
  )))
  if (length(matched) == 0) stop("No genes found matching provided patterns.")
  
  expr_mat <- GetAssayData(seurat_obj, assay = assay, slot = slot)[matched, , drop = FALSE]
  binary_mat <- (expr_mat > detection_threshold) * 1
  binary_mat <- as.matrix(binary_mat)
  genes <- rownames(binary_mat)
  
  meta <- seurat_obj@meta.data
  if (!split.by %in% colnames(meta)) stop("split.by column not found in meta.data")
  group_vec <- as.character(meta[[split.by]])
  groups <- unique(group_vec)
  
  summary_list <- lapply(groups, function(g) {
    cols <- which(group_vec == g)
    n_cells <- length(cols)
    n_detected_vec <- rowSums(binary_mat[, cols, drop = FALSE])
    data.frame(
      name_id = genes,
      treatment = g,
      n_cells = n_cells,
      n_detected = as.integer(n_detected_vec),
      percent_detected = round(100 * as.integer(n_detected_vec) / pmax(1L, n_cells), 2),
      stringsAsFactors = FALSE
    )
  })
  
  summary_df <- do.call(rbind, summary_list)
  
  summary_df <- summary_df %>%
    mutate(
      ensembl_id = sub("-.*", "", name_id),
      gene_symbol = sub(".*-", "", name_id)
    ) %>%
    select(ensembl_id, gene_symbol, name_id, treatment, n_cells, n_detected, percent_detected)
  
  list(
    raw_expr = expr_mat,
    binary_expr = binary_mat,
    summary_per_gene_treatment = summary_df
  )
}


# get_binaryExpression_multiGenes(
#   seurat_obj = integrated_pilotDexCortCtrl_filtered,
#   gene_patterns = genes
# ) -> tmp
# 
# 
# df_plot <- tmp$summary_per_gene_treatment %>%
#   filter(treatment %in% c("CTRL", "CORT")) %>% 
#   group_by(gene_symbol) %>% 
#   nest() %>% 
#   mutate(max_percent = map_dbl(data, ~ max(.x$percent_detected))) %>% 
#   unnest(cols = data) %>% 
#   ungroup() %>% 
#   filter(max_percent < 5)
# 
# # Upewniamy się, że treatment ma właściwą kolejność
# df_plot$treatment <- factor(df_plot$treatment, levels = c("CTRL", "CORT"))
# 
# ggplot(df_plot, aes(x = treatment, y = percent_detected, color = treatment, group = gene_symbol)) +
#   # geom_boxplot(fill = NA, width = 0.5, outlier.shape = NA, color = "black") +
#   geom_line(aes(group = gene_symbol),
#             color = "gray70", linewidth = 0.6, alpha = 0.6,
#             position = position_dodge(width = 0.15)) +
#   geom_jitter(width = 0.15, size = 2, alpha = 0.9) +
#   scale_color_manual(values = c("CTRL" = "gray40", "CORT" = "firebrick")) +
#   theme_minimal(base_size = 14) +
#   theme(
#     legend.position = "none",
#     axis.title.x = element_blank(),
#     axis.text.x = element_text(size = 14, face = "bold"),
#     axis.title.y = element_text(size = 14)
#   ) +
#   labs(y = "Percent detected per gene (<5%)")
# 
# tmp$summary_per_gene_treatment %>%
#   # dołącz kolumnę z różnicą (CORT - CTRL) policzoną per gen
#   left_join(
#     tmp$summary_per_gene_treatment %>%
#       filter(treatment %in% c("CTRL", "CORT")) %>%
#       select(gene_symbol, treatment, percent_detected) %>%
#       pivot_wider(names_from = treatment, values_from = percent_detected) %>%
#       mutate(diff_CORT_CTRL = CORT - CTRL) %>%
#       select(gene_symbol, diff_CORT_CTRL),
#     by = "gene_symbol"
#   ) %>% 
#   arrange(desc(diff_CORT_CTRL))
