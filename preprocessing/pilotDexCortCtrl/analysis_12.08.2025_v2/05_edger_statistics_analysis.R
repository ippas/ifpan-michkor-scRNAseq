integrated_pilotDexCortCtrl <- seurat_merge_all_counts_layers(
  seurat_obj = integrated_pilotDexCortCtrl,
  assay = "RNA",
  remove_old = FALSE
)

agg_list <- aggregate_counts_by_cluster_sample(
  seurat_obj  = integrated_pilotDexCortCtrl,
  cluster_col = "seurat_clusters_res0.04",
  sample_col  = "sample_id",
  assay       = "RNA",
  layer       = "counts",  # <- Seurat v5
  agg         = "sum",
  return      = "list",
  progress    = TRUE
)

sample_map <- integrated_pilotDexCortCtrl@meta.data %>% distinct(sample_id, treatment)

# Wszystkie klastry, z komunikatami i paskiem postępu, bez gatekeep:
res <- edger_pseudobulk_anova_per_cluster(
  agg_list,
  sample_meta  = sample_map,
  group_levels = c("CTRL","DEX","CORT"),
  verbose      = TRUE,
  progress     = TRUE,
  gatekeep     = F
)



res$cluster_0$f_test_global %>% filter(FDR < 0.01) %>% 
  filter(logFC.treatmentDEX > 1 | logFC.treatmentCORT > 1)

res$cluster_0$f_test_global %>% head

res %>%
  keep(~ "f_test_global" %in% names(.x)) %>%                
  imap(~ .x$f_test_global %>%
         filter(FDR < 0.05 & (logFC.treatmentDEX > 1 | logFC.treatmentCORT > 1))) %>% 
  write_cluster_list_to_xlsx(file_path = "data/pilotDexCortCtrl/r_outputs/statistics_results/edger_analysis_15.08.2025_v1/edgeR_AllTreatments_res0.04_FDR0.01log2FC1.xlsx")


res$cluster_0$CORT_vs_DEX %>% filter(gene_symbol == "FKBP5")

res %>% 
  keep(~ "CORT_vs_CTRL" %in% names(.x)) %>%                
  imap(~ .x$CORT_vs_CTRL %>%
         filter(FDR < 0.05) %>% 
         filter(logFC > 1 | logFC < -1)) %>% 
  write_cluster_list_to_xlsx(file_path = "data/pilotDexCortCtrl/r_outputs/statistics_results/edger_analysis_15.08.2025_v1/edgeR_CTRLvsCORT_res0.04_FDR0.01log2FC1.xlsx")

res %>% 
  keep(~ "DEX_vs_CTRL" %in% names(.x)) %>%                
  imap(~ .x$DEX_vs_CTRL %>%
         filter(FDR < 0.05) %>% 
         filter(logFC > 1 | logFC < -1)) %>% 
  write_cluster_list_to_xlsx(file_path = "data/pilotDexCortCtrl/r_outputs/statistics_results/edger_analysis_15.08.2025_v1/edgeR_CTRLvsDEX_res0.04_FDR0.01log2FC1.xlsx")

res %>% 
  keep(~ "CORT_vs_DEX" %in% names(.x)) %>%                
  imap(~ .x$CORT_vs_DEX %>%
         filter(FDR < 0.05) %>% 
         filter(logFC > 1 | logFC < -1)) %>% 
  write_cluster_list_to_xlsx(file_path = "data/pilotDexCortCtrl/r_outputs/statistics_results/edger_analysis_15.08.2025_v1/edgeR_CORTvsDEX_res0.04_FDR0.01log2FC1.xlsx")

# ##############################################################################
# ---- for res 0.13 ----
# ##############################################################################
agg_list <- aggregate_counts_by_cluster_sample(
  seurat_obj  = integrated_pilotDexCortCtrl,
  cluster_col = "seurat_clusters_res0.13",
  sample_col  = "sample_id",
  assay       = "RNA",
  layer       = "counts",  # <- Seurat v5
  agg         = "sum",
  return      = "list",
  progress    = TRUE
)

sample_map <- integrated_pilotDexCortCtrl@meta.data %>% distinct(sample_id, treatment)

# Wszystkie klastry, z komunikatami i paskiem postępu, bez gatekeep:
res <- edger_pseudobulk_anova_per_cluster(
  agg_list,
  sample_meta  = sample_map,
  group_levels = c("CTRL","DEX","CORT"),
  verbose      = TRUE,
  progress     = TRUE,
  gatekeep     = FALSE
)

res$cluster_0$f_test_global %>% filter(FDR < 0.01) %>% 
  filter(logFC.treatmentDEX > 1 | logFC.treatmentCORT > 1)


res$cluster_2$combined %>% 
  filter(g_FDR < 0.05) %>% 
  filter( DEX_CTRL_logFC > 1 |  CORT_CTRL_logFC > 1 |  CORT_CTRL_logFC < -1 | DEX_CTRL_logFC < -1) 

res %>%
  keep(~ "f_test_global" %in% names(.x)) %>%                
  imap(~ .x$f_test_global %>%
         filter(FDR < 0.05 & (logFC.treatmentDEX > 1 | logFC.treatmentCORT > 1))) %>% 
  write_cluster_list_to_xlsx(file_path = "data/pilotDexCortCtrl/r_outputs/statistics_results/edger_analysis_15.08.2025_v1/edgeR_AllTreatments_res0.13_FDR0.01log2FC1.xlsx")


res$cluster_0$CORT_vs_DEX %>% filter(gene_symbol == "FKBP5")

res %>% 
  keep(~ "CORT_vs_CTRL" %in% names(.x)) %>%                
  imap(~ .x$CORT_vs_CTRL %>%
         filter(FDR < 0.05) %>% 
         filter(logFC > 1 | logFC < -1)) %>% 
  write_cluster_list_to_xlsx(file_path = "data/pilotDexCortCtrl/r_outputs/statistics_results/edger_analysis_15.08.2025_v1/edgeR_CTRLvsCORT_res0.13_FDR0.01log2FC1.xlsx")

res %>% 
  keep(~ "DEX_vs_CTRL" %in% names(.x)) %>%                
  imap(~ .x$DEX_vs_CTRL %>%
         filter(FDR < 0.05) %>% 
         filter(logFC > 1 | logFC < -1)) %>% 
  write_cluster_list_to_xlsx(file_path = "data/pilotDexCortCtrl/r_outputs/statistics_results/edger_analysis_15.08.2025_v1/edgeR_CTRLvsDEX_res0.13_FDR0.01log2FC1.xlsx")

res %>% 
  keep(~ "CORT_vs_DEX" %in% names(.x)) %>%                
  imap(~ .x$CORT_vs_DEX %>%
         filter(FDR < 0.05) %>% 
         filter(logFC > 1 | logFC < -1)) %>% 
  write_cluster_list_to_xlsx(file_path = "data/pilotDexCortCtrl/r_outputs/statistics_results/edger_analysis_15.08.2025_v1/edgeR_CORTvsDEX_res0.13_FDR0.01log2FC1.xlsx")


res_global <- edger_pseudobulk_anova_global(
  agg          = agg_list,       # or a single matrix genes x samples
  sample_meta  = sample_map,
  group_col    = "treatment",
  sample_col   = "sample_id",
  group_levels = c("CTRL","DEX","CORT"),
  robust       = TRUE,
  verbose      = TRUE
)

res_global$combined

res_global$f_test_global %>% 
  filter(FDR < 0.05) %>%
  filter(logFC.treatmentCORT > 1 | logFC.treatmentCORT < -1 | logFC.treatmentDEX > 1 | logFC.treatmentDEX < -1) %>%
  .$name_id %>% unique() %>% length()
  

res_global$f_test_global


# ---- correlation calculateion ----
tab <- res_global$f_test_global

# filter only genes that passed filtering (optional, e.g. FDR < 0.05)
sig <- subset(tab, FDR <= 0.05)




dotplot_log2FC_correlation_v2(
  res$cluster_2$f_test_global,
  FDR_threshold   = 0.05,
  thresholds      = c(-1, 1),      # linie przerywane
  threshold_x_text = 2.5,      # etykiety, jeśli |log2FC.x| >= 1.5
  threshold_y_text = 4.5,         # etykiety, jeśli |log2FC.y| >= 2
  axis_text_size = 12,
  axis_title_size = 12,
  title_size = 16,
  legend_text_size = 12,
  color_sig_points = "firebrick",
  # color_sig_line   = "firebrick",
  # color_sig_labels = "black",
  label_size = 12,
  # x_limits = c(-4, 10.5),
  # y_limits = c(-4, 10.5)
) -> p1

p1

dotplot_log2FC_correlation_v2(
  res$cluster_5$f_test_global,
  FDR_threshold   = 0.05,
  thresholds      = c(-1, 1),      # linie przerywane
  threshold_x_text = 3,      # etykiety, jeśli |log2FC.x| >= 1.5
  threshold_y_text = 3,         # etykiety, jeśli |log2FC.y| >= 2
  axis_text_size = 12,
  axis_title_size = 12,
  title_size = 16,
  legend_text_size = 12,
  color_sig_points = "firebrick",
  # color_sig_line   = "firebrick",
  # color_sig_labels = "black",
  label_size = 12,
  # x_limits = c(-7, 10.5),
  # y_limits = c(-6, 10.5),
  title_prefix = "cluster 5"
) -> p2

dotplot_log2FC_correlation_v2(
  res$cluster_2$f_test_global,
  col_x = "logFC.treatmentCORT",
  col_y = "logFC.treatmentDEX",
  FDR_threshold   = 0.05,
  thresholds      = c(-1, 1),      # linie przerywane
  threshold_x_text = 2.5,      # etykiety, jeśli |log2FC.x| >= 1.5
  threshold_y_text = 5,         # etykiety, jeśli |log2FC.y| >= 2
  axis_text_size = 12,
  axis_title_size = 12,
  title_size = 16,
  legend_text_size = 12,
  color_sig_points = "firebrick",
  # color_sig_line   = "firebrick",
  # color_sig_labels = "black",
  label_size = 12,
  # x_limits = c(-7, 10.5),
  # y_limits = c(-6, 10.5),
  title_prefix = "cluster 2"
) -> p1

dev.off()
svg("results/tmp/figure2poster_v2.svg", width = 10, height = 7.5)
p1$plot
dev.off()


res$cluster_2$f_test_global %>% 
  filter(FDR < 0.05, logFC.treatmentCORT > 1) %>% 
  filter(grepl("COL", gene_symbol))


p2$plot  + p1$plot



p1$plot +  p3$plot

dotplot_log2FC_correlation_wrap <- function(
    res,
    df_accessor = function(x) x$f_test_global,
    cluster_prefix = "cluster_",
    ncol = 3,
    ...
){
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Pakiet 'patchwork' jest wymagany do połączenia wykresów.")
  }
  
  clusters <- names(res)
  clusters <- clusters[grepl(paste0("^", cluster_prefix), clusters)]
  if (length(clusters) == 0) stop("Nie znaleziono klastrów z prefiksem ", cluster_prefix)
  
  plots <- lapply(clusters, function(cl){
    df <- df_accessor(res[[cl]])
    res_one <- dotplot_log2FC_correlation(df, ...)
    
    # policz n_sig
    dots <- list(...)
    FDR_threshold <- if ("FDR_threshold" %in% names(dots)) dots$FDR_threshold else NULL
    fdr_col       <- if ("fdr_col"       %in% names(dots)) dots$fdr_col       else "FDR"
    thresholds    <- if ("thresholds"    %in% names(dots)) dots$thresholds    else c(-1,1)
    col_x         <- if ("col_x"         %in% names(dots)) dots$col_x         else "logFC.treatmentDEX"
    col_y         <- if ("col_y"         %in% names(dots)) dots$col_y         else "logFC.treatmentCORT"
    
    n_sig <- NA_integer_
    if (!is.null(FDR_threshold) && fdr_col %in% colnames(df)) {
      sig_idx <- df[[fdr_col]] <= FDR_threshold &
        (abs(df[[col_x]]) >= max(thresholds) | abs(df[[col_y]]) >= max(thresholds))
      n_sig <- sum(sig_idx, na.rm = TRUE)
    }
    
    # nowy tytuł
    new_title <- sprintf("%s | cor=%.2f | n_sig=%d", cl, res_one$cor_all, n_sig)
    
    res_one$plot + ggplot2::ggtitle(new_title)
  })
  
  patchwork::wrap_plots(plots, ncol = ncol, guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
}


dotplot_log2FC_correlation_wrap(
  res,
  FDR_threshold    = 0.05,
  thresholds       = c(-1, 1),
  threshold_x_text = 2.5,
  threshold_y_text = 3,
  axis_text_size   = 12,
  axis_title_size  = 12,
  title_size       = 16,
  legend_text_size = 12,
  color_sig_points = "firebrick",
  color_sig_line   = "firebrick",
  color_sig_labels = "black",
  label_size       = 12
)


out

# ##############################################################################
# ---- for res 0.2 ----
# ##############################################################################
agg_list <- aggregate_counts_by_cluster_sample(
  seurat_obj  = integrated_pilotDexCortCtrl,
  cluster_col = "seurat_clusters_res0.2",
  sample_col  = "sample_id",
  assay       = "RNA",
  layer       = "counts",  # <- Seurat v5
  agg         = "sum",
  return      = "list",
  progress    = TRUE
)

sample_map <- integrated_pilotDexCortCtrl@meta.data %>% distinct(sample_id, treatment)

# Wszystkie klastry, z komunikatami i paskiem postępu, bez gatekeep:
res <- edger_pseudobulk_anova_per_cluster(
  agg_list,
  sample_meta  = sample_map,
  group_levels = c("CTRL","DEX","CORT"),
  verbose      = TRUE,
  progress     = TRUE,
  gatekeep     = FALSE
)

res$cluster_2$f_test_global %>% filter(FDR < 0.01) %>% 
  filter(logFC.treatmentDEX > 1 | logFC.treatmentCORT > 1) %>% 
  .$gene_symbol %>% cat(sep = "\n")


res %>%
  keep(~ "f_test_global" %in% names(.x)) %>%                
  imap(~ .x$f_test_global %>%
         filter(FDR < 0.05 & (logFC.treatmentDEX > 1 | logFC.treatmentCORT > 1))) %>% 
  write_cluster_list_to_xlsx(file_path = "data/pilotDexCortCtrl/r_outputs/statistics_results/edger_analysis_15.08.2025_v1/edgeR_AllTreatments_res0.2_FDR0.01log2FC1.xlsx")


res$cluster_0$CORT_vs_DEX %>% filter(gene_symbol == "FKBP5")

res %>% 
  keep(~ "CORT_vs_CTRL" %in% names(.x)) %>%                
  imap(~ .x$CORT_vs_CTRL %>%
         filter(FDR < 0.05) %>% 
         filter(logFC > 1 | logFC < -1)) %>% 
  write_cluster_list_to_xlsx(file_path = "data/pilotDexCortCtrl/r_outputs/statistics_results/edger_analysis_15.08.2025_v1/edgeR_CTRLvsCORT_res0.2_FDR0.01log2FC1.xlsx")

res %>% 
  keep(~ "DEX_vs_CTRL" %in% names(.x)) %>%                
  imap(~ .x$DEX_vs_CTRL %>%
         filter(FDR < 0.05) %>% 
         filter(logFC > 1 | logFC < -1)) %>% 
  write_cluster_list_to_xlsx(file_path = "data/pilotDexCortCtrl/r_outputs/statistics_results/edger_analysis_15.08.2025_v1/edgeR_CTRLvsDEX_res0.2_FDR0.01log2FC1.xlsx")

res %>% 
  keep(~ "CORT_vs_DEX" %in% names(.x)) %>%                
  imap(~ .x$CORT_vs_DEX %>%
         filter(FDR < 0.05) %>% 
         filter(logFC > 1 | logFC < -1)) %>% 
  write_cluster_list_to_xlsx(file_path = "data/pilotDexCortCtrl/r_outputs/statistics_results/edger_analysis_15.08.2025_v1/edgeR_CORTvsDEX_res0.2_FDR0.01log2FC1.xlsx")




dotplot_log2FC_correlation_v2 <- function(
    df,
    col_x = "logFC.treatmentDEX",
    col_y = "logFC.treatmentCORT",
    method = "pearson",
    fdr_col = "FDR",
    FDR_threshold = NULL,         
    thresholds = c(-1, 1),        
    threshold_x_text = NULL,      
    threshold_y_text = NULL,      
    point_alpha = 0.6,
    point_size  = 1.5,
    label_genes = TRUE,
    # visual parameters
    axis_text_size = 12,
    axis_title_size = 14,
    title_size = 16,
    legend_text_size = 12,
    label_size = 12,
    x_limits = NULL,
    y_limits = NULL,
    # colors
    color_sig_points = "red",
    # new
    title_prefix = ""   # <--- dodane
) {
  if (label_genes && !requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required for labeling.")
  }
  stopifnot(all(c(col_x, col_y) %in% colnames(df)))
  
  label_size_conv <- label_size / 3.88
  
  # klasyfikacja istotności (do koloru/etykiet)
  df$highlight <- "Not significant"
  sig_idx <- rep(FALSE, nrow(df))
  if (!is.null(FDR_threshold) && fdr_col %in% colnames(df)) {
    abs_thr <- if (length(thresholds)) max(thresholds) else 0
    sig_idx <- df[[fdr_col]] <= FDR_threshold &
      (abs(df[[col_x]]) >= abs_thr | abs(df[[col_y]]) >= abs_thr)
    df$highlight[sig_idx] <- "Significant"
  }
  n_sig <- sum(sig_idx, na.rm = TRUE)
  
  # globalna korelacja
  cor_all <- cor(df[[col_x]], df[[col_y]], method = method, use = "pairwise.complete.obs")
  
  # legenda
  legend_labels <- c(
    "Not significant" = "No significant",
    "Significant" = if (!is.null(FDR_threshold)) {
      sprintf("Significant (FDR < %.2f & |log2FC| ≥ %.1f)", FDR_threshold, if (length(thresholds)) max(thresholds) else 0)
    } else {
      "Significant"
    }
  )
  legend_labels <- gsub("0+$", "", legend_labels)
  legend_labels <- gsub("\\.$", "", legend_labels)
  
  # wykres
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[col_x]], y = .data[[col_y]])) +
    ggplot2::geom_point(ggplot2::aes(color = highlight), alpha = point_alpha, size = point_size) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, color = "black", fullrange = TRUE) +
    ggplot2::scale_color_manual(
      values = c("Not significant" = "grey70", "Significant" = color_sig_points),
      labels = legend_labels
    ) +
    ggplot2::labs(
      x = paste0("log2FC (", gsub("logFC\\.", "", col_x), ")"),
      y = paste0("log2FC (", gsub("logFC\\.", "", col_y), ")"),
      title = paste0(
        if (nzchar(title_prefix)) paste0(title_prefix, " | ") else "",
        sprintf("Correlation (%s): r = %.3f | n significant = %d", method, cor_all, n_sig)
      ),
      color = ""
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "bottom",
      axis.text  = ggplot2::element_text(size = axis_text_size,  color = "black"),
      axis.title = ggplot2::element_text(size = axis_title_size, color = "black"),
      plot.title = ggplot2::element_text(size = title_size, hjust = 0.5, color = "black"),
      legend.text = ggplot2::element_text(size = legend_text_size, color = "black")
    )
  
  # linie przerywane dla progów
  if (!is.null(thresholds) && length(thresholds) > 0) {
    for (thr in thresholds) {
      p <- p +
        ggplot2::geom_vline(xintercept = c(thr, -thr), linetype = "dashed", color = "grey50") +
        ggplot2::geom_hline(yintercept = c(thr, -thr), linetype = "dashed", color = "grey50")
    }
  }
  
  if (!is.null(x_limits)) p <- p + ggplot2::xlim(x_limits)
  if (!is.null(y_limits)) p <- p + ggplot2::ylim(y_limits)
  
  # etykiety genów
  if (label_genes && !is.null(FDR_threshold) && "gene_symbol" %in% colnames(df)) {
    label_idx <- sig_idx &
      (( !is.null(threshold_x_text) & abs(df[[col_x]]) >= threshold_x_text ) |
         ( !is.null(threshold_y_text) & abs(df[[col_y]]) >= threshold_y_text ))
    df_labels <- subset(df, label_idx & !is.na(gene_symbol))
    if (nrow(df_labels) > 0) {
      p <- p + ggrepel::geom_text_repel(
        data = df_labels,
        ggplot2::aes(label = gene_symbol),
        size = label_size_conv,
        max.overlaps = Inf
      )
      message(sprintf("Labeling %d genes that passed FDR and text thresholds.", nrow(df_labels)))
    }
  }
  
  return(list(cor_all = cor_all, n_sig = n_sig, plot = p))
}
