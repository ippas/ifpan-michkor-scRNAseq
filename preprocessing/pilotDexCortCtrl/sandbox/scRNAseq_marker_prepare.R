# ##############################################################################
# ---- prepare simple functions ----
# ##############################################################################
run_enrichr <- function(gene_list, database) {
  results <- enrichR::enrichr(gene_list, database)[[1]]
  
  results %>% 
    dplyr::select(!c("Old.P.value", "Old.Adjusted.P.value")) %>% 
    mutate(n_genes = as.integer(str_split_fixed(Overlap, "/", 2)[, 1]),
           set_size = as.integer(str_split_fixed(Overlap, "/", 2)[, 2])) %>%
    set_colnames(c("term", "overlap", "pvalue", "fdr", "odds_ratio", "combined_score", "genes", "n_genes", "set_size")) %>% 
    dplyr::select(c(term, overlap, n_genes, set_size, everything())) -> results
  
  return(results)
}


run_enrichr_for_all_clusters <- function(markers_list, database, gene_col = "gene_symbol") {
  res_list <- lapply(
    X = names(markers_list),
    FUN = function(cl) {
      res <- run_enrichr(
        gene_list = markers_list[[cl]][[gene_col]],
        database  = database
      )
      dplyr::mutate(res, cluster = cl)
    }
  )
  names(res_list) <- names(markers_list)
  return(res_list)
}



barplot_enrichr_results <- function(
    df,
    min_genes = 2,
    n_top = 10,
    title = "Enrichr results",
    axis_text_size = 12,
    axis_title_size = 14,
    legend_text_size = 10,
    legend_title_size = 12,
    plot_title_size = 16,
    label_text_size = 3.5,
    fill_limits = NULL,       # np. c(0, 10) lub NULL
    squish_to_limits = FALSE, # NEW: TRUE => wartości poza zakresem są „dociśnięte” do granic
    symmetric_limits = FALSE  # opcjonalnie: TRUE => automatyczne symetryczne limity ±max(|log2_or|)
) {
  df2 <- df %>%
    dplyr::filter(n_genes >= min_genes) %>%
    dplyr::arrange(fdr, dplyr::desc(odds_ratio)) %>%
    dplyr::slice_head(n = n_top) %>%
    dplyr::mutate(
      term = factor(term, levels = rev(term)),
      log2_or = log2(odds_ratio)
    )
  
  # Jeśli chcesz automatycznie symetryczne limity (przydatne dla gradient2 z midpoint=0)
  if (isTRUE(symmetric_limits) && is.null(fill_limits)) {
    m <- max(abs(df2$log2_or), na.rm = TRUE)
    fill_limits <- c(-m, m)
  }
  
  oob_fun <- if (isTRUE(squish_to_limits)) scales::squish else scales::censor
  
  ggplot(df2, aes(x = term, y = -log10(fdr), fill = log2_or)) +
    geom_col() +
    geom_text(aes(label = overlap), hjust = -0.2, size = label_text_size, color = "black") +
    coord_flip() +
    scale_fill_gradient2(
      low = scales::muted("blue"),
      mid = "white",
      high = scales::muted("red"),
      midpoint = 0,
      limits = fill_limits,
      oob = oob_fun              # <-- klucz: „squish” zamiast domyślnego „censor”
    ) +
    labs(
      title = title,
      x = NULL,
      y = expression(-log[10](FDR)),
      fill = expression(log[2]~"(Odds Ratio)")
    ) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = legend_title_size, hjust = 0.5, color = "black"),
      legend.text = element_text(size = legend_text_size, color = "black"),
      axis.text = element_text(size = axis_text_size, color = "black"),
      axis.title = element_text(size = axis_title_size, color = "black"),
      plot.title = element_text(size = plot_title_size, hjust = 0.5, color = "black")
    ) +
    guides(fill = guide_colorbar(title.position = "left", title.hjust = 0.5)) +
    expand_limits(y = max(-log10(df$fdr), na.rm = TRUE) * 1.1)
}

# ##############################################################################
# ---- define gene markers ----
# ##############################################################################

combined_BG60BLUE.markers_Minpct0.5wilcoxLogfc0.5 <- FindAllMarkers(
  object = combined_BG60BLUE,
  only.pos = FALSE,
  min.pct = 0.5, 
  test.use = "wilcox"
)


combined_BG60BLUE.markers_Minpct0.5wilcoxLogfc0.5 <- combined_BG60BLUE.markers_Minpct0.5wilcoxLogfc0.5 %>% 
  rownames_to_column(var = "rowname") %>% 
  select(-rowname) %>% 
  filter(p_val_adj < 0.05) %>% 
  mutate(name_id = gene) %>% 
  select(name_id, everything()) %>% 
  separate(gene, into = c("gene_id", "gene_symbol"), sep = "-", remove = TRUE, extra = "merge") %>%
  mutate(cluster = paste0("cluster_", cluster)) 
  


combined_BG60BLUE.markers_Minpct0.5wilcoxLogfc0.5 %>% 
  filter(avg_log2FC < -0.5) %>% 
  group_by(cluster) %>% 
  split(.$cluster) -> combined_BG60BLUE.markers_downMinpct0.5wilcoxLogfc0.5


combined_BG60BLUE.markers_Minpct0.5wilcoxLogfc0.5 %>% 
  filter(avg_log2FC > 0.5) %>% 
  group_by(cluster) %>% 
  split(.$cluster) -> combined_BG60BLUE.markers_upMinpct0.5wilcoxLogfc0.5

# ##############################################################################
# ---- summary gene markers ----
# ##############################################################################
combined_BG60BLUE.markers_upMinpct0.5wilcoxLogfc0.5 %>% 
  lapply(., nrow) %>% 
  unlist %>% summary


combined_BG60BLUE.markers_downMinpct0.5wilcoxLogfc0.5 %>% 
  lapply(., nrow) %>% 
  unlist %>% summary


# ##############################################################################
# ---- visualization gene markers ----
# ##############################################################################
combined_BG60BLUE.markers_upMinpct0.5wilcoxLogfc0.5 %>%
  .$cluster_0 %>%
  head(9) %>%
  .$name_id -> vector_name_id

FeaturePlot(combined_BG60BLUE, pt.size = 0.5,
            features = vector_name_id)


top9_markersUp_plots_list <- combined_BG60BLUE.markers_upMinpct0.5wilcoxLogfc0.5 %>%
  map(~ .x %>%
        slice_head(n = 9) %>%
        pull(name_id) %>%
        ((genes) FeaturePlot(
          combined_BG60BLUE,
          features = genes,
          pt.size = 0.5
        ))()
  ) %>%
  set_names(names(.))


top9_markersDown_plots_list <- combined_BG60BLUE.markers_downMinpct0.5wilcoxLogfc0.5 %>%
  map(~ .x %>%
        slice_head(n = 9) %>%
        pull(name_id) %>%
        (\(genes) FeaturePlot(
          combined_BG60BLUE,
          features = genes,
          pt.size = 0.5
        ))()
  ) %>%
  set_names(names(.))


top9_markersUp_plots_list$cluster_0 


# ##############################################################################
# ---- first save to file ----
# ##############################################################################
combined_BG60BLUE.markers_upMinpct0.5wilcoxLogfc0.5 %>% 
  write_cluster_list_to_xlsx(., file_path = "data/pilotDexCortCtrl/r_outputs/cluster_markers/clusterMarkers_combinedBG60BLUE_vstFeatures2000dims20res0.5_upFDR0.05log2FC0.5.xlsx")
  
combined_BG60BLUE.markers_downMinpct0.5wilcoxLogfc0.5 %>% 
  write_cluster_list_to_xlsx(., file_path = "data/pilotDexCortCtrl/r_outputs/cluster_markers/clusterMarkers_combinedBG60BLUE_vstFeatures2000dims20res0.5_downFDR0.05log2FC0.5.xlsx")


# ##############################################################################
# ---- second save to file ----
# ##############################################################################
write_cluster_listAndPlots_to_xlsx(
  cluster_list = combined_BG60BLUE.markers_upMinpct0.5wilcoxLogfc0.5,
  plot_list    = top9_markersUp_plots_list,
  file_path    = "data/pilotDexCortCtrl/r_outputs/cluster_markers/clusterMarkers_combinedBG60BLUE_up_withPlots.xlsx",
  img_width    = 12, img_height = 10, img_dpi = 300
)

# ##############################################################################
# ---- enrichr analysysis ----
# ##############################################################################

# Lista baz do analizy
db_list <- c(
  "CellMarker_2024",
  "CellMarker_Augmented_2021",
  "Human_Gene_Atlas",
  "ARCHS4_Cell-lines"
)

# Uruchomienie analizy
enrichrResults_perDatabase_perCluster <- lapply(
  db_list,
  function(db) run_enrichr_for_all_clusters(
    combined_BG60BLUE.markers_upMinpct0.5wilcoxLogfc0.5,
    database = db,
    gene_col = "gene_symbol"
  )
)

# Nadaj nazwy baz jako klucze listy
names(enrichrResults_perDatabase_perCluster) <- db_list



# ##############################################################################
# ---- prepare plots from enrichr ----combined_BG60BLUE.markers_upMinpct0.5wilcoxLogfc0.5
# ##############################################################################

clusters <- names(enrichrResults_perDatabase_perCluster[[ db_list[1] ]])

plots_perCluster_perDatabase <- lapply(clusters, function(cl) {
  wrap_plots(
    plotlist = lapply(db_list, function(db) {
      barplot_enrichr_results(
        df = enrichrResults_perDatabase_perCluster[[db]][[cl]],
        min_genes = 2,
        n_top = 10,
        title = paste("Cluster", cl, "–", db),
        axis_text_size = 12,
        axis_title_size = 12,
        legend_text_size = 10,
        legend_title_size = 12,
        plot_title_size = 14,
        label_text_size = 4,
        fill_limits = c(0, 10),
        squish_to_limits = TRUE
      )
    }),
    ncol = 1
  ) +
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
})


names(plots_perCluster_perDatabase) <- clusters


# ##############################################################################
# ---- third save to xlsx ----
# ##############################################################################

write_cluster_listAndTwoPlots_to_xlsx(
  cluster_list = combined_BG60BLUE.markers_upMinpct0.5wilcoxLogfc0.5,
  plot_list_a  = top9_markersUp_plots_list,
  plot_list_b  = plots_perCluster_perDatabase,
  file_path    = "data/pilotDexCortCtrl/r_outputs/cluster_markers/clusterMarkersAndPlotsExamplesEnrichr_combinedBG60BLUE_vstFeatures2000dims20res0.5_upFDR0.05log2FC0.5.xlsx",
  img_width_a  = 11, img_height_a = 9,  # rozmiar 1. wykresu
  img_width_b  = 12, img_height_b = 13,  # rozmiar 2. wykresu
  gap_rows     = 45
)


# ##############################################################################
# ---- test enrichr ----
# ##############################################################################

# enrichrResults_perDatabase_perCluster
# 
# # Pojedynczy wykres – np. top 10 termów z min. 2 genami
# barplot_enrichr_results(
#   df = enrichr_cellmarker2024$cluster_0,
#   min_genes = 2,
#   n_top = 10,
#   title = "Cluster 11 – CellMarker 2024",
#   axis_text_size = 16,
#   axis_title_size = 16,
#   legend_text_size = 12,
#   legend_title_size = 16,
#   plot_title_size = 16,
#   label_text_size = 5,
#   fill_limits = c(0, 10)
# ) 
# 



