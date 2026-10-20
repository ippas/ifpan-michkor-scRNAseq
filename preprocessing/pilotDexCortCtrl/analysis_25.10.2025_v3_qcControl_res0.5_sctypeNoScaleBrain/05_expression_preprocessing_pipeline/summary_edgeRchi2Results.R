summarize_upregulated_grGenes_logFC_CPM_p_long <- function(
    df,
    contrast     = c("DEX", "CORT"),
    gr_col       = "is_grSignature",
    cpm_col      = "max_CPM",
    gene_col     = "gene_symbol"
) {
  
  contrast <- match.arg(contrast)
  
  if (contrast == "DEX") {
    p_col     <- "PValue_edgerDEX"
    logfc_col <- "log2FC_edgerDEX"
  } else {
    p_col     <- "PValue_edgerCORT"
    logfc_col <- "log2FC_edgerCORT"
  }
  
  # remove NA
  df <- df %>%
    dplyr::filter(
      !is.na(.data[[p_col]]),
      !is.na(.data[[logfc_col]]),
      !is.na(.data[[cpm_col]])
    )
  
  p_thresholds     <- c(1, 0.05, 0.01, 0.001, 0.0001, 0.00001)
  logfc_thresholds <- c(0, 0.5, 1, 1.5, 2)
  cpm_thresholds   <- c(0, 10, 25, 50, 100)
  
  out_list <- list()
  
  for (cpm_thr in cpm_thresholds) {
    
    df_cpm <- df[df[[cpm_col]] >= cpm_thr, , drop = FALSE]
    
    all_genes_total <- nrow(df_cpm)
    gr_total        <- sum(df_cpm[[gr_col]], na.rm = TRUE)
    
    for (p_thr in p_thresholds) {
      
      df_p <- df_cpm[df_cpm[[p_col]] < p_thr, , drop = FALSE]
      
      for (fc in logfc_thresholds) {
        
        # ======================
        # UP-REGULATION (logFC >= threshold)
        # ======================
        df_up <- df_p[df_p[[logfc_col]] >= fc, , drop = FALSE]
        
        n_all_up <- nrow(df_up)
        n_gr_up  <- sum(df_up[[gr_col]], na.rm = TRUE)
        
        out_list <- append(out_list, list(
          tibble::tibble(
            contrast          = contrast,
            cpm_threshold     = cpm_thr,
            p_value_threshold = p_thr,
            logFC_threshold   = fc,
            direction         = "UP",
            class             = "all",
            count             = n_all_up,
            all_genes         = all_genes_total,
            all_GR_genes      = gr_total,
            genes             = list(df_up[[gene_col]])
          )
        ))
        
        out_list <- append(out_list, list(
          tibble::tibble(
            contrast          = contrast,
            cpm_threshold     = cpm_thr,
            p_value_threshold = p_thr,
            logFC_threshold   = fc,
            direction         = "UP",
            class             = "GR",
            count             = n_gr_up,
            all_genes         = all_genes_total,
            all_GR_genes      = gr_total,
            genes             = list(df_up[[gene_col]][df_up[[gr_col]] == 1])
          )
        ))
        
        # ======================
        # DOWN-REGULATION (logFC <= -threshold)
        # ======================
        df_down <- df_p[df_p[[logfc_col]] <= -fc, , drop = FALSE]
        
        n_all_down <- nrow(df_down)
        n_gr_down  <- sum(df_down[[gr_col]], na.rm = TRUE)
        
        out_list <- append(out_list, list(
          tibble::tibble(
            contrast          = contrast,
            cpm_threshold     = cpm_thr,
            p_value_threshold = p_thr,
            logFC_threshold   = fc,
            direction         = "DOWN",
            class             = "all",
            count             = n_all_down,
            all_genes         = all_genes_total,
            all_GR_genes      = gr_total,
            genes             = list(df_down[[gene_col]])
          )
        ))
        
        out_list <- append(out_list, list(
          tibble::tibble(
            contrast          = contrast,
            cpm_threshold     = cpm_thr,
            p_value_threshold = p_thr,
            logFC_threshold   = fc,
            direction         = "DOWN",
            class             = "GR",
            count             = n_gr_down,
            all_genes         = all_genes_total,
            all_GR_genes      = gr_total,
            genes             = list(df_down[[gene_col]][df_down[[gr_col]] == 1])
          )
        ))
        
      }
    }
  }
  
  dplyr::bind_rows(out_list)
}

purrr::map_dfr(
  names(combined_summaryEdgerChi2_res0.3RNAnoScaledBrain),
  function(cl) {
    
    df_cluster <- combined_summaryEdgerChi2_res0.3RNAnoScaledBrain[[cl]] %>%
      dplyr::mutate(is_grSignature = gene_symbol %in% c(genes_global, genes_neural, genes_GlobalDown, genes_neuralDown))
    
    summarize_upregulated_grGenes_logFC_CPM_p_long(
      df = df_cluster,
      contrast = "DEX",
      cpm_col  = "meanCPM_DEX"
    ) %>% dplyr::mutate(cluster = cl)
  }
) -> summary_all_clusters_DEX

purrr::map_dfr(
  names(combined_summaryEdgerChi2_res0.3RNAnoScaledBrain),
  function(cl) {
    
    df_cluster <- combined_summaryEdgerChi2_res0.3RNAnoScaledBrain[[cl]] %>%
      dplyr::mutate(is_grSignature = gene_symbol %in% c(genes_global, genes_neural, genes_GlobalDown, genes_neuralDown))
    
    summarize_upregulated_grGenes_logFC_CPM_p_long(
      df = df_cluster,
      contrast = "CORT",
      cpm_col  = "meanCPM_CORT"
    ) %>% dplyr::mutate(cluster = cl)
  }
) -> summary_all_clusters_CORT



  
rbind(summary_all_clusters_DEX, summary_all_clusters_CORT) %>% 
  filter(
    cluster %in% c("Astrocytes", "Radial glial cells", "Oligodendrocyte precursor cells",
                   "Oligodendrocytes", "Immature neurons", "Mature neurons"),
    class == "GR",
    cpm_threshold == 50,
    p_value_threshold %in% c(0.05, 0.01, 0.001, 0.0001, 0.00001),
    contrast %in% c("CORT", "DEX")
  ) %>%
  ggplot(
    aes(
      x = factor(logFC_threshold),
      y = count,
      group = cluster,
      color = cluster
    )
  ) +
  
  geom_line(linewidth = 1) +
  geom_point(size = 4) +
  
  geom_text_repel(
    aes(label = count),
    size = 5,
    nudge_y = 40,
    min.segment.length = 0,
    box.padding = 0.5,
    point.padding = 0.5,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  
  facet_grid(p_value_threshold ~ contrast) +
  
  theme_bw() +
  labs(
    title = "GR-dependent activated genes across logFC thresholds\nCPM ≥ 100",
    x = "log2FC threshold",
    y = "Number of GR genes"
  ) +
  scale_color_brewer(palette = "Dark2") +
  
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    
    strip.text = element_text(size = 18, face = "bold"),  # facet labels: p-value & contrast
    plot.title = element_text(size = 22, face = "bold")
  )
# ##############################################################################
rbind(summary_all_clusters_DEX, summary_all_clusters_CORT) %>% 
  filter(cluster %in% c(
    "Astrocytes",
    "Immature neurons", 
    "Mature neurons", 
    "Oligodendrocyte precursor cells", 
    "Oligodendrocytes",
    "Radial glial cells"
  )) %>% 
  filter(cpm_threshold == 50) %>% 
  filter(p_value_threshold == 0.00001) %>% 
  filter(class == "all") %>%
  # filter(contrast == "DEX") %>% 
  filter(logFC_threshold == 1) %>% .$genes %>% unlist %>% unique %>% length()
  filter(cluster == "Mature neurons") %>% 




# ##############################################################################
# ---- summary log2FC results ----
# ##############################################################################

create_logfc_summary_table <- function(
    summary_all_clusters_DEX,
    summary_all_clusters_CORT,
    contrast_filter = "CORT",
    cpm_thr = 50,
    p_thr = 0.01,
    celltypes = c(
      "Radial glial cells",
      "Astrocytes",
      "Oligodendrocyte precursor cells",
      "Oligodendrocytes",
      "Immature neurons",
      "Mature neurons"
    )
) {
  library(dplyr)
  library(tidyr)
  library(purrr)
  
  # ============================
  # 0) USTAWIENIA
  # ============================

  
  logfc_levels <- c(0, 0.5, 1, 1.5, 2)
  logfc_labels <- paste0("log2FC > ", logfc_levels)
  
  # ============================
  # 1) SCAL DANE I SPRAWDŹ KOLUMNY
  # ============================
  df_raw <- dplyr::bind_rows(summary_all_clusters_DEX, summary_all_clusters_CORT)
  
  # kolumna z progami logFC
  lfc_col <- dplyr::case_when(
    "logFC_threshold"  %in% names(df_raw) ~ "logFC_threshold",
    "logfc_threshold"  %in% names(df_raw) ~ "logfc_threshold",
    TRUE ~ NA_character_
  )
  if (is.na(lfc_col)) {
    stop("Nie znalazłem kolumny z progami logFC (logFC_threshold / logfc_threshold) w summary_all_clusters_*")
  }
  
  if (!"cluster" %in% names(df_raw)) {
    stop("Brakuje kolumny 'cluster' w summary_all_clusters_* (typ komórki).")
  }
  if (!"class" %in% names(df_raw)) {
    stop("Brakuje kolumny 'class' (all/GR) w summary_all_clusters_*.")
  }
  if (!"count" %in% names(df_raw)) {
    stop("Brakuje kolumny 'count' w summary_all_clusters_*.")
  }
  
  has_genes <- "genes" %in% names(df_raw)
  
  # ============================
  # 2) FILTRUJEMY DANE
  # ============================
  df <- df_raw %>%
    dplyr::filter(
      cluster %in% celltypes,
      contrast == contrast_filter,
      cpm_threshold == cpm_thr,
      p_value_threshold == p_thr
    ) %>%
    dplyr::mutate(
      logFC_thr_num = as.numeric(.data[[lfc_col]])
    ) %>%
    dplyr::filter(logFC_thr_num %in% logfc_levels)
  
  # jeśli po filtrach nic nie zostało:
  if (nrow(df) == 0) {
    stop("Po filtrach (contrast/cpm/pvalue) df jest puste – sprawdź progi.")
  }
  
  # ============================
  # 3) TABELA GŁÓWNA: klaster × logFC
  # ============================
  main_tab <- df %>%
    dplyr::group_by(cluster, logFC_thr_num, class) %>%
    dplyr::summarise(count = sum(count), .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from  = class,
      values_from = count,
      values_fill = 0
    ) %>%
    dplyr::mutate(value = sprintf("%d (%d)", all, GR)) %>%
    dplyr::select(cluster, logFC_thr_num, value) %>%
    tidyr::pivot_wider(
      names_from  = logFC_thr_num,
      values_from = value
    )
  
  # upewnij się, że kolumny są w kolejności logfc_levels
  # (mogą mieć nazwy "0", "0.5", "1", ...)
  logfc_colnames <- intersect(as.character(logfc_levels), names(main_tab))
  main_tab <- main_tab %>%
    dplyr::select(cluster, dplyr::all_of(logfc_colnames))
  
  # nadaj ładne etykiety kolumnom log2FC >
  names(main_tab) <- c(
    "cluster",
    logfc_labels[match(logfc_colnames, as.character(logfc_levels))]
  )
  
  # uporządkuj kolejność wierszy (typy komórek)
  main_tab <- main_tab %>%
    dplyr::arrange(factor(cluster, levels = celltypes))
  
  # ============================
  # 4) WIERSZ 7: UNIKALNE GENY PER LOGFC
  # ============================
  if (has_genes) {
    # df ma kolumnę genes (list-column z symbolami)
    summary_genes <- df %>%
      dplyr::group_by(logFC_thr_num, class) %>%
      dplyr::summarise(
        genes_union = list(unique(unlist(genes))),
        .groups = "drop"
      ) %>%
      dplyr::mutate(n_genes = lengths(genes_union)) %>%
      dplyr::select(logFC_thr_num, class, n_genes) %>%
      tidyr::pivot_wider(
        names_from  = class,
        values_from = n_genes,
        values_fill = 0
      ) %>%
      dplyr::mutate(value = sprintf("%d (%d)", all, GR))
    
    # wektor w tej samej kolejności progów co w tabeli
    gene_row_values <- sapply(
      logfc_levels,
      function(thr) {
        v <- summary_genes %>%
          dplyr::filter(logFC_thr_num == thr) %>%
          dplyr::pull(value)
        if (length(v) == 0) "" else v
      }
    )
    
    # budujemy wiersz
    summary_row <- tibble::tibble(
      cluster = "Total unique genes",
      !!!setNames(as.list(gene_row_values), logfc_labels)
    )
    
    # dopasuj kolumny do main_tab (na wypadek braków)
    missing_cols <- setdiff(names(main_tab), names(summary_row))
    if (length(missing_cols) > 0) {
      summary_row[missing_cols] <- ""
    }
    summary_row <- summary_row %>%
      dplyr::select(dplyr::all_of(names(main_tab)))
    
    final_tab <- dplyr::bind_rows(main_tab, summary_row)
    
  } else {
    # brak kolumny genes – robimy tylko tabelę główną
    final_tab <- main_tab
  }
  
  final_tab %>% 
    dplyr::rename(celltype = "cluster") -> final_tab
  
  return(final_tab)
}


summary_table <- create_logfc_summary_table(
  summary_all_clusters_DEX ,
  summary_all_clusters_CORT,
  contrast_filter = "CORT",
  cpm_thr = 50,
  p_thr = 0.01
)


summary_table <- create_logfc_summary_table(
  summary_all_clusters_DEX,
  summary_all_clusters_CORT,
  contrast_filter = "BOTH",
  cpm_thr = 50,
  p_thr = 0.01
)


summary_table

