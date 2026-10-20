edgeR_anovaOneFactor_treatment_pseudobulk_res0.5_clusters13_results$cluster_Astrocytes$f_test_global %>% 
  filter(FDR < 0.05) %>% head


simple_summary_edger <- function(
    data,
    fdr_threshold = 0.05,
    logfc_threshold = 1
) {
  # 🔹 Filtracja dla DEX_vs_CTRL
  dex_genes <- data$DEX_vs_CTRL %>% 
    filter(FDR < fdr_threshold) %>% 
    filter(abs(logFC) > logfc_threshold) %>% 
    pull(gene_symbol)
  
  # 🔹 Filtracja dla CORT_vs_CTRL
  cort_genes <- data$CORT_vs_CTRL %>% 
    filter(FDR < fdr_threshold) %>% 
    filter(abs(logFC) > logfc_threshold) %>% 
    pull(gene_symbol)
  
  # 🔹 Obliczenia
  n_unique_genes <- length(unique(c(cort_genes, dex_genes)))
  common_part    <- intersect(dex_genes, cort_genes)
  
  # 🔹 Wynik
  tibble(
    fdr_threshold     = fdr_threshold,
    logfc_threshold   = logfc_threshold,
    dexVsCTRL_sig     = length(dex_genes),
    cortVsCTRL_sig    = length(cort_genes),
    n_common_genes    = length(common_part),
    n_unique_genes    = n_unique_genes,
    percent_common_genes = round(length(common_part) / n_unique_genes * 100, 2)
  )
}

simple_summary_edger(data = edgeR_anovaOneFactor_treatment_pseudobulk_res0.5_clusters13_results$cluster_Astrocytes)


map_dfr(
  edgeR_anovaOneFactor_treatment_pseudobulk_res0.4_clusters16_celltypes11_results ,
  simple_summary_edger,
  fdr_threshold = 0.01,
  .id = "celltype"
) %>% 
  mutate(celltype = str_replace_all(celltype, "cluster_", "")) %>% 
  arrange(desc(n_unique_genes)) %>% 
  select(-c(fdr_threshold, logfc_threshold))



edgeR_anovaOneFactor_treatment_pseudobulk_res0.5_clusters13_results %>% 
  lapply(., function(x){
    x$CORT_vs_CTRL %>% 
      filter(FDR < 0.01) %>% 
      filter(abs(logFC) > 1) %>% 
      .$gene_symbol
  }) %>% 
  unname %>% 
  unlist %>% 
  table %>% 
  as.data.frame() %>% 
  set_colnames(c("gene_symbol", "freq")) %>% 
  arrange(desc(freq)) %>% 
  nrow


# ##############################################################################
map_dfr(
  edgeR_anovaOneFactor_treatment_pseudobulk_res0.3_clusters16cells15_results,
  simple_summary_edger,
  fdr_threshold = 0.01,
  .id = "celltype"
) %>% 
  mutate(celltype = str_replace_all(celltype, "cluster_", "")) %>% 
  arrange(desc(n_unique_genes)) %>% 
  select(-c(fdr_threshold, logfc_threshold))



edgeR_anovaOneFactor_treatment_pseudobulk_res0.5_clusters13_results %>% 
  lapply(., function(x){
    x$CORT_vs_CTRL %>% 
      filter(FDR < 0.01) %>% 
      filter(abs(logFC) > 1) %>% 
      .$gene_symbol
  }) %>% 
  unname %>% 
  unlist %>% 
  table %>% 
  as.data.frame() %>% 
  set_colnames(c("gene_symbol", "freq")) %>% 
  arrange(desc(freq)) %>% 
  nrow


genes <- c(
  "IRAK2","ZNFX1","NAPEPLD","GOT1L1","LIMK2","FAM110B","FGFR2","SYBU","SNTA1","RPRML",
  "CPEB2","SOX17","KLRG2","MTHFD2","PPP1R3B","SDC4","FOS","GRAMD1B","APBB3","SPHK1",
  "NFIL3","GADD45G","AFF1","STIM1","ZBTB16","ARL4A","NET1","SLC25A33","CDKN1A","RGCC",
  "LCN2","TES","MAFF","ARL4D","HSPA2","ADRB2","TLE1","MRAP","LRRC8A","SLC16A9",
  "GDPD1","MS4A8","DCUN1D3","LNX2","SYDE2","KIAA0040","ENC1","IL15RA","PDK4","RANBP9",
  "AMPD3","ADAMTS15","RNF169","LARP4","ARRDC3","DUSP18","TIPARP","RGS2","LYVE1","RHOB",
  "TP53INP1","ACSS1","USP54","SLC10A6","MAP3K6","SLC7A8","SESN1","FAM107A","MS4A6E",
  "DUSP1","DDIT4","TSC22D3","HIF3A","ERRFI1","FKBP5"
)



# ##############################################################################

simple_summary_edger_v2 <- function(
    data,
    fdr_threshold = 0.05,
    logfc_threshold = 1,
    reference_genes = NULL
) {
  # 🔹 Filtracja dla DEX_vs_CTRL
  dex_genes <- data$DEX_vs_CTRL %>% 
    filter(FDR < fdr_threshold) %>% 
    filter(abs(logFC) > logfc_threshold) %>% 
    pull(gene_symbol)
  
  # 🔹 Filtracja dla CORT_vs_CTRL
  cort_genes <- data$CORT_vs_CTRL %>% 
    filter(FDR < fdr_threshold) %>% 
    filter(abs(logFC) > logfc_threshold) %>% 
    pull(gene_symbol)
  
  # 🔹 Obliczenia podstawowe
  n_unique_genes <- length(unique(c(cort_genes, dex_genes)))
  common_part    <- intersect(dex_genes, cort_genes)
  
  # 🔹 Nakładanie z wektorem referencyjnym
  if (!is.null(reference_genes)) {
    overlap_dex   <- intersect(dex_genes, reference_genes)
    overlap_cort  <- intersect(cort_genes, reference_genes)
    overlap_common <- intersect(common_part, reference_genes)
  } else {
    overlap_dex <- overlap_cort <- overlap_common <- character(0)
  }
  
  # 🔹 Wynik
  tibble(
    fdr_threshold     = fdr_threshold,
    logfc_threshold   = logfc_threshold,
    dexVsCTRL_sig     = length(dex_genes),
    cortVsCTRL_sig    = length(cort_genes),
    n_common_genes    = length(common_part),
    n_unique_genes    = n_unique_genes,
    percent_common_genes = round(length(common_part) / n_unique_genes * 100, 2),
    
    # 🔸 Nowe kolumny z wynikami nakładania
    dex_overlap_with_ref   = length(overlap_dex),
    cort_overlap_with_ref  = length(overlap_cort),
    common_overlap_with_ref = length(overlap_common),
    percent_ref_in_common  = ifelse(
      !is.null(reference_genes) && length(reference_genes) > 0,
      round(length(overlap_common) / length(reference_genes) * 100, 2),
      NA_real_
    )
  )
}


map_dfr(
  edgeR_anovaOneFactor_treatment_pseudobulk_res0.3_clusters16cells15_results,
  simple_summary_edger_v2,
  fdr_threshold = 0.01,
  .id = "celltype",
  reference_genes = genes
) %>% 
  mutate(celltype = str_replace_all(celltype, "cluster_", "")) %>% 
  arrange(desc(n_unique_genes)) %>% 
  select(-c(fdr_threshold, logfc_threshold))

# ##############################################################################

simple_summary_edger <- function(
    data,
    fdr_threshold = 0.05,
    logfc_threshold = 1,
    genes_of_interest = c("FKBP5", "TSC22D3", "SGK1", "ZBTB16", "PDK4")
) {
  # 🔹 Filtracja dla DEX_vs_CTRL
  dex_genes <- data$DEX_vs_CTRL %>% 
    filter(FDR < fdr_threshold, abs(logFC) > logfc_threshold) %>% 
    pull(gene_symbol)
  
  # 🔹 Filtracja dla CORT_vs_CTRL
  cort_genes <- data$CORT_vs_CTRL %>% 
    filter(FDR < fdr_threshold, abs(logFC) > logfc_threshold) %>% 
    pull(gene_symbol)
  
  # 🔹 Zbiór wszystkich genów spełniających kryteria
  all_genes <- unique(c(dex_genes, cort_genes))
  
  # 🔹 Obliczenia ogólne
  n_unique_genes <- length(all_genes)
  presence_tbl <- setNames(
    as.list(genes_of_interest %in% all_genes),
    paste0("has_", genes_of_interest)
  )
  
  # 🔹 Wynik
  tibble(
    fdr_threshold     = fdr_threshold,
    logfc_threshold   = logfc_threshold,
    dexVsCTRL_sig     = length(dex_genes),
    cortVsCTRL_sig    = length(cort_genes),
    n_common_genes    = length(common_part),
    n_unique_genes    = n_unique_genes,
    percent_common_genes = round(length(common_part) / n_unique_genes * 100, 2)
  ) %>%
    bind_cols(as_tibble(presence_tbl))
}


map_dfr(
  edgeR_anovaOneFactor_treatment_pseudobulk_res0.4_clusters16_celltypes11_results,
  simple_summary_edger,
  fdr_threshold = 0.01,
  logfc_threshold = 1,
  genes_of_interest = c("FKBP5", "TSC22D3", "SGK1", "ZBTB16",  "PER1", "DDIT4", "HIF3A", "OLIG2", "DUSP1"),
  .id = "celltype"
) %>%
  mutate(celltype = str_replace_all(celltype, "cluster_", "")) %>%
  arrange(desc(n_unique_genes)) %>%
  select(-c(fdr_threshold, logfc_threshold)) %>% as.data.frame()
