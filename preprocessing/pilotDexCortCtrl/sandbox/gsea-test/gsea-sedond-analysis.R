if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::source_url("https://raw.githubusercontent.com/mateuszzieba/ifpan-matzieb-geneTranslator/master/preprocessing/R/src/simple_biomart_tools.r")

# --- wybrane typy komórek glejowych/neuronów ---
selected_brain_celltypes <- c(
  "Astrocytes", "Bergmann_glia", "Microglia",
  "Oligodendrocyte_progenitor_cells", "Radial_glia_cells",
  "Satellite_glial_cells", "Tanycytes",
  "Neurons", "Immature_neurons", "Neural_stem/precursor_cells"
)

# --- pobranie genów białkowych z Ensembl v114 ---
genes_all <- download_biomart_genes(version = "112")
genes_proteinCoding_v114 <- genes_all %>%
  dplyr::filter(
    gene_biotype == "protein_coding",
    chromosome_name %in% c(1:22, "X", "Y"),
    external_gene_name != ""
  )

# --- geny obecne w Seurat object ---
genes_df_seuratObj <- data.frame(
  name_id = rownames(integrated_pilotDexCortCtrl[["RNA"]]),
  stringsAsFactors = FALSE
) %>%
  tidyr::separate(name_id, into = c("ensembl_id", "gene_symbol"),
                  sep = "-", extra = "drop") %>%
  dplyr::select(ensembl_id)

genes_df_proteinCoding_seuratObj <- genes_df_seuratObj %>%
  dplyr::left_join(genes_proteinCoding_v114,
                   by = c("ensembl_id" = "ensembl_gene_id")) %>%
  tidyr::drop_na() %>%
  dplyr::rename(gene_symbol = external_gene_name)

# --- agregacja średnich ekspresji per klaster ---
agg_list_mean <- aggregate_counts_by_cluster_sample(
  seurat_obj  = integrated_pilotDexCortCtrl,
  cluster_col = "seurat_clusters_res0.13",
  sample_col  = "sample_id",
  assay       = "RNA",
  layer       = "counts",
  agg         = "mean",
  return      = "list",
  progress    = TRUE
)

agg_list_mean <- filter_list_by_ensemblID(
  gene_expression_list = agg_list_mean,
  ensembl_ids = genes_df_proteinCoding_seuratObj$ensembl_id,
  return_only_geneSymbol = TRUE
)

# --- średnia ekspresja genu w każdym klastrze ---
agg_mean_expr <- sapply(agg_list_mean, Matrix::rowMeans)

# --- funkcja: max-normalizacja × tau ---
combine_specificity_mult <- function(expr_mat) {
  spec_max <- t(apply(expr_mat, 1, function(x) {
    if (max(x) == 0) return(rep(0, length(x)))
    x / max(x)
  }))
  spec_tau <- apply(expr_mat, 1, function(x) {
    if (max(x) == 0) return(NA)
    1 - (sum(x / max(x)) / length(x))
  })
  spec_combined <- spec_max * spec_tau
  list(max = spec_max, tau = spec_tau, combined = spec_combined)
}

# --- funkcja: max-normalizacja × tau ---
combine_specificity_mult <- function(expr_mat) {
  # Max-normalization
  spec_max <- t(apply(expr_mat, 1, function(x) {
    if (max(x) == 0) return(rep(0, length(x)))
    x / max(x)
  }))
  
  # Tau specificity
  spec_tau <- apply(expr_mat, 1, function(x) {
    if (max(x) == 0) return(NA)
    1 - (sum(x / max(x)) / length(x))
  })
  
  # Combined = spec_max × tau
  spec_combined <- spec_max * spec_tau
  
  # --- automatyczne nadanie nazw kolumn ---
  if (!is.null(colnames(expr_mat))) {
    colnames(spec_max)     <- colnames(expr_mat)
    colnames(spec_combined) <- colnames(expr_mat)
  }
  rownames(spec_max)      <- rownames(expr_mat)
  rownames(spec_combined) <- rownames(expr_mat)
  
  list(max = spec_max, tau = spec_tau, combined = spec_combined)
}

# --- wyliczenie specyficzności ---
tmp <- combine_specificity_mult(agg_mean_expr) %>% .$combined %>%  na.omit()

tmp %>% head

colnames(tmp) <- c("cluster_0", "cluster_1", "cluster_2", "cluster_3",
                   "cluster_4", "cluster_5", "cluster_6",
                   "cluster_7", "cluster_8", "cluster_9")

# --- fgsea na wybranych typach komórek ---
res <- fgsea(
  pathways = panglaoDB_list_BrainMarkers[selected_brain_celltypes],
  stats    = tmp %>% na.omit() %>% .[, 1],  # kolumna-klaster do analizy
  nperm    = 1000
)

res %>% na.omit() %>% arrange(desc(NES))


# ---- run analysis for all clusters  ----
library(dplyr)
library(fgsea)
library(purrr)
library(tibble)

run_fgsea_all <- function(mat, pathways, nperm = 1000) {
  # 1) usuń wiersze, gdzie wszystkie kolumny to NA/NaN/Inf
  mat <- mat %>% dplyr::filter(dplyr::if_any(dplyr::everything(), ~ is.finite(.)))
  
  # 2) funkcja pomocnicza do fgsea dla pojedynczego klastra
  run_fgsea_one <- function(vec, gene_names, pathways, nperm = 1000) {
    v  <- as.numeric(vec)
    ok <- is.finite(v)
    stats <- v[ok]
    names(stats) <- gene_names[ok]
    stats <- sort(stats, decreasing = TRUE)
    fgsea(pathways = pathways, stats = stats, nperm = nperm)
  }
  
  # 3) fgsea dla wszystkich klastrów
  res_list <- setNames(
    lapply(colnames(mat), function(cl) {
      out <- run_fgsea_one(
        vec        = mat[[cl]],
        gene_names = rownames(mat),
        pathways   = pathways,
        nperm      = nperm
      )
      
      if (nrow(out) == 0) {
        return(tibble::tibble(
          pathway = character(), pval = numeric(), padj = numeric(),
          ES = numeric(), NES = numeric(), nMoreExtreme = integer(),
          size = integer(), cluster = character(), genes = character()
        ))
      }
      
      out$cluster <- cl
      out <- out %>%
        dplyr::mutate(genes = purrr::map_chr(leadingEdge, ~ paste(.x, collapse = "|"))) %>%
        dplyr::select(-leadingEdge) %>%
        dplyr::arrange(dplyr::desc(NES))
      out
    }),
    colnames(mat)
  )
  
  # 4) wybierz pierwszy wynik (najwyższy NES) na klaster
  first_per_cluster <- res_list %>%
    purrr::imap(~ if (nrow(.x)) dplyr::slice_max(.x, NES, n = 1, with_ties = FALSE) else NULL) %>%
    purrr::discard(is.null) %>%
    dplyr::bind_rows() %>%
    dplyr::arrange(dplyr::desc(NES))
  
  # 5) zwróć listę
  list(
    results = res_list,
    first_per_cluster = first_per_cluster
  )
}


mat <- combine_specificity_mult(agg_mean_expr)$combined %>% as.data.frame()

colnames(mat) <- c("cluster_0", "cluster_1", "cluster_2", "cluster_3",
                   "cluster_4", "cluster_5", "cluster_6",
                   "cluster_7", "cluster_8", "cluster_9")

panglaoDB_all_markers %>% 
  # filter(organ == "Brain") %>% 
  .$cell.type %>% unique %>% sort

gene_list <- panglaoDB_list_allMarkers[
  c("Pericytes",   "Astrocytes", "Bergmann_glia", "Microglia",
    "Oligodendrocyte_progenitor_cells", "Radial_glia_cells",
    "Satellite_glial_cells", "Tanycytes",
    "Neurons", "Immature_neurons", "Neural_stem/precursor_cells")
]

gene_list %>% length

panglaoDB_list_allMarkers %>% names

gene_list$Pericites_PMID34018679 <- c(
  "NODAL", "KCNJ8", "GGT1", "ABCC9", "PLXDC1", "SLC16A12", "SLC6A20",
  "HIGD1B", "EGFLAM", "ARHGAP42", "LIN7A", "S1PR3", "ECM2", "GPC3",
  "PTH1R", "NR1H3", "NDUFA4L2", "PDGFRB", "COX4I2", "GGT5", "GPER",
  "TMEM45A", "ENPEP", "IFITM1", "CD248", "CTHRC1", "ITGA4", "LAMC1",
  "ACE2", "MIR143HG", "FAM176B", "MYLIP", "RGS5", "MYO1B", "COLEC12",
  "TFPI", "HEYL", "OFD3L1", "CALD1", "SLC38A11", "GJC1", "KCNE4",
  "PLCE1", "LAMA2", "VASN", "LHFP", "EDNRA", "NID1", "GJA4", "F2R",
  "SEPT11", "SIK1", "NOTCH3", "PTEN", "EPS8", "ASPN", "PHLDB2",
  "CYTH3", "LAMC3", "MRC2", "TENC1", "PEAR1"
)

fgsea_out <- run_fgsea_all(mat, gene_list, nperm = 1000)


fgsea_out <- run_fgsea_all(mat, panglaoDB_list_allMarkers, nperm = 1000)


fgsea_out$first_per_cluster

fgsea_out$first_per_cluster %>% 
  arrange(cluster)


fgsea_out$results$cluster_2



# dostęp:
# fgsea_out$results$cluster_9
# fgsea_out$first_per_cluster

