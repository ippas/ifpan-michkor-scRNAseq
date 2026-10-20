load("data/reference_celltype_expression_processed/sce_ref_MTGRNAseqAllNuclei2022.RData")


GSE174367_reference_sce <- readRDS("data/celltype_references/GSE174367/GSE174367_reference_sce.rds")



GSE174367_reference_sce %>% str
# ##############################################################################
# ---- functions ----
# ##############################################################################
make_top_expr_signatures <- function(
    sce,
    label_col    = "subclass_label",
    assay_name   = "logcounts",
    top_n        = c(200),        # może być wektor, np. c(50,100,200,500)
    filter_genes = NULL,          # np. protein-coding autosomy+X/Y
    prefix       = "",            # prefiks dopisywany do nazw typów (np. "ATLAS_")
    verbose      = TRUE
) {
  stopifnot(label_col %in% colnames(colData(sce)))
  mat    <- SummarizedExperiment::assay(sce, assay_name)
  labels <- colData(sce)[[label_col]]
  
  # diagnostyka komórek
  n_total <- ncol(mat)
  keep <- !is.na(labels)
  n_removed <- sum(!keep)
  mat    <- mat[, keep, drop = FALSE]
  labels <- droplevels(labels[keep])
  n_final <- ncol(mat)
  if (verbose) {
    message("Cells: total = ", n_total,
            ", removed_NA = ", n_removed,
            ", remaining = ", n_final)
  }
  
  # filtr genów (opcjonalny)
  if (!is.null(filter_genes)) {
    keep_genes <- intersect(rownames(mat), filter_genes)
    if (verbose) {
      message("Genes: total = ", nrow(mat),
              ", after filtering = ", length(keep_genes))
    }
    mat <- mat[keep_genes, , drop = FALSE]
  } else if (verbose) {
    message("Genes: total = ", nrow(mat), " (no gene filter applied)")
  }
  
  # policz średnie ekspresje per cell type tylko raz
  ct_levels <- levels(labels)
  avg_by_ct <- lapply(ct_levels, function(ct) {
    cells <- which(labels == ct)
    if (length(cells) == 0) return(NULL)
    Matrix::rowMeans(mat[, cells, drop = FALSE])
  })
  names(avg_by_ct) <- ct_levels
  avg_by_ct <- avg_by_ct[!vapply(avg_by_ct, is.null, logical(1))]
  
  # zastosuj prefiks do nazw typów
  ct_names_with_prefix <- paste0(prefix, names(avg_by_ct))
  names(avg_by_ct) <- ct_names_with_prefix
  
  # helper: top-N dla jednego cell type
  top_for_ct <- function(avg_vec, N) {
    nm <- names(sort(avg_vec, decreasing = TRUE))
    head(nm, min(N, length(nm)))
  }
  
  # pojedyncze top_n → lista {prefixed_ct → geny}
  if (length(top_n) == 1) {
    N <- top_n[1]
    signatures <- lapply(avg_by_ct, top_for_ct, N = N)
    return(signatures)
  }
  
  # wiele progów → lista poziomu wyżej
  signatures_multi <- lapply(top_n, function(N) {
    lapply(avg_by_ct, top_for_ct, N = N)
  })
  names(signatures_multi) <- paste0("Top", top_n)
  
  signatures_multi
}

make_singatures_avgNormTau <- function(
    sce,
    label_col    = "subclass_label",
    assay_name   = "logcounts",
    filter_genes = NULL,
    top_n        = c(100),   # może być wektor, np. c(50,100,200,500)
    min_cells    = 1,
    pseudocount  = 0,
    prefix       = "",       # <<< nowy argument
    verbose      = TRUE
) {
  stopifnot(label_col %in% colnames(colData(sce)))
  mat    <- SummarizedExperiment::assay(sce, assay_name)
  labels <- colData(sce)[[label_col]]
  
  # --- komórki: usuwanie NA + małe typy ---
  n_total <- ncol(mat)
  keep_lab <- !is.na(labels)
  mat    <- mat[, keep_lab, drop = FALSE]
  labels <- droplevels(labels[keep_lab])
  
  if (min_cells > 1) {
    tab_ct <- table(labels)
    ok_ct  <- names(tab_ct)[tab_ct >= min_cells]
    keep_ct <- labels %in% ok_ct
    mat    <- mat[, keep_ct, drop = FALSE]
    labels <- droplevels(labels[keep_ct])
  }
  
  if (verbose) {
    removed_na    <- sum(!keep_lab)
    removed_small <- if (exists("keep_ct")) sum(!keep_ct) else 0
    message("Cells: total=", n_total,
            ", removed_NA=", removed_na,
            ", removed_smallTypes=", removed_small,
            ", remaining=", ncol(mat))
  }
  
  # --- geny: filtracja ---
  if (!is.null(filter_genes)) {
    keep_genes <- intersect(rownames(mat), filter_genes)
    if (verbose) message("Genes: total=", nrow(mat),
                         ", after_filter=", length(keep_genes))
    mat <- mat[keep_genes, , drop = FALSE]
  } else if (verbose) {
    message("Genes: total=", nrow(mat), " (no gene filter applied)")
  }
  
  # --- średnia ekspresja per typ ---
  lvl <- levels(labels)
  avg_list <- vector("list", length(lvl)); names(avg_list) <- lvl
  for (ct in lvl) {
    idx <- which(labels == ct)
    if (length(idx) == 0) next
    avg_list[[ct]] <- Matrix::rowMeans(mat[, idx, drop = FALSE]) + pseudocount
  }
  keep_nonnull <- !vapply(avg_list, is.null, logical(1))
  avg_expr <- do.call(cbind, avg_list[keep_nonnull])
  
  # dopisz prefix do nazw kolumn (cell types)
  colnames(avg_expr) <- paste0(prefix, names(avg_list)[keep_nonnull])
  rownames(avg_expr) <- rownames(mat)
  
  if (verbose) message("Built avg_expr: ", nrow(avg_expr), " genes × ", ncol(avg_expr), " cell types")
  
  # --- normalizacja per gen ---
  row_max <- apply(avg_expr, 1, max, na.rm = TRUE)
  row_max[row_max == 0 | is.na(row_max)] <- 1
  norm_expr <- sweep(avg_expr, 1, row_max, "/")
  
  # --- tau ---
  k <- ncol(norm_expr)
  tau <- apply(norm_expr, 1, function(x) {
    if (all(x == 0)) return(0)
    sum(1 - x / max(x)) / (k - 1)
  })
  
  # --- score specyficzności ---
  spec_score <- norm_expr * tau
  
  # --- top-N genów dla wielu progów ---
  top_gene_sets <- lapply(top_n, function(N) {
    lapply(seq_len(ncol(spec_score)), function(j) {
      vals <- spec_score[, j]
      nm   <- names(sort(vals, decreasing = TRUE))
      head(nm, min(N, length(nm)))
    }) |> `names<-`(colnames(spec_score))
  })
  names(top_gene_sets) <- paste0("Top", top_n)
  
  list(
    avg_expr   = avg_expr,
    norm_expr  = norm_expr,
    gene_specificity_tau = tau,
    specificity_score_matrix = spec_score,
    top_marker_genes = top_gene_sets,
    diagnostics = list(
      n_cells_total   = n_total,
      n_cells_final   = ncol(mat),
      n_celltypes     = ncol(avg_expr),
      n_genes_final   = nrow(avg_expr),
      label_col       = label_col,
      assay_used      = assay_name,
      top_n           = top_n,
      prefix          = prefix
    )
  )
}

# ##############################################################################
# ---- prepare data ----
# ##############################################################################
proteinCodingGenes_autosomalSex_v112 <- download_biomart_genes(version = "v112") %>% 
  filter(gene_biotype == "protein_coding") %>% 
  filter(chromosome_name %in% c("1", "2", "3", "4", "5", "6", "7", "8","9", "10", 
                                "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22",
                                "Y", "X")) %>% 
  filter(external_gene_name != "") %>% 
  .$external_gene_name %>% unique

proteinCodingGenes_autosomalSex_v112



# ##############################################################################
# ---- analysis ----
# ##############################################################################
signatture_atals_topExpr <- make_top_expr_signatures(
  sce_ref_MTGRNAseqAllNuclei2022,
  label_col = "subclass_label",
  assay_name = "logcounts",
  top_n = c(50, 100, 200, 500),
  filter_genes = proteinCodingGenes_autosomalSex_v112   # <- podajesz wektor genów
)

signatture_atals_topExprTop200$Astrocyte


signature_atlas_avgNormTau <- make_singatures_avgNormTau(
  sce           = sce_ref_MTGRNAseqAllNuclei2022,
  label_col     = "subclass_label",
  assay_name    = "logcounts",
  filter_genes  = proteinCodingGenes_autosomalSex_v112,
  top_n         = c(50, 100, 200, 500),
  min_cells     = 1,
  pseudocount   = 0,
  prefix = "MTGRNAseqAllNuclei2022_",
  verbose       = TRUE
)


GSE174367_reference_sce$Cell.Type <- factor(GSE174367_reference_sce$Cell.Type)

signature_gse_avgNormTau <- make_singatures_avgNormTau(
  sce           = GSE174367_reference_sce,
  label_col     = "Cell.Type",     # albo "Cell.Type.Full"
  assay_name    = "logcounts",
  filter_genes  = proteinCodingGenes_autosomalSex_v112,
  top_n         = c(50,100,200,500),
  min_cells     = 10,              # np. odrzucaj typy z <10 komórek
  pseudocount   = 0,
  prefix        = "GSE174367_",    # prefiks żeby łatwo odróżniać
  verbose       = TRUE
)


rm(GSE174367_reference_sce)
rm(sce_ref_MTGRNAseqAllNuclei2022)


names(panglaoDB_list_allMarkers) <- paste0("PanglaoDB_", names(panglaoDB_list_allMarkers))

c(
  as.list(signature_gse_avgNormTau$top_marker_genes$Top100),
  as.list(panglaoDB_list_allMarkers),
  as.list(signature_atlas_avgNormTau$top_marker_genes$Top100)
)  %>% 
  (\(x) x[sapply(x, length) >= 5])() -> allGeneSignatures_panglaoDB_GSE174367_MTGRNAseqAllNuclei2022

gsea_results_panglaoDB_GSE174367_MTGRNAseqAllNuclei2022 <- run_celltype_fgseaEnrichemnt_multClusters(
  seurat_obj   = integrated_pilotDexCortCtrl,
  biomart_version = "112",
  cluster_cols = c(
    "seurat_clusters_res0.04",
    "seurat_clusters_res0.13",
    "seurat_clusters_res0.2",
    "seurat_clusters_res0.5"
  ),
  gene_lists =  allGeneSignatures_panglaoDB_GSE174367_MTGRNAseqAllNuclei2022,
  nperm        = 1000,
  do_fisher    = TRUE,
  export_xlsx  = TRUE,
  # select_samples = c(
  #   "BG60BLUE_A1", "BG60BLUE_A2", "BG60BLUE_A3", "BG60BLUE_A4",
  #   "BG60GREEN_A1", "BG60GREEN_A2", "BG60GREEN_A3", "BG60GREEN_A4"
  # ),
  xlsx_path    = "data/pilotDexCortCtrl/r_outputs/gsea_results/gsea_PanglaoDB_GSE174367_MTGRNAseqAllNuclei2022_{cluster}.xlsx"
)

save(gsea_results_panglaoDB_GSE174367_MTGRNAseqAllNuclei2022, "data/pilotDexCortCtrl/r_outputs/gsea_results/gsea_results_panglaoDB_GSE174367_MTGRNAseqAllNuclei2022.Rdata")


download_biomart_genes(version = "v112")

res_multi_aprioriSelectedCellType2$seurat_clusters_res0.13$fgsea$gene_lists_used

res_multi_aprioriSelectedCellType2$seurat_clusters_res0.13$fgsea$results$cluster_0

