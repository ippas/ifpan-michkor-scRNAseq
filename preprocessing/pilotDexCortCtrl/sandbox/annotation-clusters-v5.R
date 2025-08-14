
# Słownik pełnych nazw typów komórek
full_names <- c(
  "ODC"     = "Oligodendrocyte",
  "MG"      = "Microglia",
  "OPC"     = "Oligodendrocyte precursor cell",
  "INH"     = "Inhibitory neuron",
  "EX"      = "Excitatory neuron",
  "ASC"     = "Astrocyte",
  "PER.END" = "Pericyte / Endothelial cell"
)

# Ścieżki do danych referencyjnych GSE174367
path_data <- "~/projects/ifpan-michkor-scRNAseq/data/celltype_references/GSE174367"
file_h5   <- file.path(path_data, "GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5")
file_meta <- file.path(path_data, "GSE174367_snRNA-seq_cell_meta.csv.gz")

# Wywołanie funkcji
sce_ref <- prepare_expression_reference_singler_v1(
  expression_input = file_h5,       # plik z macierzą ekspresji
  meta_input       = file_meta,     # plik z metadanymi komórek
  label_col        = "Cell.Type",   # kolumna etykiet w meta
  full_names_map   = full_names,    # mapowanie skrótów na pełne nazwy
  assay_name       = "RNA",
  save_rds         = file.path("data/reference_celltype_expression_processed/", "GSE174367_reference_sce.rds")
)


# Przygotowanie SCE z danych pilotDexCortCtrl (query) dla SingleR
sce_pilotDexCortCtrl <- prepare_expression_query_singler_v1(
  seurat_obj = merged_pilotDexCortCtrl,        # Twój obiekt Seurat
  assay = "RNA",                               # assay z ekspresją
  data_layer_regex = "^data",                  # wybór warstwy data
  normalize_if_missing = TRUE,                 # normalizuj jeśli brak warstwy data
  gene_name_from_rownames = TRUE,              # konwersja nazw genów
  gene_sep = "-",                              # separator w nazwach genów
  agg_fun_for_duplicates = "mean",             # agregacja duplikatów
  cluster_col = "seurat_clusters_res0.17",     # kolumna z klastrami
  save_rds = "~/projects/ifpan-michkor-scRNAseq/data/reference_celltype_expression_processed/sce_pilotDexCortCtrl.rds",
  verbose = TRUE                               # logi kroków
)


res_pred <- predict_celltypes_singler_by_cluster_v1(
  sce_query = sce_pilotDexCortCtrl,
  sce_ref   = sce_ref,                   # z prepare_expression_reference_singler_v1
  label_col = "Cell.Type.Full",
  cluster_col = "cluster",
  min_common_genes = 1000,
  save_dir = NULL,
  verbose = TRUE
)

res_pred$annotation_table


# ##############################################################################
# ---- multi resolution predict ----
# ##############################################################################

# 1) Query z wieloma kolumnami klastrów
sce_pilotDexCortCtrl_multi <- prepare_expression_query_multi_singler_v1(
  seurat_obj = merged_pilotDexCortCtrl,
  assay = "RNA",
  data_layer_regex = "^data",
  normalize_if_missing = TRUE,
  gene_name_from_rownames = TRUE,
  gene_sep = "-",
  agg_fun_for_duplicates = "mean",
  cluster_cols = c( "seurat_clusters_res0.022", "seurat_clusters_res0.5", "seurat_clusters_res0.17"),
  alias_first_as = "cluster",     # SingleR-friendly alias
  save_rds = "data/reference_celltype_expression_processed//sce_query_multi.rds",
  verbose = TRUE
)



# 2) Predykcja dla każdej kolumny klastrów
res_multi <- predict_celltypes_singler_multi_v1(
  sce_query = sce_pilotDexCortCtrl_multi,
  sce_ref   = sce_ref,                   # z prepare_expression_reference_singler_v1
  label_col = "Cell.Type.Full",
  cluster_cols = c("seurat_clusters_res0.022", "seurat_clusters_res0.5", "seurat_clusters_res0.17"),
  min_common_genes = 1000,
  save_dir = "data/processed_reference_celltype_expression/singler_multi_results",
  verbose = TRUE
)



# ##############################################################################

# --------- helper: agregacja w rzadkiej macierzy po nazwach genów ----------
aggregate_sparse_matrix <- function(mat, groups, fun = c("mean","sum")) {
  stopifnot(inherits(mat, "dgCMatrix"))
  fun <- match.arg(fun)
  groups <- factor(groups, levels = unique(groups))
  G <- Matrix::sparse.model.matrix(~ groups - 1)
  summed <- Matrix::crossprod(G, mat)
  if (fun == "mean") {
    counts <- as.numeric(table(groups))
    summed <- summed / counts
  }
  rownames(summed) <- levels(groups)
  methods::as(summed, "dgCMatrix")
}

# =========================
# 1) Referencja z plików 10x + meta
# =========================
prepare_reference_from_10x <- function(
    file_h5,
    file_meta,
    label_col = "Cell.Type",
    full_names_map = NULL,       # named character vector: short -> full
    assay_name = "RNA",
    log_normalize_if_missing = TRUE,
    save_rds = NULL              # ścieżka do zapisu SCE (opcjonalnie)
) {
  stopifnot(file.exists(file_h5), file.exists(file_meta))
  counts_mat <- Read10X_h5(file_h5)
  seurat_ref <- CreateSeuratObject(counts = counts_mat, assay = assay_name)
  
  meta_df <- utils::read.csv(file_meta)
  stopifnot("Barcode" %in% colnames(meta_df))
  rownames(meta_df) <- meta_df$Barcode
  meta_df <- meta_df[colnames(seurat_ref), , drop = FALSE]
  seurat_ref <- AddMetaData(seurat_ref, meta_df)
  
  sce_ref <- as.SingleCellExperiment(seurat_ref)
  if (!label_col %in% colnames(colData(sce_ref))) {
    stop(sprintf("Kolumna etykiet '%s' nie występuje w meta referencji.", label_col))
  }
  
  # pełne nazwy (opcjonalnie)
  if (!is.null(full_names_map)) {
    sce_ref$Cell.Type.Full <- plyr::revalue(
      x = as.character(sce_ref[[label_col]]),
      replace = full_names_map,
      warn_missing = FALSE
    )
  } else {
    sce_ref$Cell.Type.Full <- as.character(sce_ref[[label_col]])
  }
  
  # logcounts
  if (log_normalize_if_missing && !"logcounts" %in% assayNames(sce_ref)) {
    sce_ref <- scuttle::logNormCounts(sce_ref)
  }
  
  if (!is.null(save_rds)) {
    saveRDS(sce_ref, file = save_rds)
  }
  sce_ref
}

# =========================
# 2) Przygotowanie danych testowych z Seurat v5
# =========================
prepare_test_from_seurat_v5 <- function(
    seurat_obj,
    assay = "RNA",
    data_layer_regex = "^data",           # które warstwy zebrać (zwykle "data")
    normalize_if_missing = TRUE,          # jeśli brak warstwy data -> NormalizeData()
    gene_name_from_rownames = TRUE,       # konwersja "ENSG...-SYMBOL" -> "SYMBOL"
    gene_sep = "-",
    agg_fun_for_duplicates = c("mean","sum"),
    cluster_col,                          # np. "seurat_clusters_res0.17"
    save_rds = NULL                       # opcjonalny zapis SCE
) {
  stopifnot(inherits(seurat_obj, "Seurat"))
  DefaultAssay(seurat_obj) <- assay
  agg_fun_for_duplicates <- match.arg(agg_fun_for_duplicates)
  
  # warstwy "data"
  data_layers <- grep(data_layer_regex, Layers(seurat_obj[[assay]]), value = TRUE)
  if (length(data_layers) == 0) {
    if (normalize_if_missing) {
      seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
      data_layers <- "data"
    } else {
      stop("Brak warstw 'data' i normalize_if_missing=FALSE.")
    }
  }
  
  # zlepienie warstw w jedną macierz GENES x CELLS
  mat_list <- lapply(
    data_layers,
    function(ly) LayerData(seurat_obj[[assay]], layer = ly)
  )
  logmat <- do.call(cbind, mat_list)
  logmat <- methods::as(logmat, "dgCMatrix")
  
  # nazwy genów -> symbole (opcjonalnie)
  if (gene_name_from_rownames) {
    rn <- rownames(logmat)
    # jeśli nazwa nie ma separatora, zostawiamy bez zmian
    gene_symbols <- ifelse(grepl(gene_sep, rn, fixed = TRUE),
                           sub(paste0("^[^", gene_sep, "]+", gene_sep), "", rn),
                           rn)
  } else {
    gene_symbols <- rownames(logmat)
  }
  
  dup_idx <- duplicated(gene_symbols) | duplicated(gene_symbols, fromLast = TRUE)
  mat_nondup <- logmat[!dup_idx, , drop = FALSE]
  rownames(mat_nondup) <- gene_symbols[!dup_idx]
  
  if (any(dup_idx)) {
    mat_dup <- logmat[dup_idx, , drop = FALSE]
    mat_dup_mean <- aggregate_sparse_matrix(mat_dup, groups = gene_symbols[dup_idx], fun = agg_fun_for_duplicates)
    logmat_avg <- rbind(mat_nondup, mat_dup_mean)
  } else {
    logmat_avg <- mat_nondup
  }
  logmat_avg <- methods::as(logmat_avg, "dgCMatrix")
  
  # SCE + klastry
  sce_test <- SingleCellExperiment(assays = list(logcounts = logmat_avg))
  
  if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
    stop(sprintf("Kolumna klastrów '%s' nie występuje w meta Seurat.", cluster_col))
  }
  if (!all(colnames(sce_test) %in% colnames(seurat_obj))) {
    stop("Kolumny SCE nie pokrywają się z kolumnami (komórkami) w obiekcie Seurat.")
  }
  sce_test$cluster <- seurat_obj[[cluster_col]][colnames(sce_test)]
  
  if (!is.null(save_rds)) saveRDS(sce_test, file = save_rds)
  sce_test
}

# =========================
# 3) Uruchomienie SingleR (per klaster) + raport
# =========================
run_singler_by_cluster <- function(
    sce_test,
    sce_ref,
    label_col = "Cell.Type.Full",
    cluster_col = "cluster",
    min_common_genes = 1000,
    assay.type.test = "logcounts",
    assay.type.ref  = "logcounts",
    save_dir = NULL      # jeśli podasz ścieżkę: zapisze predykcje i tabelę podsum.
) {
  if (!label_col %in% colnames(colData(sce_ref))) {
    stop(sprintf("W referencji brak kolumny etykiet '%s'.", label_col))
  }
  if (!cluster_col %in% colnames(colData(sce_test))) {
    stop(sprintf("W teście brak kolumny klastrów '%s'.", cluster_col))
  }
  
  common <- intersect(rownames(sce_ref), rownames(sce_test))
  if (length(common) < min_common_genes) {
    stop(sprintf("Zbyt mało wspólnych genów: %d (< %d).", length(common), min_common_genes))
  }
  
  sce_ref2  <- sce_ref [common, ]
  sce_test2 <- sce_test[common, ]
  
  pred <- SingleR(
    test            = sce_test2,
    ref             = sce_ref2,
    labels          = sce_ref2[[label_col]],
    clusters        = sce_test2[[cluster_col]],
    assay.type.test = assay.type.test,
    assay.type.ref  = assay.type.ref
  )
  
  # podsumowanie per klaster
  annot_tbl <- data.frame(
    cluster    = rownames(pred),
    label      = pred$labels,
    best_score = apply(pred$scores, 1, max),
    row.names  = NULL
  )
  
  if (!is.null(save_dir)) {
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
    saveRDS(pred, file = file.path(save_dir, "singler_predictions.rds"))
    utils::write.csv(annot_tbl, file = file.path(save_dir, "singler_annotation_table.csv"), row.names = FALSE)
  }
  
  list(pred = pred, annotation_table = annot_tbl)
}

# =========================
# 4) Dopięcie etykiet do obiektu Seurat
# =========================
attach_predictions_to_seurat <- function(
    seurat_obj,
    pred,                          # wynik SingleR (obiekt 'SingleR')
    cluster_col,                   # ta sama kolumna co użyta w prepare_test...
    new_meta_col = "celltype_pred",
    save_rds = NULL
) {
  stopifnot(inherits(seurat_obj, "Seurat"))
  if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
    stop(sprintf("W meta Seurat brak kolumny klastrów '%s'.", cluster_col))
  }
  # mapowanie: wartości klastrów -> etykiety
  cl_to_label <- setNames(pred$labels, rownames(pred))
  
  seurat_obj[[new_meta_col]] <- plyr::mapvalues(
    seurat_obj[[cluster_col]],
    from = names(cl_to_label),
    to   = unname(cl_to_label),
    warn_missing = FALSE
  )
  
  if (!is.null(save_rds)) saveRDS(seurat_obj, file = save_rds)
  seurat_obj
}
