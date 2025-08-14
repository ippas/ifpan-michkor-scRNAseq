# =========================
# 1) Słownik pełnych nazw
# =========================
full_names <- c(
  "ODC"     = "Oligodendrocyte",
  "MG"      = "Microglia",
  "OPC"     = "Oligodendrocyte precursor cell",
  "INH"     = "Inhibitory neuron",
  "EX"      = "Excitatory neuron",
  "ASC"     = "Astrocyte",
  "PER.END" = "Pericyte / Endothelial cell"
)

# =========================
# 2) Referencja GSE174367
# =========================
path_data <- "~/projects/ifpan-michkor-scRNAseq/data/celltype_references/GSE174367"
file_h5   <- file.path(path_data, "GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5")
file_meta <- file.path(path_data, "GSE174367_snRNA-seq_cell_meta.csv.gz")

counts_mat <- Read10X_h5(file_h5)
seurat_ref <- CreateSeuratObject(counts = counts_mat, assay = "RNA")

meta_df <- read.csv(file_meta)
rownames(meta_df) <- meta_df$Barcode
meta_df <- meta_df[colnames(seurat_ref), , drop = FALSE]
seurat_ref <- AddMetaData(seurat_ref, meta_df)

# -> SCE + pełne etykiety
sce_ref <- as.SingleCellExperiment(seurat_ref)
stopifnot("Cell.Type" %in% colnames(colData(sce_ref)))
sce_ref$Cell.Type.Full <- plyr::revalue(sce_ref$Cell.Type, full_names, warn_missing = FALSE)

# log-normalizacja referencji (jeśli nie ma logcounts)
if (!"logcounts" %in% assayNames(sce_ref)) {
  sce_ref <- scuttle::logNormCounts(sce_ref)
}

# =========================
# 3) TEST: przygotowanie z Seurat v5 (obj = merged_pilotDexCortCtrl)
# =========================
obj <- merged_pilotDexCortCtrl
DefaultAssay(obj) <- "RNA"

# a) Pobierz warstwę zlogowanej ekspresji (jeśli brak, znormalizuj)
data_layers <- grep("^data", Layers(obj[["RNA"]]), value = TRUE)
if (length(data_layers) == 0) {
  obj <- NormalizeData(obj, verbose = FALSE)
  data_layers <- "data"
}
# połącz warstwy w jedną macierz GENES x CELLS
mat_list <- lapply(data_layers, function(ly) LayerData(obj[["RNA"]], layer = ly))
logmat   <- do.call(cbind, mat_list)
logmat   <- as(logmat, "dgCMatrix")

# b) Zamień "ENSG...-SYMBOL" -> "SYMBOL" i uśrednij duplikaty
aggregate_sparse_matrix <- function(mat, groups, fun = c("sum","mean")) {
  fun <- match.arg(fun)
  groups <- factor(groups, levels = unique(groups))
  G <- sparse.model.matrix(~ groups - 1)
  summed <- Matrix::crossprod(G, mat)
  if (fun == "mean") {
    counts <- as.numeric(table(groups))
    summed <- summed / counts
  }
  rownames(summed) <- levels(groups)
  as(summed, "dgCMatrix")
}

gene_symbols <- sub("^[^-]+-", "", rownames(logmat))
dup_idx      <- duplicated(gene_symbols) | duplicated(gene_symbols, fromLast = TRUE)

mat_nondup <- logmat[!dup_idx, , drop = FALSE]
rownames(mat_nondup) <- gene_symbols[!dup_idx]

mat_dup_mean <- aggregate_sparse_matrix(
  logmat[dup_idx, , drop = FALSE],
  groups = gene_symbols[dup_idx],
  fun = "mean"
)

logmat_avg <- rbind(mat_nondup, mat_dup_mean)
logmat_avg <- as(logmat_avg, "dgCMatrix")

# c) SCE testowy + klastry
sce_test <- SingleCellExperiment(assays = list(logcounts = logmat_avg))
stopifnot("seurat_clusters_res0.022" %in% colnames(obj@meta.data))
sce_test$seurat_clusters <- obj$seurat_clusters_res0.17[colnames(sce_test)]

# =========================
# 4) Ujednolicenie genów i SingleR per klaster
# =========================
common <- intersect(rownames(sce_ref), rownames(sce_test))
stopifnot(length(common) > 1000)  # sanity check

sce_ref2  <- sce_ref [common, ]
sce_test2 <- sce_test[common, ]

pred <- SingleR(
  test            = sce_test2,
  ref             = sce_ref2,
  labels          = sce_ref2$Cell.Type.Full,
  clusters        = sce_test2$seurat_clusters,
  assay.type.test = "logcounts",
  assay.type.ref  = "logcounts"
)

# =========================
# 5) Dopinanie etykiet do Seurat + szybki podgląd
# =========================
# pred wiersze = nazwy klastrów
obj$celltype_GSE174367 <- plyr::mapvalues(
  obj$seurat_clusters_res0.022,
  from = rownames(pred),
  to   = pred$labels,
  warn_missing = FALSE
)

# Podgląd i prosta tabela jakości
print(table(obj$celltype_GSE174367, useNA = "ifany"))

annot_tbl <- data.frame(
  cluster    = rownames(pred),
  label      = pred$labels,
  best_score = apply(pred$scores, 1, max),
  row.names  = NULL
)
annot_tbl[order(annot_tbl$cluster), ]
