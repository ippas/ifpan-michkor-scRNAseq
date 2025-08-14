full_names <- c(
  "ODC"     = "Oligodendrocyte",
  "MG"      = "Microglia",
  "OPC"     = "Oligodendrocyte precursor cell",
  "INH"     = "Inhibitory neuron",
  "EX"      = "Excitatory neuron",
  "ASC"     = "Astrocyte",
  "PER.END" = "Pericyte / Endothelial cell"
)


# 1. Ścieżki do plików
path_data <- "~/projects/ifpan-michkor-scRNAseq/data/celltype_references/GSE174367"
file_h5   <- file.path(path_data, "GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5")
file_meta <- file.path(path_data, "GSE174367_snRNA-seq_cell_meta.csv.gz")


# 2. Wczytanie macierzy ekspresji
counts_mat <- Read10X_h5(file_h5)
seurat_obj <- CreateSeuratObject(counts = counts_mat)


# 3. Wczytanie metadanych
meta_df <- read.csv(file_meta)
rownames(meta_df) <- meta_df$Barcode


# 4. Dopasowanie metadanych do obiektu Seurat
meta_df <- meta_df[colnames(seurat_obj), , drop = FALSE]
seurat_obj <- AddMetaData(seurat_obj, meta_df)

# 5. Konwersja do SingleCellExperiment
sce <- as.SingleCellExperiment(seurat_obj)

sce@metadata

sce$Cell.Type.Full <- recode(sce$Cell.Type, !!!full_names)

table(sce$Cell.Type.Full)

# saveRDS(sce, file.path(path_data, "GSE174367_reference_sce.rds"))


# ##############################################################################



# --- 0) Referencja: 'sce' już wczytane z GSE174367 ---
# Załóżmy, że masz: sce$Cell.Type.Full z pełnymi nazwami
stopifnot("Cell.Type.Full" %in% colnames(colData(sce)))

# --- 1) Przygotuj dane testowe z Seurat v5 ---
obj <- merged_pilotDexCortCtrl
DefaultAssay(obj) <- "RNA"

# Połącz wszystkie warstwy 'data.*' (log-normalized)
data_layers <- grep("^data", Layers(obj[["RNA"]]), value = TRUE)
if (length(data_layers) == 0) {
  obj <- NormalizeData(obj, verbose = FALSE)
  data_layers <- "data"
}
mat_list <- lapply(data_layers, function(ly) LayerData(obj[["RNA"]], layer = ly))
logmat   <- do.call(cbind, mat_list)                  # GENES x CELLS
logmat   <- as(logmat, "dgCMatrix")

# Zamień "ENSG...-SYMBOL" -> "SYMBOL" i zrób średnią dla duplikatów symboli
gene_symbols <- sub("^[^-]+-", "", rownames(logmat))
dup_mask <- duplicated(gene_symbols) | duplicated(gene_symbols, fromLast = TRUE)
if (any(dup_mask)) {
  grp <- factor(gene_symbols)
  logmat_sum <- rowsum(logmat, group = grp)                 # SYMBOL x CELLS
  cnt <- as.numeric(table(grp))
  logmat_avg <- logmat_sum / cnt
} else {
  logmat_avg <- logmat
  rownames(logmat_avg) <- gene_symbols
}


# ##############################################################################


aggregate_sparse_matrix <- function(mat, groups, fun = c("sum","mean")) {
  fun <- match.arg(fun)
  groups <- factor(groups, levels = unique(groups))        # długość = nrow(mat)
  
  # Macierz designu: (n_genes x n_groups)
  G <- sparse.model.matrix(~ groups - 1)                   # kolumny = poziomy grup
  
  # Agregacja po wierszach: (n_groups x n_cells)
  summed <- Matrix::crossprod(G, mat)                      # t(G) %*% mat
  
  if (fun == "mean") {
    counts <- as.numeric(table(groups))
    summed <- summed / counts
  }
  rownames(summed) <- levels(groups)
  as(summed, "dgCMatrix")
}

# --- użycie ---
gene_symbols <- sub("^[^-]+-", "", rownames(logmat))
dup_idx <- duplicated(gene_symbols) | duplicated(gene_symbols, fromLast = TRUE)

# 1) nieduplikaty
mat_nondup <- logmat[!dup_idx, , drop = FALSE]
rownames(mat_nondup) <- gene_symbols[!dup_idx]

# 2) duplikaty (średnia po symbolu)
mat_dup_mean <- aggregate_sparse_matrix(
  logmat[dup_idx, , drop = FALSE],
  groups = gene_symbols[dup_idx],
  fun = "mean"
)

# 3) połącz
logmat_avg <- rbind(mat_nondup, mat_dup_mean)
logmat_avg <- as(logmat_avg, "dgCMatrix")

# ##############################################################################


# 1) Zbuduj obiekt testowy SCE z Twoich danych (logmat_avg = GENES x CELLS)
sce_test <- SingleCellExperiment(assays = list(logcounts = as(logmat_avg, "dgCMatrix")))
colData(sce_test)$seurat_clusters <- merged_pilotDexCortCtrl$seurat_clusters[colnames(sce_test)]

# 2) Ujednolić geny z referencją GSE174367 (sce)
stopifnot("Cell.Type.Full" %in% colnames(colData(sce)))
common <- intersect(rownames(sce), rownames(sce_test))
length(common)  # kontrolnie, powinno być dużo (tysiące)

sce_ref   <- sce[common, ]
sce_test2 <- sce_test[common, ]

# 3) SingleR per klaster
pred <- SingleR(
  test     = sce_test2,
  ref      = sce_ref,
  labels   = sce_ref$Cell.Type.Full,
  clusters = sce_test2$seurat_clusters
)

# 4) Dopnij etykiety do Seurat
merged_pilotDexCortCtrl$celltype_GSE174367 <- plyr::mapvalues(
  merged_pilotDexCortCtrl$seurat_clusters,
  from = rownames(pred),
  to   = pred$labels
)

# 5) Szybki podgląd jakości i zapis wyników
print(table(merged_pilotDexCortCtrl$celltype_GSE174367, useNA = "ifany"))

# (opcjonalnie) mapa „pewności” przypisań
plotScoreHeatmap(pred)

# zapisz tabelę klaster → typ + „pewność”
annot_tbl <- data.frame(
  cluster    = rownames(pred),
  label      = pred$labels,
  best_score = apply(pred$scores, 1, max),
  stringsAsFactors = FALSE
)


# ##############################################################################
if (!requireNamespace("scuttle", quietly = TRUE)) BiocManager::install("scuttle")
library(scuttle)

# dodaj logcounts do referencji i testu
sce_ref   <- scuttle::logNormCounts(sce_ref)    # doda assay "logcounts"
sce_test2 <- scuttle::logNormCounts(sce_test2)  # doda assay "logcounts"

assayNames(sce_ref)    # sprawdź, że jest "logcounts"
assayNames(sce_test2)

# uruchom SingleR
pred <- SingleR(
  test     = sce_test2,
  ref      = sce_ref,
  labels   = sce_ref$Cell.Type.Full,
  clusters = sce_test2$seurat_clusters,
  assay.type.test = "logcounts",
  assay.type.ref  = "logcounts"
)

# Zliczenie przypisanych typów komórek
table(pred$labels)
# ##############################################################################




# 1) Predykcja PER KOMÓRKA
pred_cells <- SingleR(
  test            = sce_test2,
  ref             = sce_ref,
  labels          = sce_ref$Cell.Type.Full,
  assay.type.test = "logcounts",
  assay.type.ref  = "logcounts"
)

# 2) Szybki podgląd (to są LICZBY KOMÓREK)
table(pred_cells$labels)

# 3) Dopnij do Seurat (po nazwach komórek)
merged_pilotDexCortCtrl$celltype_GSE174367_cell <- pred_cells$labels[colnames(sce_test2)]

# 4) (opcjonalnie) bardziej konserwatywne etykiety – tylko pewne
merged_pilotDexCortCtrl$celltype_GSE174367_pruned <- pred_cells$pruned.labels[colnames(sce_test2)]

# 5) (opcjonalnie) UMAP z etykietami
DimPlot(merged_pilotDexCortCtrl, group.by = "celltype_GSE174367_cell",
        label = TRUE, repel = TRUE) + NoLegend()

# 6) (opcjonalnie) heatmapa „pewności” na próbce 500 komórek
plotScoreHeatmap(pred_cells[sample(seq_len(nrow(pred_cells)), min(500, nrow(pred_cells))), ])