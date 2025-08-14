# 1) instalacja jednorazowo (jeśli nie masz)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("SingleCellExperiment","SummarizedExperiment"))

# install.packages(c("Seurat","SingleR","celldex","SummarizedExperiment","plyr"))  # jeśli potrzeba
library(Seurat) 
library(SingleR) 
library(celldex) 
library(SummarizedExperiment) 
library(plyr)
library(SingleCellExperiment)

obj <- merged_pilotDexCortCtrl
DefaultAssay(obj) <- "RNA"

# 1) Połącz wszystkie warstwy log-normalizowane
data_layers <- grep("^data", Layers(obj[["RNA"]]), value = TRUE)
mat_list <- lapply(data_layers, function(ly) LayerData(obj[["RNA"]], layer = ly))
logmat <- do.call(cbind, mat_list)
logmat <- as(logmat, "dgCMatrix")

# 2) Wyciągnij same symbole
gene_symbols <- sub("^[^-]+-", "", rownames(logmat))

# 3) Identyfikuj duplikaty
dup_symbols <- gene_symbols[duplicated(gene_symbols) | duplicated(gene_symbols, fromLast = TRUE)]

if (length(dup_symbols) > 0) {
  message("Liczba genów z duplikatami symboli: ", length(unique(dup_symbols)))
  
  # Agreguj tylko duplikaty
  dup_idx <- gene_symbols %in% dup_symbols
  nondup_idx <- !dup_idx
  
  # rowsum dla duplikatów
  dup_groups <- factor(gene_symbols[dup_idx])
  logmat_dup_avg <- rowsum(logmat[dup_idx, ], group = dup_groups) / as.numeric(table(dup_groups))
  
  # Połącz z nieduplikatami
  logmat_avg <- rbind(logmat[nondup_idx, ], logmat_dup_avg)
  
} else {
  message("Brak duplikatów symboli — używam oryginalnej macierzy")
  logmat_avg <- logmat
}

# 4) Ustaw czyste symbole jako rownames
rownames(logmat_avg) <- c(gene_symbols[!gene_symbols %in% dup_symbols],
                          unique(dup_symbols))

# 5) Zbuduj SCE
sce <- SingleCellExperiment(
  assays = list(logcounts = as(logmat_avg, "dgCMatrix")),
  colData = DataFrame(seurat_clusters = obj$seurat_clusters[colnames(logmat_avg)])
)

# 6) Referencja HPCA do mózgu
hpca <- HumanPrimaryCellAtlasData()
brain_idx <- grepl("astro|neur|oligo|microglia|ependym|neural",
                   hpca$label.main, ignore.case = TRUE)
hpca_brain <- hpca[, brain_idx]

# 7) Sprawdź wspólne geny
message("Wspólnych genów: ", length(intersect(rownames(sce), rownames(hpca_brain))))

# 8) SingleR na klastrach
pred <- SingleR(
  test     = sce,
  ref      = hpca_brain,
  labels   = hpca_brain$label.main,
  clusters = sce$seurat_clusters
)

# 9) Przypisz etykiety
obj$celltype_HPCA_brain <- mapvalues(
  obj$seurat_clusters,
  from = rownames(pred),
  to   = pred$labels
)
