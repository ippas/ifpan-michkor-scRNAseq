install.packages("remotes")  # jeśli nie masz
remotes::install_github("satijalab/seurat", ref = "release/5.0.0")


library(Matrix)
library(dplyr)
library(Seurat)
# ##############################################################################
# ---- read data ----
# ##############################################################################
# Ścieżka do folderu DGE_filtered
data_dir <- "~/projects/ifpan-michkor-scRNAseq/data/pilotDexCortCtrl/BG60BLUE_L6-out/all-well/DGE_filtered"


counts <- readMM(file.path(data_dir, "count_matrix.mtx"))

# Wczytanie nazw genów
genes <- read.csv(file.path(data_dir, "all_genes.csv"))
rownames(counts) <- genes$gene_name[1:nrow(counts)]  # albo gene_id

# Nadanie prostych nazw komórkom
colnames(counts) <- paste0("cell_", seq_len(ncol(counts)))

# Stworzenie obiektu Seurat bez metadanych
seurat_obj <- CreateSeuratObject(counts = counts, project = "BG60_simple")


# Jeśli nie masz jeszcze normalization
seurat_obj <- NormalizeData(seurat_obj)

# Wybór najbardziej zmiennych genów (anchorów)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# (opcjonalnie) podejrzyj top geny
top10 <- head(VariableFeatures(seurat_obj), 10)
print(top10)

# Skalowanie danych (na anchorach)
seurat_obj <- ScaleData(seurat_obj)

# PCA tylko na tych zmiennych genach
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))

# Redukcja wymiarów – UMAP
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

# Sąsiedztwa i klasteryzacja
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# UMAP z klastrami
DimPlot(seurat_obj, label = TRUE)

# Prosty pipeline
# seurat_obj <- NormalizeData(seurat_obj)
# seurat_obj <- FindVariableFeatures(seurat_obj)
# seurat_obj <- ScaleData(seurat_obj)
# seurat_obj <- RunPCA(seurat_obj)
# seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
# seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10, graph.name = NULL, method = "rann", verbose = T)

seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Wizualizacja
DimPlot(seurat_obj, label = TRUE)


DimPlot(seurat_obj, reduction = "umap", label = TRUE)



# ##############################################################################
subset_obj <- subset(seurat_obj, cells = sample(colnames(seurat_obj), 3000))

subset_obj <- NormalizeData(subset_obj)
subset_obj <- FindVariableFeatures(subset_obj)
subset_obj <- ScaleData(subset_obj)
subset_obj <- RunPCA(subset_obj)
subset_obj <- RunUMAP(subset_obj, dims = 1:20)
subset_obj <- FindNeighbors(subset_obj, dims = 1:20)  # teraz powinno zadziałać
subset_obj <- FindClusters(subset_obj)
