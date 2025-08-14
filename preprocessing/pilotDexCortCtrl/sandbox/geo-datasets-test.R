if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("SeuratDisk", quietly = TRUE)) install.packages("SeuratDisk")
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
  BiocManager::install("SingleCellExperiment")
}


install.packages("hdf5r")
library(hdf5r)


# Ścieżki do plików
path_data <- "~/projects/ifpan-michkor-scRNAseq/data/celltype_references/GSE174367"
file_h5   <- file.path(path_data, "GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5")
file_meta <- file.path(path_data, "GSE174367_snRNA-seq_cell_meta.csv.gz")




# 1. Wczytanie macierzy ekspresji
sce_seurat <- Read10X_h5(file_h5)

# Jeśli w pliku jest kilka assayów, możesz wybrać jeden:
sce_seurat <- CreateSeuratObject(counts = sce_seurat)

# 2. Wczytanie metadanych
meta_df <- read.csv(file_meta)

# 3. Dopasowanie metadanych do komórek
# Zakładam, że w meta_df jest kolumna 'barcode' lub podobna z ID komórek
# i że odpowiada ona kolumnom w Seurat object
barcode_col <- "barcode" # zmień, jeśli kolumna nazywa się inaczej
meta_df <- as.data.frame(meta_df)
rownames(meta_df) <- meta_df[[barcode_col]]

# Upewniamy się, że kolejność metadanych odpowiada komórkom w Seurat
meta_df <- meta_df[colnames(sce_seurat), , drop = FALSE]

# 4. Dodanie metadanych do obiektu Seurat
sce_seurat <- AddMetaData(sce_seurat, metadata = meta_df)

# 5. Konwersja do SingleCellExperiment (opcjonalne — wymagane przez SingleR)
sce <- as.SingleCellExperiment(sce_seurat)

# Teraz w 'sce' masz dane ekspresji + metadane gotowe do filtrowania
sce

"Cell.Type"
