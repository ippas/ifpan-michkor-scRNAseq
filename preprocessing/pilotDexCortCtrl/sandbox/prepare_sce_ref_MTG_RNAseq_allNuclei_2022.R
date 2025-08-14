sour

sce_ref_MTGRNAseqAllNuclei2022 <- readH5AD(
  "data/celltype_references/human_allenBrainAtlas/MTG_10xSEA-AD_2022/Reference_MTG_RNAseq_all-nuclei.2022-06-07.h5ad"
)

# 1. Kopiuj X do counts
assay(sce_ref_MTGRNAseqAllNuclei2022, "counts") <- assay(sce_ref_MTGRNAseqAllNuclei2022, "X")

# 2. Usuń X
assays(sce_ref_MTGRNAseqAllNuclei2022) <- assays(sce_ref_MTGRNAseqAllNuclei2022)[setdiff(assayNames(sce_ref_MTGRNAseqAllNuclei2022), "X")]


# --- 1. Kopia oryginalnego obiektu ---

orig <- sce_ref_MTGRNAseqAllNuclei2022
ids  <- sub("\\.\\d+$", "", rownames(orig))  # bez sufiksów

# agregacja counts po symbolu; zwraca SCE i zachowuje colData
sce_ref_clean <- aggregateAcrossFeatures(
  orig,
  ids            = ids,
  use.assay.type = "counts",  # w nowszych scuttle jest 'use.assay.type'; w starszych 'assay.type'
  average        = FALSE      # sumuj counts
)

# log-normalizacja po agregacji
sce_ref_clean <- scuttle::logNormCounts(sce_ref_clean)

assayNames(sce_ref_clean)

assay(sce_ref_clean, "counts") <- as(assay(sce_ref_clean, "counts"), "dgCMatrix")
assay(sce_ref_clean, "logcounts") <- as(assay(sce_ref_clean, "logcounts"), "dgCMatrix")

saveHDF5SummarizedExperiment(
  sce_ref_clean,
  dir = "data/reference_celltype_expression_processed/sce_ref_clean_MTGRNAseqAllNuclei2022.hdf5",
  replace = TRUE
)


sce_ref_clean_MTGRNAseqAllNuclei2022_hdf5



sce_ref_MTGRNAseqAllNuclei2022


dim(assay(sce_ref_MTGRNAseqAllNuclei2022, "X")) 

colnames(colData(sce_ref_MTGRNAseqAllNuclei2022))    # nazwy kolumn z metadanych
head(colData(sce_ref_MTGRNAseqAllNuclei2022))     


table(sce_ref_MTGRNAseqAllNuclei2022$CA_subclass_label)

summary(as.vector(assay(sce_ref_MTGRNAseqAllNuclei2022, "X")[, 1:100]))


# 1. Kopiuj X do counts
assay(sce_ref_MTGRNAseqAllNuclei2022, "counts") <- assay(sce_ref_MTGRNAseqAllNuclei2022, "X")

# 2. Usuń X
assays(sce_ref_MTGRNAseqAllNuclei2022) <- assays(sce_ref_MTGRNAseqAllNuclei2022)[setdiff(assayNames(sce_ref_MTGRNAseqAllNuclei2022), "X")]

# 3. Zrób logcounts
sce_ref_MTGRNAseqAllNuclei2022 <- logNormCounts(sce_ref_MTGRNAseqAllNuclei2022)

# 4. Sprawdź
assayNames(sce_ref_MTGRNAseqAllNuclei2022)

head(rowData(sce_ref_MTGRNAseqAllNuclei2022))

sce_ref_MTGRNAseqAllNuclei2022 %>% str

rownames(sce_ref_MTGRNAseqAllNuclei2022) %>%
  as.data.frame() %>% dim



# Skopiuj, żeby nie modyfikować oryginału
sce_ref_clean <- sce_ref_MTGRNAseqAllNuclei2022

# Usuń suffix .digit
gene_names_clean <- sub("\\.\\d+$", "", rownames(sce_ref_clean))

# Agregacja duplikatów przez sumowanie
library(Matrix)
library(SingleCellExperiment)

if (anyDuplicated(gene_names_clean)) {
  counts_mat <- assay(sce_ref_clean, "counts")
  counts_agg <- as.matrix(rowsum(as.matrix(counts_mat), group = gene_names_clean))
  sce_ref_clean <- SingleCellExperiment(
    assays = list(counts = counts_agg)
  )
  rownames(sce_ref_clean) <- unique(gene_names_clean)
} else {
  rownames(sce_ref_clean) <- gene_names_clean
}

rm(counts_agg, counts_mat)
rm(sce_ref_MTGRNAseqAllNuclei2022)


assayNames(sce_ref_clean)


rownames(sce_ref_clean) %>%
  as.data.frame()


table(sce_ref_clean$CA_subclass_label)


predict_test_pilotDexCortCtrl <- predict_celltypes_singler_multi_v1(
  sce_query = sce_integrated_pilotDexCortCtrl_multi,
  sce_ref   = sce_ref_clean,                   # z prepare_expression_reference_singler_v1
  label_col = "CA_subclass_label",
  cluster_cols = c( "seurat_clusters_res0.04", "seurat_clusters_res0.13", "seurat_clusters_res0.2"),
  min_common_genes = 1000,
  save_dir = "data/processed_reference_celltype_expression/test",
  verbose = TRUE
)

