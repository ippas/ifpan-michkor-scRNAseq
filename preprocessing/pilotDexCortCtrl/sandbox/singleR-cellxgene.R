BiocManager::install("zellkonverter")


library(zellkonverter)


sce_ref <- readH5AD("data/celltype_references/cellxgene/86afbf8c-e49c-4e71-9ba3-e07b571f1acf.h5ad")

assay(sce_ref, "logcounts") <- assay(sce_ref, "X")

dict <- build_human_mouse_gene_dictionary(
  ensembl_version = 110
)


predict_test_pilotDexCortCtrl <- predict_celltypes_singler_multi_v1(
  sce_query = sce_integrated_pilotDexCortCtrl_multi,
  sce_ref   = sce_ref,                   # z prepare_expression_reference_singler_v1
  label_col = "sub_celltype",
  cluster_cols = c( "seurat_clusters_res0.04", "seurat_clusters_res0.13", "seurat_clusters_res0.2"),
  min_common_genes = 1000,
  save_dir = "data/processed_reference_celltype_expression/test",
  verbose = TRUE
)

sce_ref2 <- sce_ref

colnames(colData(sce_ref))

head(colData(sce_ref))



assayNames(sce_ref)


sce_ref$sub_celltype %>% table

summary(as.vector(assay(sce_ref, "X")[,1:100]))


# Podmiana Ensembl -> symbole
gene_symbols <- toupper(as.character(rowData(sce_ref)$gene_symbols))

# Upewnij się, że nie są puste/NA
ok <- !is.na(gene_symbols) & gene_symbols != ""
sce_ref <- sce_ref[ok, ]
rownames(sce_ref) <- gene_symbols[ok]

# Usuń duplikaty symboli (jeśli są)
sce_ref <- sce_ref[!duplicated(rownames(sce_ref)), ]

# Sprawdź, czy są wspólne geny z query
length(intersect((rownames(sce_ref)), rownames(sce_integrated_pilotDexCortCtrl_multi)))
