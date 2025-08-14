
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


# ##############################################################################
# ---- multi resolution predict ----
# ##############################################################################

# 1) Query z wieloma kolumnami klastrów
sce_integrated_pilotDexCortCtrl_multi <- prepare_expression_query_multi_singler_v1(
  seurat_obj = integrated_pilotDexCortCtrl,
  assay = "RNA",
  data_layer_regex = "^data",
  normalize_if_missing = TRUE,
  gene_name_from_rownames = TRUE,
  gene_sep = "-",
  agg_fun_for_duplicates = "mean",
  cluster_cols = c( "seurat_clusters_res0.04", "seurat_clusters_res0.13", "seurat_clusters_res0.2"),
  alias_first_as = "cluster",     # SingleR-friendly alias
  save_rds = "data/reference_celltype_expression_processed/sce_integrated_pilotDexCortCtrl_res_0.04_0.13_0.2.rds",
  verbose = TRUE
)



# 2) Predykcja dla każdej kolumny klastrów
predict_GSE174367_pilotDexCortCtrl <- predict_celltypes_singler_multi_v1(
  sce_query = sce_integrated_pilotDexCortCtrl_multi,
  sce_ref   = sce_ref,                   # z prepare_expression_reference_singler_v1
  label_col = "Cell.Type.Full",
  cluster_cols = c( "seurat_clusters_res0.04", "seurat_clusters_res0.13", "seurat_clusters_res0.2"),
  min_common_genes = 1000,
  save_dir = "data/processed_reference_celltype_expression/singler_GSE174367_pilotDexCortCtrl",
  verbose = TRUE
)


predict_GSE174367_pilotDexCortCtrl$results_per_clustering$seurat_clusters_res0.04$annotation_table


# ##############################################################################
install.packages(
  "cellxgene.census",
  repos = c("https://chanzuckerberg.r-universe.dev", "https://cloud.r-project.org")
)

library(cellxgene.census)

# Pobierz obiekt Seurat z datasetu
# (z dokumentacji: get_seurat lub get_sce)
sce_ref <- get_sce(
  collection_id = "6d7d23d0-237d-4430-9200-92858abba2d8",
  dataset_id    = "86afb8fc-e49c-4e71-9ba3-e07b57f1facf"
)

# Sprawdź dostępne kolumny w metadanych
names(colData(sce_ref))

