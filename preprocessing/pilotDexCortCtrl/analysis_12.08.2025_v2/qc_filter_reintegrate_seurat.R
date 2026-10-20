# ==============================================================================
# QC FILTER + REINTEGRATION + CLUSTERING SCRIPT
# ==============================================================================

message("[STEP 0] Ładowanie pakietów i ustawień...")
source("preprocessing/src/setup_packages.R")

# ------------------------------------------------------------------------------
# Prepare barcode list for filtering
# ------------------------------------------------------------------------------
message("[STEP 1] Przygotowanie listy genów mitochondrialnych...")
mitochondrial_name_id <- integrated_pilotDexCortCtrl[["RNA"]] %>%
  rownames() %>%
  as.data.frame() %>%
  set_colnames("name_id") %>%
  filter(grepl("-MT-", name_id)) %>%
  .$name_id

message("[STEP 2] Pobieranie macierzy zliczeń RNA...")
rna_counts <- GetAssayData(integrated_pilotDexCortCtrl, assay = "RNA", slot = "counts")

message("[STEP 3] Identyfikacja genów mitochondrialnych w macierzy...")
mt_genes <- intersect(rownames(rna_counts), mitochondrial_name_id)

message("[STEP 4] Obliczanie % mitochondrialnego RNA dla każdej komórki...")
percent_mt <- (colSums(rna_counts[mt_genes, , drop = FALSE]) /
                 colSums(rna_counts)) * 100

message("[STEP 5] Tworzenie tabeli QC i przypisywanie wyników...")
qc_summary_df <- integrated_pilotDexCortCtrl@meta.data %>%
  rownames_to_column("cell_id") %>%
  mutate(percent_MTDNA = percent_mt) %>%
  mutate(
    fail_reason_500 = paste(
      ifelse(percent_MTDNA > 10, "high_mt", ""),
      ifelse(nCount_RNA <= 500, "low_counts", ""),
      ifelse(nCount_RNA >= 10000, "high_counts", ""),
      ifelse(nFeature_RNA <= 200, "low_genes", ""),
      ifelse(nFeature_RNA >= 4000, "high_genes", ""),
      sep = ";"
    ),
    qc_pass_500 = fail_reason_500 == ";;;;",
    fail_reason_1000 = paste(
      ifelse(percent_MTDNA > 10, "high_mt", ""),
      ifelse(nCount_RNA <= 1000, "low_counts", ""),
      ifelse(nCount_RNA >= 20000, "high_counts", ""),
      ifelse(nFeature_RNA <= 200, "low_genes", ""),
      ifelse(nFeature_RNA >= 6000, "high_genes", ""),
      sep = ";"
    ),
    qc_pass_1000 = fail_reason_1000 == ";;;;"
  )

message("[STEP 6] Filtrowanie komórek na podstawie QC...")
cells_keep <- qc_summary_df %>%
  filter(qc_pass_500) %>%
  pull(cell_id)

message(paste0("[INFO] Zachowano ", length(cells_keep), " komórek po filtracji QC."))

# ------------------------------------------------------------------------------
# Subset + Reintegrate
# ------------------------------------------------------------------------------
message("[STEP 7] Tworzenie nowego obiektu Seurat po filtracji...")
integrated_pilotDexCortCtrl_filtered <- subset(integrated_pilotDexCortCtrl, cells = cells_keep)

if ("integrated" %in% names(Assays(integrated_pilotDexCortCtrl_filtered))) {
  message("[INFO] Usuwam stary assay 'integrated' przed ponowną integracją.")
  integrated_pilotDexCortCtrl_filtered[["integrated"]] <- NULL
}

message("[STEP 8] Ustawienie domyślnego assay na 'RNA'...")
DefaultAssay(integrated_pilotDexCortCtrl_filtered) <- "RNA"

message("[STEP 9] Rozpoczęcie ponownej integracji danych...")
integrated_pilotDexCortCtrl_filtered <- seuratV5_precluster_preprocessing_integrate(
  seurat_obj      = integrated_pilotDexCortCtrl_filtered,
  split_by        = "sample_id",
  method          = "LogNormalize",
  nfeatures       = 2000,
  pca_npcs        = 50,
  parallel        = FALSE,
  n_cores         = 10,
  future_plan     = "multisession",
  seed            = 777
)

message("[STEP 10] Zapisywanie zintegrowanego obiektu do pliku RData (v1)...")
save(
  integrated_pilotDexCortCtrl_filtered,
  file = "data/pilotDexCortCtrl/r_outputs/integrated_pilotDexCortCtrl_filtered_14.09.2025_v1.RData"
)

# ------------------------------------------------------------------------------
# Clustering + UMAP
# ------------------------------------------------------------------------------
message("[STEP 11] Rozpoczęcie klasteryzacji i UMAP...")
integrated_pilotDexCortCtrl_filtered <- seurat_multiRes_clustering_umap(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  dims = 1:20,
  resolution = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1),
  k.param = 20,
  nn.method = "rann",
  algorithm = 1,
  reduction = "pca",
  do_umap = TRUE,
  umap_n_neighbors = 30,
  umap_min_dist = 0.3,
  umap_spread = 1.0,
  umap_metric = "cosine",
  verbose = TRUE,
  seed = 777
)

message("[STEP 12] Zapisywanie obiektu po klasteryzacji i UMAP (v2)...")
save(
  integrated_pilotDexCortCtrl_filtered,
  file = "data/pilotDexCortCtrl/r_outputs/integrated_pilotDexCortCtrl_filtered_14.09.2025_v2.RData"
)

message("[DONE] Skrypt zakończony pomyślnie ✅")
