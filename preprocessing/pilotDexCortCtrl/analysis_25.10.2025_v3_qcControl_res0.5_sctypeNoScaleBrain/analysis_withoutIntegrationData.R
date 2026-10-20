# ======================================================================
# UMAP bez integracji (na dokładnie tych samych komórkach co w obiekcie filtrowanym)
# ======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

# ---- 1. Wczytaj surowe dane (niezintegrowane) ----
dge_dir_blue  <- "data/pilotDexCortCtrl/combined_BG60BLUE/all-well/DGE_filtered"
dge_dir_green <- "data/pilotDexCortCtrl/combined_BG60GREEN/all-well/DGE_filtered"

combined_BG60BLUE  <- read_parse_to_seurat(dge_dir = dge_dir_blue,  project = "ParsePilot", use_gene_symbol = FALSE)
combined_BG60GREEN <- read_parse_to_seurat(dge_dir = dge_dir_green, project = "ParsePilot", use_gene_symbol = FALSE)

noIntegrated_raw <- merge_multiple_seurat_objects(
  seurat_list = list(combined_BG60BLUE, combined_BG60GREEN),
  project_name = "mergedParsePilot_noInt",
  prefixes = c("BG60BLUE", "BG60GREEN"),
  metadata_column_name = "sublibrary"
)

# ---- 2. Dopasuj meta dane i filtry ----
# ważne: pobieramy identyczne komórki jak w obiekcie zintegrowanym
filtered_cells <- rownames(integrated_pilotDexCortCtrl_filtered@meta.data)
noIntegrated_raw <- subset(noIntegrated_raw, cells = filtered_cells)

# ---- 3. Przypisz treatment i sample_id tak jak wcześniej ----
noIntegrated_raw$treatment <- case_when(
  noIntegrated_raw$bc1_well %in% c("A1", "A2", "A3", "A4")  ~ "CTRL",
  noIntegrated_raw$bc1_well %in% c("A5", "A6", "A7", "A8")  ~ "DEX",
  noIntegrated_raw$bc1_well %in% c("A9", "A10", "A11", "A12") ~ "CORT",
  TRUE ~ NA_character_
)

noIntegrated_raw$sample_id <- with(noIntegrated_raw@meta.data,
                                   paste0(trimws(sublibrary), "_", trimws(bc1_well))
)

# ---- 4. Klasyczne przetwarzanie bez integracji ----
DefaultAssay(noIntegrated_raw) <- "RNA"

noIntegrated_raw <- NormalizeData(noIntegrated_raw, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
noIntegrated_raw <- FindVariableFeatures(noIntegrated_raw, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
noIntegrated_raw <- ScaleData(noIntegrated_raw, features = VariableFeatures(noIntegrated_raw), verbose = FALSE)
noIntegrated_raw <- RunPCA(noIntegrated_raw, features = VariableFeatures(noIntegrated_raw), npcs = 50, reduction.name = "pca_noInt", seed.use = 777)
noIntegrated_raw <- RunUMAP(noIntegrated_raw, reduction = "pca_noInt", dims = 1:20, reduction.name = "umap_noInt", seed.use = 777)

# ---- 5. Zapisz wynik ----
save(noIntegrated_raw, file = "data/pilotDexCortCtrl/r_outputs/noIntegrated_pilotDexCortCtrl-filteredCells-27.10.2025.RData")

# ---- 6. Podgląd (opcjonalnie) ----
DimPlot(noIntegrated_raw, reduction = "umap_noInt", group.by = "treatment", pt.size = 0.3)

noIntegrated_raw@reductions$umap <- noIntegrated_raw@reductions$umap_noInt
umap <- Embeddings(noIntegrated_raw, reduction = "umap_noInt")
colnames(umap) <- c("umap_1", "umap_2")
noIntegrated_raw@reductions$umap <- CreateDimReducObject(
  embeddings = umap,
  key = "umap_",
  assay = DefaultAssay(noIntegrated_raw)
)


plot_umap_cluster_faceted(
  noIntegrated_raw,
  cluster_var = "sample",
  plot_title = "",
  facet_var = "sample",
  label_clusters = TRUE,
  text_size_axis_title = 18,
  text_size_axis_text = 18,
  text_size_legend_title = 18,
  text_size_legend_text = 18,
  text_size_facet = 18,
  text_size_plot_title = 18,
  text_size_label = 7,
  n_row_legend = 4
) -> p1


plot_umap_cluster_faceted(
  noIntegrated_raw,
  cluster_var = "treatment",
  plot_title = "",
  facet_var = "sample",
  label_clusters = TRUE,
  text_size_axis_title = 18,
  text_size_axis_text = 18,
  text_size_legend_title = 18,
  text_size_legend_text = 18,
  text_size_facet = 18,
  text_size_plot_title = 18,
  text_size_label = 7,
  n_row_legend = 4
) -> p1


plot_umap_cluster_faceted(
  integrated_pilotDexCortCtrl_filtered,
  cluster_var = "treatment",
  plot_title = "",
  facet_var = "sample",
  label_clusters = TRUE,
  text_size_axis_title = 18,
  text_size_axis_text = 18,
  text_size_legend_title = 18,
  text_size_legend_text = 18,
  text_size_facet = 18,
  text_size_plot_title = 18,
  text_size_label = 7,
  n_row_legend = 4
) -> p2

p1 + p2


# Ścieżka do folderu docelowego
out_dir <- "/home/mateusz/projects/ifpan-michkor-scRNAseq/results/results_analysis_25.10.2025_v3_qcControl_res0.5_sctypeNoScaleBrain/presentation-29.10.2025"

# Upewnij się, że folder istnieje
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Pełna ścieżka do pliku wynikowego
out_file <- file.path(out_dir, "p1_p2_combined.svg")

# Zapis wykresu do pliku SVG
svg(filename = out_file, width = 12, height = 8)  # szerokość ~2× wysokość
print(p1 + p2)
dev.off()

message("✅ Zapisano plik: ", out_file)

