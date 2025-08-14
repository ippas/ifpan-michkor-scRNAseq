
print("# ##############################################################################")
print("# ---- 1. import packages and functions ----")
print("# ##############################################################################")

source("preprocessing/src/setup_packages.R")

# ##############################################################################
# ---- read data ----
# ##############################################################################


print("# ##############################################################################")
print("# ---- 2. preprocessing data ----")
print("# ##############################################################################")

dge_dir <- "data/pilotDexCortCtrl/combined_BG60BLUE/all-well/DGE_filtered"
combined_BG60BLUE <- read_parse_to_seurat(dge_dir = dge_dir, project = "ParsePilot", use_gene_symbol = FALSE)

dge_dir <- "data/pilotDexCortCtrl/combined_BG60GREEN/all-well/DGE_filtered"
combined_BG60GREEN <- read_parse_to_seurat(dge_dir = dge_dir, project = "ParsePilot", use_gene_symbol = FALSE)

sample_metadata <- data.frame(
  sample_id = sprintf("%02d", 1:12),  # czyli "01", "02", ..., "12"
  well_label = paste0("A", 1:12),
  condition = c(
    rep("CTRL", 4),   # A1–A4
    rep("DEX", 4),    # A5–A8
    rep("CORT", 4)    # A9–A12
  )
)
# ##############################################################################
# ---- merged sublibraries ----
# ##############################################################################
prefixes <- c("BG60BLUE", "BG60GREEN")

integrated_pilotDexCortCtrl <- merge_multiple_seurat_objects(
  seurat_list = list(combined_BG60BLUE, combined_BG60GREEN),
  project_name = "mergedParsePilot",
  prefixes = prefixes, 
  metadata_column_name = "sublibrary"
)


integrated_pilotDexCortCtrl$treatment <- case_when(
  integrated_pilotDexCortCtrl$bc1_well %in% c("A1", "A2", "A3", "A4")  ~ "CTRL",
  integrated_pilotDexCortCtrl$bc1_well %in% c("A5", "A6", "A7", "A8")  ~ "DEX",
  integrated_pilotDexCortCtrl$bc1_well %in% c("A9", "A10", "A11", "A12") ~ "CORT",
  TRUE ~ NA_character_
)


integrated_pilotDexCortCtrl$sample_id <- with(integrated_pilotDexCortCtrl@meta.data,
                                          paste0(trimws(sublibrary), "_", trimws(bc1_well)))

# ##############################################################################
# ---- run seruat analysis ----
# ##############################################################################
# Preprocessing step

print("# ##############################################################################")
print("# ---- 3. Integrate data ----")
print("# ##############################################################################")

integrated_pilotDexCortCtrl <- seuratV5_precluster_preprocessing_integrate(
  seurat_obj      = integrated_pilotDexCortCtrl,
  split_by        = "sample_id",          # albo "sublibrary" itd.
  method          = "LogNormalize",       # lub "SCT"
  nfeatures       = 2000,
  pca_npcs        = 50,
  vars_to_regress = NULL,                 # np. c("nCount_RNA")
  parallel        = TRUE,
  n_cores         = 10,
  future_plan     = "multisession",
  maxsize_gb      = 100,
  blas_threads    = NULL,
  verbose         = TRUE,
  seed            = 777
)


print("# ##############################################################################")
print("# ---- 4. Save backup integrate data ----")
print("# ##############################################################################")

save(
  integrated_pilotDexCortCtrl,
  file = "data/pilotDexCortCtrl/r_outputs/integrated_pilotDexCortCtrl_backupFirstIntegration.RData"
)


print("# ##############################################################################")
print("# ---- 5. Clustering data ----")
print("# ##############################################################################")

# Clustering + UMAP step with multiple resolutions
integrated_pilotDexCortCtrl <- seurat_multiRes_clustering_umap(
  seurat_obj = integrated_pilotDexCortCtrl,
  dims = 1:20,
  resolution = c(0.01, 0.02, 0.021, 0.022, 0.025, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.13, 0.15, 0.17, 0.2, 0.3, 0.5, 0.4, 0.8, 1),
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


print("# ##############################################################################")
print("# ---- 6. Save integrate and clustering data----")
print("# ##############################################################################")

save(
  integrated_pilotDexCortCtrl,
  file = "data/pilotDexCortCtrl/r_outputs/integrated_pilotDexCortCtrl.RData"
)


print("# ##############################################################################")
print("# ---- 7. End analysis----")
print("# ##############################################################################")


