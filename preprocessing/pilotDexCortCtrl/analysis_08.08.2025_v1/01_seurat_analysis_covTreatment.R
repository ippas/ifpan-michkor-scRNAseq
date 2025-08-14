
# ##############################################################################
# ---- read data ----
# ##############################################################################
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

merged_pilotDexCortCtrl_covTreatment <- merge_multiple_seurat_objects(
  seurat_list = list(combined_BG60BLUE, combined_BG60GREEN),
  project_name = "mergedParsePilot",
  prefixes = prefixes, 
  metadata_column_name = "sublibrary"
)


merged_pilotDexCortCtrl_covTreatment$treatment <- case_when(
  merged_pilotDexCortCtrl_covTreatment$bc1_well %in% c("A1", "A2", "A3", "A4")  ~ "CTRL",
  merged_pilotDexCortCtrl_covTreatment$bc1_well %in% c("A5", "A6", "A7", "A8")  ~ "DEX",
  merged_pilotDexCortCtrl_covTreatment$bc1_well %in% c("A9", "A10", "A11", "A12") ~ "CORT",
  TRUE ~ NA_character_
)


# ##############################################################################
# ---- run seruat analysis ----
# ##############################################################################
# Preprocessing step
merged_pilotDexCortCtrl_covTreatment <- seurat_precluster_preprocessing(
  seurat_obj = merged_pilotDexCortCtrl_covTreatment,
  n_variable_features = 2000,
  selection_method = "vst",
  vars_to_regress = "treatment",
  verbose = TRUE
)


merged_pilotDexCortCtrl_covTreatment@assays$RNA@

rm(merged_pilotDexCortCtrl)

# Clustering + UMAP step with multiple resolutions
merged_pilotDexCortCtrl_covTreatment <- seurat_multiRes_clustering_umap(
  seurat_obj = merged_pilotDexCortCtrl_covTreatment,
  dims = 1:20,
  resolution = c(0.26),
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


save(
  merged_pilotDexCortCtrl_covTreatment,
  file = "data/pilotDexCortCtrl/r_outputs/merged_pilotDexCortCtrl_covTreatment.RData"
)



