# ##############################################################################
# ---- import packages, functions and data ----
# ##############################################################################

source("preprocessing/src/setup_packages.R")

# ##############################################################################
# ---- run analysis ----
# ##############################################################################

# Start total timer
script_start_time <- Sys.time()

# Parameters
cluster_cols <- c("seurat_clusters_res0.04", "seurat_clusters_res0.13", "seurat_clusters_res0.2")
min_pct_vals <- c(0.9, 0.8, 0.75, 0.5)

p_adj_thresh       <- 0.01
logfc_thresh       <- 1
top_n_featureplots <- 9
enrichr_direction  <- "up"
enrichr_databases  <- c(
  "CellMarker_2024",
  "CellMarker_Augmented_2021",
  "Human_Gene_Atlas",
  "GO_Biological_Process_2025",
  "GO_Cellular_Component_2025",
  "Jensen_DISEASES_Experimental_2025"
)

# Output folder
out_dir <- "data/pilotDexCortCtrl/r_outputs/cluster_markers"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Input Seurat object
seurat_obj <- integrated_pilotDexCortCtrl

# Total number of iterations
total_iterations <- length(cluster_cols) * length(min_pct_vals)
iteration <- 0

# Loop
for (cl_col in cluster_cols) {
  for (mp in min_pct_vals) {
    iteration <- iteration + 1
    
    # Start iteration timer
    iter_start_time <- Sys.time()
    
    # Build filenames
    base_filename <- glue(
      "clusterMarkersAndPlotsTopAndEnrichr_{cl_col}_fdr{p_adj_thresh}logFC{logfc_thresh}minpct{mp}"
    )
    xlsx_path  <- file.path(out_dir, paste0(base_filename, ".xlsx"))
    rdata_path <- file.path(out_dir, paste0(base_filename, ".RData"))
    
    # Log start
    cat("\n", strrep("─", 60), "\n", sep = "")
    cat(glue("[Iteration {iteration}/{total_iterations}] START\n"))
    cat(glue("Cluster column  : {cl_col}\n"))
    cat(glue("min.pct         : {mp}\n"))
    cat(strrep("─", 60), "\n", sep = "")
    
    # Run analysis
    cluster_markers_enrichr_results <- analyze_markers_enrichr_to_xlsx(
      seurat_obj         = seurat_obj,
      xlsx_path          = xlsx_path,
      cluster_col        = cl_col,
      assay              = "integrated",
      min.pct            = mp,
      p_adj_thresh       = p_adj_thresh,
      logfc_thresh       = logfc_thresh,
      top_n_featureplots = top_n_featureplots,
      enrichr_direction  = enrichr_direction,
      enrichr_databases  = enrichr_databases
    )
    
    # Save RData
    save(cluster_markers_enrichr_results, file = rdata_path)
    
    # Log end
    iter_end_time <- Sys.time()
    iter_time <- round(as.numeric(difftime(iter_end_time, iter_start_time, units = "secs")), 2)
    cat(strrep("─", 60), "\n", sep = "")
    cat(glue("[Iteration {iteration}/{total_iterations}] DONE in {iter_time} sec\n"))
    cat(glue("XLSX  → {xlsx_path}\n"))
    cat(glue("RData → {rdata_path}\n"))
    cat(strrep("─", 60), "\n", sep = "")
  }
}

# End total timer
script_end_time <- Sys.time()
total_time <- round(as.numeric(difftime(script_end_time, script_start_time, units = "mins")), 2)

cat("\n", strrep("═", 60), "\n", sep = "")
cat(glue("All {total_iterations} iterations completed in {total_time} min\n"))
cat(strrep("═", 60), "\n", sep = "")

# outs <- analyze_markers_enrichr_to_xlsx(
#   seurat_obj         = integrated_pilotDexCortCtrl,
#   xlsx_path          = "data/pilotDexCortCtrl/r_outputs/cluster_markers/clusterMarkersAndPlotsExampleGenesEnrichr_integratedPilotDexCortCtrl_res0.04_fdr0.01logFC1minpct0.95.xlsx",
#   cluster_col        = "seurat_clusters_res0.04",
#   assay              = "integrated",
#   min.pct            = 0.9,
#   p_adj_thresh       = 0.01,
#   logfc_thresh       = 1,            # symetryczny próg FC
#   top_n_featureplots = 9,
#   enrichr_direction  = "up",         # "up", "down" lub "both"
#   enrichr_databases  = c(
#     "CellMarker_2024",
#     "CellMarker_Augmented_2021",
#     "Human_Gene_Atlas",
#     "GO_Biological_Process_2025",
#     "GO_Cellular_Component_2025",
#     "Jensen_DISEASES_Experimental_2025"
#   )
# )
# 

