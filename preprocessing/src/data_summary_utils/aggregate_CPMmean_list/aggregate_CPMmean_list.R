aggregate_CPMmean_list <- function(
    seurat_obj,
    cluster_col   = "seurat_clusters",
    sample_col    = "sample_id",
    assay         = "RNA",
    layer         = "counts",
    genes_keep     = NULL,
    select_samples = NULL,
    scale_factor   = 1e6,
    progress       = TRUE
) {
  
  # ---- 1: aggregate SUM counts per cluster × sample ----
  sum_list <- aggregate_counts_by_cluster_sample(
    seurat_obj     = seurat_obj,
    cluster_col    = cluster_col,
    sample_col     = sample_col,
    assay          = assay,
    layer          = layer,
    agg            = "sum",   # <--- SUMA
    return         = "list",
    genes_keep     = genes_keep,
    select_samples = select_samples,
    progress       = progress
  )
  
  # ---- 2: compute CPM from SUM counts ----
  cpm_sum_list <- normalize_CPM_perCluster(
    count_list   = sum_list,
    scale_factor = scale_factor,
    progress     = progress
  )
  
  # ---- 3: mean CPM per cell (divide CPM by number of cells) ----
  meta <- seurat_obj@meta.data
  clusters <- meta[[cluster_col]]
  samples  <- meta[[sample_col]]
  grp <- interaction(clusters, samples, drop = TRUE)
  cell_counts <- table(grp)  # number of cells in each cluster x sample
  
  cpm_mean_list <- lapply(names(cpm_sum_list), function(cl_nm) {
    m <- cpm_sum_list[[cl_nm]]   # CPM from sum
    samp_names <- colnames(m)
    cl_id <- sub("^cluster_", "", cl_nm)
    
    # number of cells per sample in this cluster
    n_cells <- as.numeric(cell_counts[paste0(cl_id, ".", samp_names)])
    
    # divide CPM by cell count => mean CPM per cell
    sweep(m, 2, n_cells, "/")
  })
  
  names(cpm_mean_list) <- names(cpm_sum_list)
  return(cpm_mean_list)
}
