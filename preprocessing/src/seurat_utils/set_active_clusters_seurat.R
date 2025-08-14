set_active_clusters_seurat <- function(seurat_obj, cluster_var = "seurat_clusters") {
  if (!cluster_var %in% colnames(seurat_obj@meta.data)) {
    stop("Cluster column not found in metadata")
  }
  Idents(seurat_obj) <- seurat_obj[[cluster_var]]
  return(seurat_obj)
}
