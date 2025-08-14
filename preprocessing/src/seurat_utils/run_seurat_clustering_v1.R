# src/seurat_utils/run_seurat_clustering.R

run_seurat_clustering <- function(seurat_obj,
                                  dims_pca = 1:20,
                                  resolution = 0.5,
                                  run_umap = TRUE,
                                  dims_umap = 1:10) {
  seurat_obj <- FindNeighbors(seurat_obj, dims = dims_pca)
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
  
  if (run_umap) {
    seurat_obj <- RunUMAP(seurat_obj, dims = dims_umap)
  }
  
  return(seurat_obj)
}