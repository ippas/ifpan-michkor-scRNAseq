seurat_merge_all_counts_layers <- function(
    seurat_obj,
    assay = "RNA",
    remove_old = FALSE
) {
  # Check if assay exists
  if (!assay %in% names(seurat_obj@assays)) {
    stop(paste0("Assay '", assay, "' not found in the Seurat object."))
  }
  
  assay_obj <- seurat_obj[[assay]]
  
  # Identify all count layers (matching pattern 'counts.*')
  count_layers <- grep("^counts(\\.|$)", Layers(assay_obj), value = TRUE)
  if (length(count_layers) == 0) {
    stop("No count layers found.")
  }
  
  message("Merging ", length(count_layers), " layers: ", paste(count_layers, collapse = ", "))
  
  # Extract each count layer as a matrix and merge
  merged_counts <- do.call(
    cbind,
    lapply(count_layers, function(layer_name) {
      mat <- LayerData(assay_obj, layer = layer_name)
      # Ensure proper column naming
      if (is.null(colnames(mat))) {
        stop(paste("Layer", layer_name, "has no column names."))
      }
      mat
    })
  )
  
  # Reorder columns to match the Seurat object
  merged_counts <- merged_counts[, colnames(seurat_obj)]
  
  # Assign merged counts to a new "counts" layer
  LayerData(seurat_obj[[assay]], layer = "counts") <- merged_counts
  
  # Optionally remove the old count layers
  if (remove_old) {
    for (layer_name in count_layers) {
      assay_obj@layers[[layer_name]] <- NULL
    }
    seurat_obj[[assay]] <- assay_obj
  }
  
  return(seurat_obj)
}
