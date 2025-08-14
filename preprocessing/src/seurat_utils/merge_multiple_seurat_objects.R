merge_multiple_seurat_objects <- function(
    seurat_list,
    project_name = "merged_project",
    prefixes = NULL,
    metadata_column_name = "sublibrary"
) {
  stopifnot(length(seurat_list) >= 2)
  
  if (!is.null(prefixes)) {
    stopifnot(length(prefixes) == length(seurat_list))
  } else {
    prefixes <- paste0("sub", seq_along(seurat_list))
  }
  
  # Dodaj prefixy do nazw komórek + dodaj wskazaną kolumnę do meta.data
  seurat_list <- Map(function(obj, prefix) {
    obj <- RenameCells(obj, add.cell.id = prefix)
    obj$orig.ident <- prefix  # można też zostawić oryginalne jeśli wolisz
    obj[[metadata_column_name]] <- prefix  # <--- dynamiczna kolumna
    return(obj)
  }, seurat_list, prefixes)
  
  # Merge wszystkie obiekty
  merged <- Reduce(function(x, y) merge(x, y = y, project = project_name), seurat_list)
  
  return(merged)
}
