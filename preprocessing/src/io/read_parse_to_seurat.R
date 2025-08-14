read_parse_to_seurat <- function(
    dge_dir,
    project = "ParseProject",
    use_gene_symbol = TRUE
) {
  # Wczytaj macierz zliczeń (i przetransponuj, bo geny są w kolumnach)
  count_matrix <- readMM(file.path(dge_dir, "count_matrix.mtx")) %>% t()
  
  # Wczytaj geny
  genes <- fread(file.path(dge_dir, "all_genes.csv"))
  
  if (use_gene_symbol) {
    gene_ids <- genes$gene_name
  } else {
    gene_ids <- paste0(genes$gene_id, "_", genes$gene_name)
  }
  
  rownames(count_matrix) <- gene_ids
  
  # Wczytaj metadane komórek
  cell_metadata <- fread(file.path(dge_dir, "cell_metadata.csv"))
  cell_ids <- cell_metadata$bc_wells
  colnames(count_matrix) <- cell_ids
  
  # Utwórz obiekt Seurat
  seurat_obj <- CreateSeuratObject(
    counts = count_matrix,
    meta.data = as.data.frame(cell_metadata),
    project = project
  )
  
  return(seurat_obj)
}
