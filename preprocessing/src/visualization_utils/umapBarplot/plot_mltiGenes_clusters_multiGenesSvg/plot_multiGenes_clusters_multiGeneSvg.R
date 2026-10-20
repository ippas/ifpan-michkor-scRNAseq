plot_multiGenes_clusters_multiGenesSvg <- function(
    summary_list,
    cluster_names,
    genes_vector,
    output_dir,
    width  = 16,
    height = 20,
    force_create_dir = TRUE,
    gene_title_size = 14,
    ...
) {
  require(purrr)
  require(ggplot2)
  
  # ---- ensure folder exists ----
  if (!dir.exists(output_dir)) {
    if (force_create_dir) {
      dir.create(output_dir, recursive = TRUE)
      message("📁 Directory created: ", output_dir)
    } else {
      stop("❌ Output directory does not exist. Set force_create_dir = TRUE to create it automatically.")
    }
  }
  
  total <- length(genes_vector)
  message("🚀 Starting processing of ", total, " genes...")
  
  # ---- loop over genes ----
  for (i in seq_along(genes_vector)) {
    gene <- genes_vector[i]
    message("▶ ", i, "/", total, "  ", gene)
    
    tryCatch({
      filename <- paste0(gene, ".svg")   # <--- no suffix
      filepath <- file.path(output_dir, filename)
      
      svg(filepath, width = width, height = height)
      
      p <- plot_gene_multiClusters(
        summary_list  = summary_list,
        cluster_names = cluster_names,
        gene_symbol   = gene,
        gene_title_size = gene_title_size,
        ...
      )
      
      print(p)
      dev.off()
      
    }, error = function(e) {
      message("⚠ Error for gene ", gene, ": ", e$message)
      dev.off()
    })
  }
  
  message("🎉 DONE: Output saved to → ", output_dir)
}
