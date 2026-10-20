seurat_featurePlot_binary_clusterSummary_generateMultiSVG <- function(
    seurat_obj,
    genes_vector,
    output_dir,
    split.by = "treatment",
    split_levels = c("CTRL", "CORT", "DEX"),
    cluster_col = "res0.3_RNA_noScaled_brain",
    cluster_summary_figure = c(
      "Astrocytes" = "Astrocytes",
      "Oligodendrocytes" = "Oligodendrocytes",
      "Oligodendrocyte precursor cells" = "OPCs",
      "Radial glial cells" = "RGCs"
    ),
    force_create_dir = TRUE,
    verbose = TRUE
) {
  # 🧩 Require needed package
  if (!requireNamespace("progress", quietly = TRUE)) {
    install.packages("progress")
  }
  library(progress)
  
  # 📁 Create directory if needed
  if (!dir.exists(output_dir)) {
    if (force_create_dir) {
      dir.create(output_dir, recursive = TRUE)
      if (verbose) message(sprintf("📂 Created output directory: %s", output_dir))
    } else {
      stop(sprintf("❌ Directory does not exist: %s", output_dir))
    }
  }
  
  # 📊 Progress bar (only if not verbose)
  if (!verbose) {
    pb <- progress_bar$new(
      format = "Generating plots [:bar] :percent | ETA: :eta",
      total = length(genes_vector),
      width = 70
    )
  }
  
  # 🔁 Iterate over genes
  for (gene in genes_vector) {
    if (verbose) {
      message(sprintf("▶️ Generating plot for gene: %s", gene))
    } else {
      pb$tick()
    }
    
    output_file <- file.path(output_dir, paste0(gene, "_cellSummary_binary.svg"))
    svg(filename = output_file, width = 15, height = 6.5)
    
    tryCatch({
      seurat_featurePlot_binary_clusterSummary(
        seurat_obj = seurat_obj,
        gene_pattern = gene,
        split.by = split.by,
        split_levels = split_levels,
        cluster_col = cluster_col,
        label_size = 12,
        point_size = 0.6,
        cluster_summary_figure = cluster_summary_figure
      ) %>% print()
    }, error = function(e) {
      message(sprintf("⚠️ Error while plotting %s: %s", gene, e$message))
    })
    
    dev.off()
  }
  
  if (verbose) {
    message("✅ All SVG plots generated successfully.")
  } else {
    pb$terminate()
    message("\n✅ All SVG plots generated successfully.")
  }
}