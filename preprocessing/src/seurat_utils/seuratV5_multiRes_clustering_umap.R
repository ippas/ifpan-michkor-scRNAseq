# v5-friendly: NN/SNN po wybranym assay, PCA z podanej redukcji
seuratV5_multiRes_clustering_umap <- function(
    seurat_obj,
    dims = 1:20,
    resolution = NULL,              # jeśli NULL -> 0.8
    k.param = 20,
    nn.method = "rann",             # możesz użyć "annoy" dla szybkości
    algorithm = 1,
    reduction = "pca",
    do_umap = TRUE,
    umap_n_neighbors = 30,
    umap_min_dist = 0.3,
    umap_spread = 1.0,
    umap_metric = "cosine",
    umap_n_threads = NULL,
    seed = 0,
    verbose = TRUE,                 # nasz krótki log
    verbose2 = FALSE,               # gadatliwość funkcji Seurat
    recompute_neighbors = FALSE,
    recompute_umap = FALSE,
    set_ident_to_first_res = TRUE,
    assay_for_graphs = NULL         # np. "integrated" (po integracji) albo "RNA"
) {
  if (is.null(resolution)) resolution <- 0.8
  if (!reduction %in% names(seurat_obj@reductions)) {
    stop(sprintf("Reduction '%s' not found. Run PCA first.", reduction))
  }
  if (is.null(assay_for_graphs)) assay_for_graphs <- DefaultAssay(seurat_obj)
  
  # dims sanity
  max_dim <- ncol(Embeddings(seurat_obj, reduction = reduction))
  if (length(dims) == 0 || max(dims) > max_dim) dims <- seq_len(max_dim)
  
  # nazwy grafów wynikowe (v5 i v4 zgodne)
  nn_name  <- paste0(assay_for_graphs, "_nn")
  snn_name <- paste0(assay_for_graphs, "_snn")
  
  say <- function(x) if (isTRUE(verbose)) message(x)
  
  # 1) FindNeighbors (po wybranym assay; redukcja = PCA)
  have_graph <- snn_name %in% names(seurat_obj@graphs)
  if (recompute_neighbors || !have_graph) {
    say(sprintf("FindNeighbors on reduction='%s', assay='%s'", reduction, assay_for_graphs))
    seurat_obj <- FindNeighbors(
      seurat_obj,
      reduction = reduction,
      dims      = dims,
      k.param   = k.param,
      nn.method = nn.method,
      assay     = assay_for_graphs,   # <--- KLUCZOWE w v5
      verbose   = verbose2
    )
  } else {
    say(sprintf("Skipping FindNeighbors(): graph %s already exists.", snn_name))
  }
  
  # 2) FindClusters tylko dla brakujących rozdzielczości
  fmt <- function(x) sub("\\.?0+$", "", formatC(x, format = "f", digits = 3))
  target_cols <- paste0("seurat_clusters_res", fmt(resolution))
  have_cols <- colnames(seurat_obj@meta.data)
  missing_idx <- which(!(target_cols %in% have_cols))
  
  if (length(missing_idx) > 0) {
    res_to_run <- resolution[missing_idx]
    say(sprintf("FindClusters on graph '%s' for resolutions: %s",
                snn_name, paste(fmt(res_to_run), collapse = ", ")))
    seurat_obj <- FindClusters(
      seurat_obj,
      resolution  = res_to_run,
      algorithm   = algorithm,
      random.seed = seed,
      graph.name  = snn_name,         # <--- jawnie po właściwym grafie
      verbose     = verbose2
    )
    # ujednolić nazwy kolumn
    for (res in res_to_run) {
      old1 <- paste0(snn_name, "_res.", fmt(res))
      old2 <- paste0(snn_name, "_res.", res)
      new  <- paste0("seurat_clusters_res", fmt(res))
      if (old1 %in% colnames(seurat_obj@meta.data)) {
        colnames(seurat_obj@meta.data)[colnames(seurat_obj@meta.data) == old1] <- new
      } else if (old2 %in% colnames(seurat_obj@meta.data)) {
        colnames(seurat_obj@meta.data)[colnames(seurat_obj@meta.data) == old2] <- new
      }
    }
  } else {
    say(sprintf("All requested resolutions already present: %s", paste(target_cols, collapse = ", ")))
  }
  
  # 3) Ustaw Idents na pierwszą rozdzielczość
  if (set_ident_to_first_res) {
    first_col <- target_cols[1]
    if (first_col %in% colnames(seurat_obj@meta.data)) {
      Idents(seurat_obj) <- seurat_obj@meta.data[[first_col]]
      seurat_obj$seurat_clusters <- Idents(seurat_obj)
    }
  }
  
  # 4) UMAP (na PCA)
  if (do_umap) {
    if (recompute_umap || !("umap" %in% names(seurat_obj@reductions))) {
      say("RunUMAP on PCA")
      args <- list(
        object      = seurat_obj,
        reduction   = reduction,
        dims        = dims,
        n.neighbors = umap_n_neighbors,
        min.dist    = umap_min_dist,
        spread      = umap_spread,
        metric      = umap_metric,
        verbose     = verbose2,
        seed.use    = seed
      )
      if (!is.null(umap_n_threads)) args$n.threads <- umap_n_threads
      seurat_obj <- do.call(RunUMAP, args)
    } else {
      say("Skipping RunUMAP(): 'umap' already exists.")
    }
  }
  
  seurat_obj
}
