#' Perform clustering for one or multiple resolutions and optional UMAP visualization
#'
#' Runs FindNeighbors, FindClusters (supports multiple resolutions) and optionally RunUMAP.
#' For multiple resolutions, metadata columns are renamed to seurat_clusters_resX format.
seurat_multiRes_clustering_umap <- function(
    seurat_obj,
    dims = 1:20,
    resolution = NULL,         # NULL = Seurat default 0.8
    k.param = 20,
    nn.method = "rann",
    algorithm = 1,
    reduction = "pca",
    do_umap = TRUE,
    umap_n_neighbors = 30,
    umap_min_dist = 0.3,
    umap_spread = 1.0,
    umap_metric = "cosine",
    seed = 777,
    verbose = TRUE
) {
  # Check if the specified reduction exists
  if (!reduction %in% names(seurat_obj@reductions)) {
    stop(sprintf("Reduction '%s' not found. Run PCA or provide another reduction.", reduction))
  }
  
  # Adjust dims if the requested max dimension exceeds available components
  max_dim <- ncol(Embeddings(seurat_obj, reduction = reduction))
  if (length(dims) == 0 || max(dims) > max_dim) {
    dims <- seq_len(max_dim)
  }
  
  # Set default resolution if NULL
  if (is.null(resolution)) {
    resolution <- 0.8
  }
  
  # 1) Find nearest neighbors
  seurat_obj <- FindNeighbors(
    seurat_obj,
    dims = dims,
    reduction = reduction,
    k.param = k.param,
    nn.method = nn.method,
    verbose = verbose
  )
  
  # 2) Cluster cells for one or multiple resolutions
  seurat_obj <- FindClusters(
    seurat_obj,
    resolution = resolution,
    algorithm = algorithm,
    random.seed = seed,
    verbose = verbose
  )
  
  # 3) Rename cluster metadata columns
  for (res in resolution) {
    old_name <- paste0(DefaultAssay(seurat_obj), "_snn_res.", res)
    new_name <- paste0("seurat_clusters_res", res)
    if (old_name %in% colnames(seurat_obj@meta.data)) {
      colnames(seurat_obj@meta.data)[colnames(seurat_obj@meta.data) == old_name] <- new_name
    }
  }
  
  # 4) Run UMAP if requested
  if (do_umap) {
    seurat_obj <- RunUMAP(
      seurat_obj,
      dims = dims,
      reduction = reduction,
      n.neighbors = umap_n_neighbors,
      min.dist = umap_min_dist,
      spread = umap_spread,
      metric = umap_metric,
      verbose = verbose,
      seed.use = seed
    )
  }
  
  return(seurat_obj)
}



#' Clustering for one/many resolutions + optional UMAP (idempotent)
seurat_multiRes_clustering_umap <- function(
    seurat_obj,
    dims = 1:20,
    resolution = NULL,          # NULL -> 0.8
    k.param = 20,
    nn.method = "rann",
    algorithm = 1,
    reduction = "pca",
    do_umap = TRUE,
    umap_n_neighbors = 30,
    umap_min_dist = 0.3,
    umap_spread = 1.0,
    umap_metric = "cosine",
    seed = 0,
    verbose = TRUE,
    recompute_neighbors = FALSE,  # <- skip FindNeighbors if graphs already exist
    recompute_umap = FALSE,       # <- skip RunUMAP if UMAP already exists
    set_ident_to_first_res = TRUE # <- set Idents() to the first resolution present
) {
  if (is.null(resolution)) resolution <- 0.8
  if (!reduction %in% names(seurat_obj@reductions)) {
    stop(sprintf("Reduction '%s' not found. Run PCA (or specify another).", reduction))
  }
  
  # normalize resolution labels e.g. 0.05 -> "0.05", 1 -> "1"
  fmt <- function(x) sub("\\.?0+$", "", formatC(x, format = "f", digits = 3))
  
  # ensure dims are valid
  max_dim <- ncol(Embeddings(seurat_obj, reduction = reduction))
  if (length(dims) == 0 || max(dims) > max_dim) dims <- seq_len(max_dim)
  
  # 1) FindNeighbors (recompute only if asked or missing)
  snn_name <- paste0(DefaultAssay(seurat_obj), "_snn")
  graphs_have <- snn_name %in% names(seurat_obj@graphs)
  if (recompute_neighbors || !graphs_have) {
    seurat_obj <- FindNeighbors(
      seurat_obj,
      dims = dims,
      reduction = reduction,
      k.param = k.param,
      nn.method = nn.method,
      verbose = verbose
    )
  } else if (verbose) {
    message("Skipping FindNeighbors(): existing graph found (", snn_name, ").")
  }
  
  # 2) Run FindClusters ONLY for missing resolutions
  target_cols <- paste0("seurat_clusters_res", fmt(resolution))
  have_cols <- colnames(seurat_obj@meta.data)
  missing_idx <- which(!(target_cols %in% have_cols))
  
  if (length(missing_idx) > 0) {
    res_to_run <- resolution[missing_idx]
    # run clustering (Seurat will write default names like RNA_snn_res.0.3)
    seurat_obj <- FindClusters(
      seurat_obj,
      resolution = res_to_run,
      algorithm = algorithm,
      random.seed = seed,
      verbose = verbose
    )
    # rename newly created columns to our target naming
    for (res in res_to_run) {
      old1 <- paste0(DefaultAssay(seurat_obj), "_snn_res.", fmt(res))
      old2 <- paste0(DefaultAssay(seurat_obj), "_snn_res.", res) # fallback
      new  <- paste0("seurat_clusters_res", fmt(res))
      if (old1 %in% colnames(seurat_obj@meta.data)) {
        colnames(seurat_obj@meta.data)[colnames(seurat_obj@meta.data) == old1] <- new
      } else if (old2 %in% colnames(seurat_obj@meta.data)) {
        colnames(seurat_obj@meta.data)[colnames(seurat_obj@meta.data) == old2] <- new
      }
    }
  } else if (verbose) {
    message("All requested resolutions already present (", paste(target_cols, collapse = ", "), ").")
  }
  
  # 3) Set Idents to the first requested/available resolution
  if (set_ident_to_first_res) {
    first_col <- target_cols[1]
    if (first_col %in% colnames(seurat_obj@meta.data)) {
      Idents(seurat_obj) <- seurat_obj@meta.data[[first_col]]
    }
  }
  
  # 4) UMAP (recompute only if asked or missing)
  if (do_umap) {
    if (recompute_umap || !("umap" %in% names(seurat_obj@reductions))) {
      seurat_obj <- RunUMAP(
        seurat_obj,
        dims = dims,
        reduction = reduction,
        n.neighbors = umap_n_neighbors,
        min.dist = umap_min_dist,
        spread = umap_spread,
        metric = umap_metric,
        verbose = verbose,
        seed.use = seed
      )
    } else if (verbose) {
      message("Skipping RunUMAP(): reduction 'umap' already exists.")
    }
  }
  
  seurat_obj
}



seurat_multiRes_clustering_umap <- function(
    seurat_obj,
    dims = 1:20,
    resolution = NULL,          # NULL -> 0.8
    k.param = 20,
    nn.method = "rann",
    algorithm = 1,
    reduction = "pca",
    do_umap = TRUE,
    umap_n_neighbors = 30,
    umap_min_dist = 0.3,
    umap_spread = 1.0,
    umap_metric = "cosine",
    umap_n_threads = NULL,      # <--- nowość: wątki UMAP
    seed = 0,
    verbose = TRUE,
    recompute_neighbors = FALSE,
    recompute_umap = FALSE,
    set_ident_to_first_res = TRUE,
    assay_for_graphs = NULL     # <--- nowość: np. "integrated" albo "RNA"
) {
  if (is.null(resolution)) resolution <- 0.8
  if (!reduction %in% names(seurat_obj@reductions)) {
    stop(sprintf("Reduction '%s' not found. Run PCA (or specify another).", reduction))
  }
  
  # Ustal prefix grafów wg wybranego assay (domyślnie bierze się z DefaultAssay)
  if (is.null(assay_for_graphs)) assay_for_graphs <- DefaultAssay(seurat_obj)
  nn_names <- paste0(assay_for_graphs, c("_nn", "_snn"))  # c("..._nn","..._snn")
  
  # upewnij się, że zakres dims jest OK
  max_dim <- ncol(Embeddings(seurat_obj, reduction = reduction))
  if (length(dims) == 0 || max(dims) > max_dim) dims <- seq_len(max_dim)
  
  # 1) FindNeighbors (recompute tylko gdy trzeba)
  snn_name <- nn_names[2]
  graphs_have <- snn_name %in% names(seurat_obj@graphs)
  if (recompute_neighbors || !graphs_have) {
    seurat_obj <- FindNeighbors(
      seurat_obj,
      dims = dims,
      reduction = reduction,
      k.param = k.param,
      nn.method = nn.method,
      graph.name = nn_names,       # <--- jawnie nazwij grafy
      verbose = verbose
    )
  } else if (verbose) {
    message("Skipping FindNeighbors(): existing graph found (", snn_name, ").")
  }
  
  # 2) FindClusters tylko dla brakujących rozdzielczości
  fmt <- function(x) sub("\\.?0+$", "", formatC(x, format = "f", digits = 3))
  target_cols <- paste0("seurat_clusters_res", fmt(resolution))
  have_cols <- colnames(seurat_obj@meta.data)
  missing_idx <- which(!(target_cols %in% have_cols))
  
  if (length(missing_idx) > 0) {
    res_to_run <- resolution[missing_idx]
    seurat_obj <- FindClusters(
      seurat_obj,
      resolution = res_to_run,
      algorithm = algorithm,
      random.seed = seed,
      graph.name = snn_name,       # <--- klastruj po tym konkretnym grafie
      verbose = verbose
    )
    # Seurat tworzy kolumny w stylu "<assay>_snn_res.X" – przenieśmy je do jednolitej nazwy
    for (res in res_to_run) {
      old1 <- paste0(snn_name, "_res.", fmt(res))
      old2 <- paste0(snn_name, "_res.", res) # fallback
      new  <- paste0("seurat_clusters_res", fmt(res))
      if (old1 %in% colnames(seurat_obj@meta.data)) {
        colnames(seurat_obj@meta.data)[colnames(seurat_obj@meta.data) == old1] <- new
      } else if (old2 %in% colnames(seurat_obj@meta.data)) {
        colnames(seurat_obj@meta.data)[colnames(seurat_obj@meta.data) == old2] <- new
      }
    }
  } else if (verbose) {
    message("All requested resolutions already present (", paste(target_cols, collapse = ", "), ").")
  }
  
  # 3) Ustaw Idents na pierwszą rozdzielczość (opcjonalnie podbij też seurat_clusters)
  if (set_ident_to_first_res) {
    first_col <- target_cols[1]
    if (first_col %in% colnames(seurat_obj@meta.data)) {
      Idents(seurat_obj) <- seurat_obj@meta.data[[first_col]]
      seurat_obj$seurat_clusters <- Idents(seurat_obj)
    }
  }
  
  # 4) UMAP (tylko gdy trzeba)
  if (do_umap) {
    if (recompute_umap || !("umap" %in% names(seurat_obj@reductions))) {
      args <- list(
        object      = seurat_obj,
        dims        = dims,
        reduction   = reduction,
        n.neighbors = umap_n_neighbors,
        min.dist    = umap_min_dist,
        spread      = umap_spread,
        metric      = umap_metric,
        verbose     = verbose,
        seed.use    = seed
      )
      if (!is.null(umap_n_threads)) args$n.threads <- umap_n_threads
      seurat_obj <- do.call(RunUMAP, args)
    } else if (verbose) {
      message("Skipping RunUMAP(): reduction 'umap' already exists.")
    }
  }
  
  seurat_obj
}