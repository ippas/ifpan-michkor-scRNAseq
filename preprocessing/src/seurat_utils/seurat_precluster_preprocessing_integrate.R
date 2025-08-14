# Integracja (LogNormalize/SCT) -> ScaleData -> RunPCA
# Parallel & RAM limit ogarniam wewnątrz (future), a na końcu przywracam poprzedni plan.
seurat_precluster_preprocessing_integrate <- function(
    seurat_obj,
    split_by,                         # np. "sample_id" albo "sample"
    method = c("LogNormalize", "SCT"),
    nfeatures = 2000,
    pca_npcs = 50,
    vars_to_regress = NULL,           # po integracji; np. c("nCount_RNA")
    parallel = TRUE,
    n_cores = 40,
    future_plan = c("multisession","multicore","sequential"),
    maxsize_gb = 50,                  # limit dla future.globals.maxSize
    blas_threads = NULL,              # np. 40; wymaga RhpcBLASctl
    verbose = TRUE,
    seed = 777
) {
  stopifnot(split_by %in% colnames(seurat_obj@meta.data))
  method <- match.arg(method)
  future_plan <- match.arg(future_plan)
  set.seed(seed)
  
  # --- parallel setup (z przywróceniem stanu po wyjściu) ---
  old_plan <- NULL; old_max <- getOption("future.globals.maxSize")
  if (parallel) {
    if (!requireNamespace("future", quietly = TRUE)) stop("Zainstaluj pakiet 'future'.")
    old_plan <- future::plan()
    on.exit({
      if (!is.null(old_plan)) future::plan(old_plan)
      options(future.globals.maxSize = old_max)
    }, add = TRUE)
    options(future.globals.maxSize = maxsize_gb * 1024^3)
    if (future_plan == "multicore") {
      future::plan(future::multicore, workers = n_cores)
    } else if (future_plan == "multisession") {
      future::plan(future::multisession, workers = n_cores)
    } else {
      future::plan(future::sequential)
    }
  }
  
  if (!is.null(blas_threads) && requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(blas_threads)
    RhpcBLASctl::omp_set_num_threads(blas_threads)
  }
  
  # --- 1) Split ---
  obj_list <- Seurat::SplitObject(seurat_obj, split.by = split_by)
  
  # --- 2) Per-sample normalizacja + HVG i integracja ---
  use_flapply <- parallel && requireNamespace("future.apply", quietly = TRUE)
  apply_fun <- if (use_flapply) future.apply::future_lapply else lapply
  
  if (method == "LogNormalize") {
    if (verbose) message("[Integrate] LogNormalize + VST HVG per sample")
    obj_list <- apply_fun(obj_list, function(x) {
      x <- Seurat::NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = verbose)
      x <- Seurat::FindVariableFeatures(x, selection.method = "vst", nfeatures = nfeatures, verbose = verbose)
      x
    })
    features <- Seurat::SelectIntegrationFeatures(object.list = obj_list, nfeatures = nfeatures)
    anchors  <- Seurat::FindIntegrationAnchors(object.list = obj_list, anchor.features = features, verbose = verbose)
    obj_int  <- Seurat::IntegrateData(anchorset = anchors, verbose = verbose)
    
  } else { # SCT
    if (verbose) message("[Integrate] SCTransform per sample + SCT integration")
    obj_list <- apply_fun(obj_list, function(x) Seurat::SCTransform(x, verbose = verbose))
    features <- Seurat::SelectIntegrationFeatures(object.list = obj_list, nfeatures = nfeatures)
    obj_list <- Seurat::PrepSCTIntegration(object.list = obj_list, anchor.features = features, verbose = verbose)
    anchors  <- Seurat::FindIntegrationAnchors(object.list = obj_list,
                                               normalization.method = "SCT",
                                               anchor.features = features, verbose = verbose)
    obj_int  <- Seurat::IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = verbose)
  }
  
  # --- 3) Po integracji: ScaleData (+ ewentualna regresja) i PCA (BEZ NN/UMAP/klastrów) ---
  Seurat::DefaultAssay(obj_int) <- "integrated"
  
  # helper: one-hot dla czynników do regresji
  prepare_regressors <- function(obj, vars) {
    if (is.null(vars) || length(vars) == 0) return(list(obj = obj, regressors = NULL, temp_cols = character(0)))
    vars <- as.character(vars)
    miss <- vars[!(vars %in% colnames(obj@meta.data))]
    if (length(miss) > 0) stop("vars_to_regress nie ma w meta.data: ", paste(miss, collapse = ", "))
    temp_cols <- character(0); regs <- character(0)
    for (v in vars) {
      col <- obj@meta.data[[v]]
      if (is.numeric(col)) {
        regs <- c(regs, v)
      } else if (is.factor(col) || is.character(col)) {
        f <- if (is.factor(col)) col else factor(col)
        mm <- stats::model.matrix(~ f - 1)
        base <- paste0(v, "_oh_")
        newn <- paste0(base, gsub("^f", "", colnames(mm)))
        obj@meta.data[, newn] <- mm
        temp_cols <- c(temp_cols, newn); regs <- c(regs, newn)
        if (verbose) message("One-hot: ", v, " -> ", paste(newn, collapse = ", "))
      } else warning("Pomijam '", v, "': nie numeric/factor/character.")
    }
    list(obj = obj, regressors = unique(regs), temp_cols = temp_cols)
  }
  
  prep <- prepare_regressors(obj_int, vars_to_regress)
  obj_int <- prep$obj
  
  all_genes <- rownames(obj_int)
  if (is.null(prep$regressors)) {
    obj_int <- Seurat::ScaleData(obj_int, features = all_genes, verbose = verbose)
  } else {
    obj_int <- Seurat::ScaleData(obj_int, features = all_genes, vars.to.regress = prep$regressors, verbose = verbose)
    if (length(prep$temp_cols) > 0) obj_int@meta.data[, prep$temp_cols] <- NULL
  }
  
  vfeat <- Seurat::VariableFeatures(obj_int)
  if (length(vfeat) > 0) {
    obj_int <- Seurat::RunPCA(obj_int, features = vfeat, npcs = pca_npcs, verbose = verbose, seed.use = seed)
  } else {
    obj_int <- Seurat::RunPCA(obj_int, npcs = pca_npcs, verbose = verbose, seed.use = seed)
  }
  
  return(obj_int)
}


# Integracja (LogNormalize/SCT) -> ScaleData -> RunPCA
# Parallel & RAM limit wewnątrz; na końcu przywracamy ustawienia. Z pomiarem czasu.
seurat_precluster_preprocessing_integrate <- function(
    seurat_obj,
    split_by,                         # np. "sample_id" albo "sample"
    method = c("LogNormalize", "SCT"),
    nfeatures = 2000,
    pca_npcs = 50,
    vars_to_regress = NULL,           # po integracji; np. c("nCount_RNA")
    parallel = TRUE,
    n_cores = 40,
    future_plan = c("multisession","multicore","sequential"),
    maxsize_gb = 50,                  # limit dla future.globals.maxSize
    blas_threads = NULL,              # np. 40; wymaga RhpcBLASctl
    verbose = TRUE,
    seed = 777,
    record_timing = TRUE              # <--- NOWE: zapis i log czasu
) {
  stopifnot(split_by %in% colnames(seurat_obj@meta.data))
  method <- match.arg(method)
  future_plan <- match.arg(future_plan)
  set.seed(seed)
  
  t0 <- Sys.time()  # start
  
  # --- parallel setup (z przywróceniem stanu po wyjściu) ---
  old_plan <- NULL; old_max <- getOption("future.globals.maxSize")
  if (parallel) {
    if (!requireNamespace("future", quietly = TRUE)) stop("Zainstaluj pakiet 'future'.")
    old_plan <- future::plan()
    on.exit({
      if (!is.null(old_plan)) future::plan(old_plan)
      options(future.globals.maxSize = old_max)
    }, add = TRUE)
    options(future.globals.maxSize = maxsize_gb * 1024^3)
    if (future_plan == "multicore") {
      future::plan(future::multicore, workers = n_cores)
    } else if (future_plan == "multisession") {
      future::plan(future::multisession, workers = n_cores)
    } else {
      future::plan(future::sequential)
    }
  }
  
  if (!is.null(blas_threads) && requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(blas_threads)
    RhpcBLASctl::omp_set_num_threads(blas_threads)
  }
  
  # --- 1) Split ---
  t_split0 <- Sys.time()
  obj_list <- Seurat::SplitObject(seurat_obj, split.by = split_by)
  t_split1 <- Sys.time()
  
  # --- 2) Per-sample normalizacja + HVG i integracja ---
  use_flapply <- parallel && requireNamespace("future.apply", quietly = TRUE)
  apply_fun <- if (use_flapply) future.apply::future_lapply else lapply
  
  t_norm0 <- Sys.time()
  if (method == "LogNormalize") {
    if (verbose) message("[Integrate] LogNormalize + VST HVG per sample")
    obj_list <- apply_fun(obj_list, function(x) {
      x <- Seurat::NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = verbose)
      x <- Seurat::FindVariableFeatures(x, selection.method = "vst", nfeatures = nfeatures, verbose = verbose)
      x
    })
    features <- Seurat::SelectIntegrationFeatures(object.list = obj_list, nfeatures = nfeatures)
    t_int0 <- Sys.time()
    anchors  <- Seurat::FindIntegrationAnchors(object.list = obj_list, anchor.features = features, verbose = verbose)
    obj_int  <- Seurat::IntegrateData(anchorset = anchors, verbose = verbose)
    t_int1 <- Sys.time()
  } else { # SCT
    if (verbose) message("[Integrate] SCTransform per sample + SCT integration")
    obj_list <- apply_fun(obj_list, function(x) Seurat::SCTransform(x, verbose = verbose))
    features <- Seurat::SelectIntegrationFeatures(object.list = obj_list, nfeatures = nfeatures)
    t_int0 <- Sys.time()
    obj_list <- Seurat::PrepSCTIntegration(object.list = obj_list, anchor.features = features, verbose = verbose)
    anchors  <- Seurat::FindIntegrationAnchors(object.list = obj_list,
                                               normalization.method = "SCT",
                                               anchor.features = features, verbose = verbose)
    obj_int  <- Seurat::IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = verbose)
    t_int1 <- Sys.time()
  }
  t_norm1 <- Sys.time()
  
  # --- 3) Po integracji: ScaleData (+ ew. regresja) i PCA (BEZ NN/UMAP/klastrów) ---
  Seurat::DefaultAssay(obj_int) <- "integrated"
  
  # helper: one-hot dla czynników do regresji
  prepare_regressors <- function(obj, vars) {
    if (is.null(vars) || length(vars) == 0) return(list(obj = obj, regressors = NULL, temp_cols = character(0)))
    vars <- as.character(vars)
    miss <- vars[!(vars %in% colnames(obj@meta.data))]
    if (length(miss) > 0) stop("vars_to_regress nie ma w meta.data: ", paste(miss, collapse = ", "))
    temp_cols <- character(0); regs <- character(0)
    for (v in vars) {
      col <- obj@meta.data[[v]]
      if (is.numeric(col)) {
        regs <- c(regs, v)
      } else if (is.factor(col) || is.character(col)) {
        f <- if (is.factor(col)) col else factor(col)
        mm <- stats::model.matrix(~ f - 1)
        base <- paste0(v, "_oh_")
        newn <- paste0(base, gsub("^f", "", colnames(mm)))
        obj@meta.data[, newn] <- mm
        temp_cols <- c(temp_cols, newn); regs <- c(regs, newn)
        if (verbose) message("One-hot: ", v, " -> ", paste(newn, collapse = ", "))
      } else warning("Pomijam '", v, "': nie numeric/factor/character.")
    }
    list(obj = obj, regressors = unique(regs), temp_cols = temp_cols)
  }
  
  prep <- prepare_regressors(obj_int, vars_to_regress)
  obj_int <- prep$obj
  
  t_scale0 <- Sys.time()
  all_genes <- rownames(obj_int)
  if (is.null(prep$regressors)) {
    obj_int <- Seurat::ScaleData(obj_int, features = all_genes, verbose = verbose)
  } else {
    obj_int <- Seurat::ScaleData(obj_int, features = all_genes, vars.to.regress = prep$regressors, verbose = verbose)
    if (length(prep$temp_cols) > 0) obj_int@meta.data[, prep$temp_cols] <- NULL
  }
  t_scale1 <- Sys.time()
  
  t_pca0 <- Sys.time()
  vfeat <- Seurat::VariableFeatures(obj_int)
  if (length(vfeat) > 0) {
    obj_int <- Seurat::RunPCA(obj_int, features = vfeat, npcs = pca_npcs, verbose = verbose, seed.use = seed)
  } else {
    obj_int <- Seurat::RunPCA(obj_int, npcs = pca_npcs, verbose = verbose, seed.use = seed)
  }
  t_pca1 <- Sys.time()
  
  t1 <- Sys.time()  # end
  
  if (record_timing) {
    timing <- list(
      start   = t0,
      end     = t1,
      seconds = as.numeric(difftime(t1, t0, units = "secs")),
      breakdown = list(
        split_seconds   = as.numeric(difftime(t_split1, t_split0, units = "secs")),
        per_sample_norm_hvg_seconds = as.numeric(difftime(t_norm1, t_norm0, units = "secs")),
        integrate_seconds = as.numeric(difftime(t_int1, t_int0, units = "secs")),
        scale_seconds   = as.numeric(difftime(t_scale1, t_scale0, units = "secs")),
        pca_seconds     = as.numeric(difftime(t_pca1, t_pca0, units = "secs"))
      ),
      method  = method,
      split_by = split_by,
      nfeatures = nfeatures,
      pca_npcs = pca_npcs,
      parallel = parallel,
      n_cores  = n_cores,
      future_plan = future_plan,
      maxsize_gb  = maxsize_gb
    )
    if (is.null(obj_int@misc)) obj_int@misc <- list()
    obj_int@misc$precluster_integrate_timing <- timing
    
    if (verbose) {
      pretty_min <- function(s) sprintf("%.1f min", s/60)
      message(sprintf(
        "[Timing] total=%s | split=%s | per-sample=%s | integrate=%s | scale=%s | pca=%s",
        pretty_min(timing$seconds),
        pretty_min(timing$breakdown$split_seconds),
        pretty_min(timing$breakdown$per_sample_norm_hvg_seconds),
        pretty_min(timing$breakdown$integrate_seconds),
        pretty_min(timing$breakdown$scale_seconds),
        pretty_min(timing$breakdown$pca_seconds)
      ))
    }
  }
  
  return(obj_int)
}


# Integracja w Seurat v5 z warstwami
# Seurat v5: Integracja (LogNormalize/SCT) -> ScaleData -> RunPCA
seuratV5_precluster_preprocessing_integrate <- function(
    seurat_obj, split_by, assay = NULL,
    method = c("LogNormalize","SCT"),
    nfeatures = 2000, pca_npcs = 50,
    vars_to_regress = NULL,
    parallel = TRUE, n_cores = 32, future_plan = c("multisession","multicore","sequential"),
    maxsize_gb = 100, blas_threads = NULL, dt_threads = NULL,
    verbose = TRUE, 
    verbose2 = TRUE, seed = 777
) {
  method <- match.arg(method); future_plan <- match.arg(future_plan)
  if (is.null(assay)) assay <- DefaultAssay(seurat_obj)
  set.seed(seed)
  
  say <- function(msg) if (verbose) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), msg))
  
  # parallel setup (jak wcześniej)
  old_plan <- NULL; old_max <- getOption("future.globals.maxSize")
  if (parallel) {
    if (!requireNamespace("future", quietly = TRUE)) stop("Zainstaluj 'future'")
    old_plan <- future::plan()
    on.exit({ if (!is.null(old_plan)) future::plan(old_plan); options(future.globals.maxSize = old_max) }, add = TRUE)
    options(future.globals.maxSize = maxsize_gb * 1024^3)
    if (future_plan == "multicore") future::plan(future::multicore, workers = n_ores <- n_cores) else
      if (future_plan == "multisession") future::plan(future::multisession, workers = n_cores) else
        future::plan(future::sequential)
  }
  if (!is.null(blas_threads) && requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(blas_threads); RhpcBLASctl::omp_set_num_threads(blas_threads)
  }
  if (!is.null(dt_threads) && requireNamespace("data.table", quietly = TRUE)) {
    data.table::setDTthreads(dt_threads)
  }
  
  # 1) Split
  say(sprintf("SplitObject by '%s'", split_by))
  obj_list <- Seurat::SplitObject(seurat_obj, split.by = split_by)
  obj_list <- lapply(obj_list, function(x) { DefaultAssay(x) <- assay; x })
  
  # 1a) JoinLayers
  say("JoinLayers (if multiple layers present)")
  obj_list <- lapply(obj_list, function(x) {
    ax <- x[[assay]]
    if (inherits(ax, "Assay5")) {
      ln <- SeuratObject::Layers(ax)   # <-- tu jest zmiana
      if (length(ln) > 1) x <- Seurat::JoinLayers(x, assay = assay)
    }
    x
  })
  
  # 2) Per-sample normalize + HVG
  use_flapply <- parallel && requireNamespace("future.apply", quietly = TRUE)
  apply_fun <- if (use_flapply) future.apply::future_lapply else lapply
  
  if (method == "LogNormalize") {
    say("Per-sample: LogNormalize + HVG (layer='counts')")
    obj_list <- apply_fun(obj_list, function(x) {
      x <- Seurat::NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 1e4,
                                 layer = "counts", verbose = verbose)
      x <- Seurat::FindVariableFeatures(x, selection.method = "vst", nfeatures = nfeatures,
                                        layer = "counts", verbose = verbose)
      x
    })
    features <- Seurat::SelectIntegrationFeatures(object.list = obj_list, nfeatures = nfeatures)
    say("FindIntegrationAnchors ...")
    anchors  <- Seurat::FindIntegrationAnchors(object.list = obj_list, anchor.features = features, verbose = verbose)
    say("IntegrateData ...")
    obj_int  <- Seurat::IntegrateData(anchorset = anchors, verbose = verbose)
    
  } else { # SCT
    say("Per-sample: SCTransform (layer='counts') + SCT integration")
    obj_list <- apply_fun(obj_list, function(x) Seurat::SCTransform(x, layer = "counts", verbose = verbose2))
    features <- Seurat::SelectIntegrationFeatures(object.list = obj_list, nfeatures = nfeatures)
    obj_list <- Seurat::PrepSCTIntegration(object.list = obj_list, anchor.features = features, verbose = verbose2)
    anchors  <- Seurat::FindIntegrationAnchors(object.list = obj_list, normalization.method = "SCT",
                                               anchor.features = features, verbose = verbose2)
    obj_int  <- Seurat::IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = verbose2)
  }
  
  # 3) Scale + PCA
  Seurat::DefaultAssay(obj_int) <- "integrated"
  
  # (regresja jak wcześniej)
  prepare_regressors <- function(obj, vars) {
    if (is.null(vars) || length(vars) == 0) return(list(obj = obj, regressors = NULL, temp_cols = character(0)))
    vars <- as.character(vars); miss <- vars[!(vars %in% colnames(obj@meta.data))]
    if (length(miss) > 0) stop("vars_to_regress nie ma w meta.data: ", paste(miss, collapse = ", "))
    temp_cols <- character(0); regs <- character(0)
    for (v in vars) {
      col <- obj@meta.data[[v]]
      if (is.numeric(col)) regs <- c(regs, v) else if (is.factor(col) || is.character(col)) {
        f <- if (is.factor(col)) col else factor(col)
        mm <- stats::model.matrix(~ f - 1); base <- paste0(v, "_oh_")
        newn <- paste0(base, gsub("^f", "", colnames(mm)))
        obj@meta.data[, newn] <- mm; temp_cols <- c(temp_cols, newn); regs <- c(regs, newn)
        say(sprintf("One-hot: %s -> %s", v, paste(newn, collapse = ", ")))
      } else warning("Pomijam '", v, "': nie numeric/factor/character.")
    }
    list(obj = obj, regressors = unique(regs), temp_cols = temp_cols)
  }
  prep <- prepare_regressors(obj_int, vars_to_regress); obj_int <- prep$obj
  
  say("ScaleData on integrated")
  all_genes <- rownames(obj_int)
  if (is.null(prep$regressors)) {
    obj_int <- Seurat::ScaleData(obj_int, features = all_genes, verbose = verbose2)
  } else {
    obj_int <- Seurat::ScaleData(obj_int, features = all_genes, vars.to.regress = prep$regressors, verbose = verbose2)
    if (length(prep$temp_cols) > 0) obj_int@meta.data[, prep$temp_cols] <- NULL
  }
  
  say(sprintf("RunPCA (npcs=%d)", pca_npcs))
  vfeat <- Seurat::VariableFeatures(obj_int)
  if (length(vfeat) > 0) obj_int <- Seurat::RunPCA(obj_int, features = vfeat, npcs = pca_npcs, verbose = verbose2, seed.use = seed)
  else                   obj_int <- Seurat::RunPCA(obj_int, npcs = pca_npcs, verbose = verbose2, seed.use = seed)
  
  obj_int
}
