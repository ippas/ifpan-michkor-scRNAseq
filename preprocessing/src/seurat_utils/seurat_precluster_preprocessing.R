# src/seurat_utils/preprocess_seurat_data_v1.R

# preprocess_seurat_data_v1 <- function(seurat_obj,
#                                       n_variable_features = 2000,
#                                       selection_method = "vst") {
#   seurat_obj <- NormalizeData(seurat_obj)
#   seurat_obj <- FindVariableFeatures(seurat_obj,
#                                      selection.method = selection_method,
#                                      nfeatures = n_variable_features)
#   all.genes <- rownames(seurat_obj)
#   seurat_obj <- ScaleData(seurat_obj, features = all.genes)
#   seurat_obj <- RunPCA(seurat_obj,
#                        features = VariableFeatures(seurat_obj))
#   return(seurat_obj)
# }

seurat_precluster_preprocessing <- function(
    seurat_obj,
    n_variable_features = 2000,
    selection_method = "vst",
    normalization_method = "LogNormalize",
    scale_factor = 1e4,
    vars_to_regress = NULL,
    verbose = FALSE,
    seed = 777
) {
  set.seed(seed)
  seurat_obj <- NormalizeData(
    seurat_obj,
    normalization.method = normalization_method,
    scale.factor = scale_factor,
    verbose = verbose
  )
  
  seurat_obj <- FindVariableFeatures(
    seurat_obj,
    selection.method = selection_method,
    nfeatures = n_variable_features,
    verbose = verbose
  )
  
  seurat_obj <- ScaleData(
    seurat_obj,
    features = rownames(seurat_obj),
    vars.to.regress = vars_to_regress,
    verbose = verbose
  )
  
  seurat_obj <- RunPCA(
    seurat_obj,
    features = VariableFeatures(seurat_obj),
    verbose = verbose
  )
  
  return(seurat_obj)
}


seurat_precluster_preprocessing <- function(
    seurat_obj,
    n_variable_features = 2000,
    selection_method = "vst",
    normalization_method = "LogNormalize",
    scale_factor = 1e4,
    vars_to_regress = NULL,   # NULL -> no regression
    verbose = FALSE,
    seed = 777
) {
  set.seed(seed)
  
  # 1) Normalize
  seurat_obj <- NormalizeData(
    seurat_obj,
    normalization.method = normalization_method,
    scale.factor = scale_factor,
    verbose = verbose
  )
  
  # 2) Variable features
  seurat_obj <- FindVariableFeatures(
    seurat_obj,
    selection.method = selection_method,
    nfeatures = n_variable_features,
    verbose = verbose
  )
  
  # Helper: prepare regression variables
  prepare_regressors <- function(obj, vars) {
    if (is.null(vars) || length(vars) == 0) {
      return(list(obj = obj, regressors = NULL, temp_cols = character(0)))
    }
    vars <- as.character(vars)
    
    # check presence
    missing <- vars[!(vars %in% colnames(obj@meta.data))]
    if (length(missing) > 0) {
      stop("vars_to_regress not found in meta.data: ", paste(missing, collapse = ", "))
    }
    
    temp_cols <- character(0)
    regressors <- character(0)
    
    for (v in vars) {
      col <- obj@meta.data[[v]]
      if (is.numeric(col)) {
        regressors <- c(regressors, v)
      } else if (is.factor(col) || is.character(col)) {
        # one-hot encode factor/character
        f <- if (is.factor(col)) col else factor(col)
        mm <- stats::model.matrix(~ f - 1)  # no intercept, one column per level
        base <- paste0(v, "_oh_")
        new_names <- paste0(base, gsub("^f", "", colnames(mm)))
        # append to meta.data
        obj@meta.data[, new_names] <- mm
        temp_cols <- c(temp_cols, new_names)
        regressors <- c(regressors, new_names)
        if (verbose) message("One-hot encoded '", v, "' into: ", paste(new_names, collapse = ", "))
      } else {
        warning("Variable '", v, "' is neither numeric nor factor/character; skipping.")
      }
    }
    list(obj = obj, regressors = unique(regressors), temp_cols = temp_cols)
  }
  
  prep <- prepare_regressors(seurat_obj, vars_to_regress)
  seurat_obj <- prep$obj
  
  # 3) Scale (with or without regression)
  all_genes <- rownames(seurat_obj)
  if (is.null(prep$regressors)) {
    seurat_obj <- ScaleData(
      seurat_obj,
      features = all_genes,
      verbose = verbose
    )
  } else {
    seurat_obj <- ScaleData(
      seurat_obj,
      features = all_genes,
      vars.to.regress = prep$regressors,
      verbose = verbose
    )
    # clean up temporary one-hot cols
    if (length(prep$temp_cols) > 0) {
      seurat_obj@meta.data[, prep$temp_cols] <- NULL
    }
  }
  
  # 4) PCA
  seurat_obj <- RunPCA(
    seurat_obj,
    features = VariableFeatures(seurat_obj),
    verbose = verbose
  )
  
  return(seurat_obj)
}



# seurat_precluster_preprocessing <- function(
#     seurat_obj,
#     n_variable_features = 2000,
#     selection_method = "vst",
#     normalization_method = "LogNormalize",
#     scale_factor = 1e4,
#     vars_to_regress = NULL,     # NULL -> bez regresji
#     hvg_after_regression = TRUE, # << kluczowa zmiana
#     verbose = FALSE,
#     seed = 777
# ) {
#   set.seed(seed)
#   
#   # 1) Normalize
#   seurat_obj <- NormalizeData(
#     seurat_obj,
#     normalization.method = normalization_method,
#     scale.factor = scale_factor,
#     verbose = verbose
#   )
#   
#   # Helper: prepare one-hot regressors
#   prepare_regressors <- function(obj, vars) {
#     if (is.null(vars) || length(vars) == 0) {
#       return(list(obj = obj, regressors = NULL, temp_cols = character(0)))
#     }
#     vars <- as.character(vars)
#     missing <- vars[!(vars %in% colnames(obj@meta.data))]
#     if (length(missing) > 0) stop("vars_to_regress not in meta.data: ", paste(missing, collapse = ", "))
#     
#     temp_cols <- character(0); regressors <- character(0)
#     for (v in vars) {
#       col <- obj@meta.data[[v]]
#       if (is.numeric(col)) {
#         regressors <- c(regressors, v)
#       } else if (is.factor(col) || is.character(col)) {
#         f <- if (is.factor(col)) col else factor(col)
#         mm <- stats::model.matrix(~ f - 1)
#         base <- paste0(v, "_oh_")
#         new_names <- paste0(base, gsub("^f", "", colnames(mm)))
#         obj@meta.data[, new_names] <- mm
#         temp_cols <- c(temp_cols, new_names)
#         regressors <- c(regressors, new_names)
#         if (verbose) message("One-hot '", v, "': ", paste(new_names, collapse = ", "))
#       } else {
#         warning("Variable '", v, "' not numeric/factor/character; skipping.")
#       }
#     }
#     list(obj = obj, regressors = unique(regressors), temp_cols = temp_cols)
#   }
#   
#   prep <- prepare_regressors(seurat_obj, vars_to_regress)
#   seurat_obj <- prep$obj
#   
#   if (isTRUE(hvg_after_regression) && !is.null(prep$regressors)) {
#     # 2A) Najpierw regresja (na wszystkich genach), potem wybór HVG na rezyduach
#     seurat_obj <- ScaleData(
#       seurat_obj,
#       features = rownames(seurat_obj),
#       vars.to.regress = prep$regressors,
#       verbose = verbose
#     )
#     if (length(prep$temp_cols) > 0) seurat_obj@meta.data[, prep$temp_cols] <- NULL
#     
#     # Zrób assay z rezyduów (scale.data) i wybierz HVG na nim
#     resid_assay <- Seurat::CreateAssayObject(
#       data = Seurat::GetAssayData(seurat_obj, slot = "scale.data")
#     )
#     seurat_obj[["resid"]] <- resid_assay
#     DefaultAssay(seurat_obj) <- "resid"
#     
#     seurat_obj <- FindVariableFeatures(
#       seurat_obj,
#       selection.method = selection_method,
#       nfeatures = n_variable_features,
#       verbose = verbose
#     )
#     
#     # PCA na rezyduach
#     seurat_obj <- RunPCA(
#       seurat_obj,
#       features = VariableFeatures(seurat_obj),
#       verbose = verbose
#     )
#     
#     # (opcjonalnie) wróć do RNA jako domyślnego assayu do DE
#     DefaultAssay(seurat_obj) <- "RNA"
#     
#   } else {
#     # 2B) Klasyczna ścieżka: najpierw HVG, potem (ew.) regresja i PCA na HVG z RNA
#     seurat_obj <- FindVariableFeatures(
#       seurat_obj,
#       selection.method = selection_method,
#       nfeatures = n_variable_features,
#       verbose = verbose
#     )
#     
#     if (is.null(prep$regressors)) {
#       seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj), verbose = verbose)
#     } else {
#       seurat_obj <- ScaleData(
#         seurat_obj,
#         features = rownames(seurat_obj),
#         vars.to.regress = prep$regressors,
#         verbose = verbose
#       )
#       if (length(prep$temp_cols) > 0) seurat_obj@meta.data[, prep$temp_cols] <- NULL
#     }
#     
#     seurat_obj <- RunPCA(
#       seurat_obj,
#       features = VariableFeatures(seurat_obj),
#       verbose = verbose
#     )
#   }
#   
#   return(seurat_obj)
# }
