
# 0) Kopia i identyfikator próbki
merged_pilotDexCortCtrl_intTest <- merged_pilotDexCortCtrl
merged_pilotDexCortCtrl_intTest$sample_id <- with(
  merged_pilotDexCortCtrl_intTest@meta.data,
  paste0(sublibrary, "_", bc1_well)
)

# 1) Split po sample_id
obj_list <- SplitObject(merged_pilotDexCortCtrl_intTest, split.by = "sample_id")

# (ważne) uniknij błędu future.globals.maxSize — na czas integracji działamy sekwencyjnie
plan(sequential)
options(future.globals.maxSize = 8 * 1024^3) 
plan(multisession, workers = 10)                

# 2) LogNormalize + HVG w każdej próbce
obj_list <- lapply(obj_list, function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = TRUE)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = TRUE)
  x
})

# 3) Cechy do integracji
features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 2000)

# 4) Anchory i integracja
anchors <- FindIntegrationAnchors(object.list = obj_list, anchor.features = features, verbose = TRUE)
obj_int <- IntegrateData(anchorset = anchors, verbose = TRUE)

# 5) Downstream na assay 'integrated'
DefaultAssay(obj_int) <- "integrated"

obj_int <- ScaleData(obj_int, verbose = TRUE)
obj_int <- RunPCA(obj_int, npcs = 30, verbose = TRUE)

obj_int <- FindNeighbors(obj_int, reduction = "pca", dims = 1:30, verbose = TRUE)
obj_int <- FindClusters(obj_int, resolution = 0.5, verbose = TRUE)

obj_int <- RunUMAP(obj_int, reduction = "pca", dims = 1:30, verbose = TRUE, seed.use = 777)

# 6) Podgląd
DimPlot(obj_int, reduction = "umap", group.by = "sample_id")
DimPlot(obj_int, reduction = "umap", group.by = "treatment")
DimPlot(obj_int, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
