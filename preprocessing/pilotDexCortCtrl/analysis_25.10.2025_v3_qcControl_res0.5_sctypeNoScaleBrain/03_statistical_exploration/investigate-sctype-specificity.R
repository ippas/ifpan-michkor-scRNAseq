# --- przygotowanie macierzy ekspresji dla klastrów (średnia) ---
exprMat_cluster_meanData <- prepare_proteinCodingGenesMatrix_fromSeurat_v2(
  seurat_obj       = integrated_pilotDexCortCtrl,
  assay            = "RNA",
  layer            = "data",
  biomart_version  = "v112",
  method           = "mean",
  aggregate_cluster = TRUE,
  group.by         = "seurat_clusters_res0.5"
)

# --- obliczenie macierzy wyników ScType dla zestawu "brain markers" ---
sctype_noscaled_brain <- sctype_score(
  exprMat_cluster_meanData,
  scaled_data = FALSE,     # <── noScaled
  gs = brain_markers_list, # <── brain markers
  gene_names_to_uppercase = TRUE
)


sctype_noscaled_brain[]

apply(sctype_noscaled_brain, 2, function(x) {
  rownames(sctype_noscaled_brain)[which.max(x)]
})


integrated_pilotDexCortCtrl_filtered@meta.data %>% names
