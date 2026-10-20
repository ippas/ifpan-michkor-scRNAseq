# ##############################################################################
# ---- analysis for all aggregate protein coding ----
# ##############################################################################
exprMat_cluster_meanDataRNA <- prepare_proteinCodingGenesMatrix_fromSeurat_v2(
  integrated_pilotDexCortCtrl,
  assay           = "RNA",
  layer           = "data",
  biomart_version = "v112",
  method          = "mean",
  aggregate_cluster = TRUE,           # czy robić agregację po klastrach
  group.by = "seurat_clusters_res0.5"      
)

sctype_scores_clusterMeanDataRNA_scaled <- sctype_score(
  scRNAseqData           = exprMat_cluster_meanDataRNA,
  scaled                 = TRUE,
  gs                     = all_markers_list,   # alternatively: brain_markers_list
  gs2                    = NULL,               # negative markers (optional)
  gene_names_to_uppercase = TRUE
)

sctype_scores_clusterMeanDataRNA_noScaled <- sctype_score(
  scRNAseqData           = exprMat_cluster_meanDataRNA,
  scaled                 = FALSE,
  gs                     = all_markers_list,   # alternatively: brain_markers_list
  gs2                    = NULL,               # negative markers (optional)
  gene_names_to_uppercase = TRUE
)


sctype_scores_clusterMeanDataRNA_scaled_brain <- sctype_score(
  scRNAseqData           = exprMat_cluster_meanDataRNA,
  scaled                 = TRUE,
  gs                     = brain_markers_list,   # alternatively: brain_markers_list
  gs2                    = NULL,               # negative markers (optional)
  gene_names_to_uppercase = TRUE
)

sctype_scores_clusterMeanDataRNA_noScaled_brain <- sctype_score(
  scRNAseqData           = exprMat_cluster_meanDataRNA,
  scaled                 = FALSE,
  gs                     = brain_markers_list,   # alternatively: brain_markers_list
  gs2                    = NULL,               # negative markers (optional)
  gene_names_to_uppercase = TRUE
)

apply(sctype_scores_clusterMeanDataRNA, 2, function(x) {
  rownames(sctype_scores_clusterMeanDataRNA)[which.max(x)]
})


# ##############################################################################
# ---- analysis from integrate protein coding  ----
# ##############################################################################
exprMat_cluster_meanDataIntegrated <- prepare_proteinCodingGenesMatrix_fromSeurat_v2(
  integrated_pilotDexCortCtrl,
  assay = "integrated",
  layer = "data",
  biomart_version = "v112",
  method = "mean",
  aggregate_cluster = TRUE,
  group.by = "seurat_clusters_res0.5"
)

sctype_scores_clusterMeanDataIntegrated_scaled <- sctype_score(
  scRNAseqData           = exprMat_cluster_meanDataIntegrated,
  scaled                 = TRUE,
  gs                     = all_markers_list,   # alternatively: brain_markers_list
  gs2                    = NULL,               # negative markers (optional)
  gene_names_to_uppercase = TRUE
)

sctype_scores_clusterMeanDataIntegrated_noScale <- sctype_score(
  scRNAseqData           = exprMat_cluster_meanDataIntegrated,
  scaled                 = FALSE,
  gs                     = all_markers_list,   # alternatively: brain_markers_list
  gs2                    = NULL,               # negative markers (optional)
  gene_names_to_uppercase = TRUE
)

sctype_scores_clusterMeanDataIntegrated_scaled_brain <- sctype_score(
  scRNAseqData           = exprMat_cluster_meanDataIntegrated,
  scaled                 = TRUE,
  gs                     = brain_markers_list,   # alternatively: brain_markers_list
  gs2                    = NULL,               # negative markers (optional)
  gene_names_to_uppercase = TRUE
)

sctype_scores_clusterMeanDataIntegrated_noScaled_brain <- sctype_score(
  scRNAseqData           = exprMat_cluster_meanDataIntegrated,
  scaled                 = FALSE,
  gs                     = brain_markers_list,   # alternatively: brain_markers_list
  gs2                    = NULL,               # negative markers (optional)
  gene_names_to_uppercase = TRUE
)

apply(sctype_scores_clusterMeanDataIntegrated, 2, function(x) {
  rownames(sctype_scores_clusterMeanDataIntegrated)[which.max(x)]
})

sctype_scores_clusterMeanDataIntegrated %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "celltype") %>% 
  arrange(cluster_0)

# ##############################################################################
# ---- combine results to one table ----
# ##############################################################################
# --- helper do wyciągnięcia top celltype ---
get_top_celltype <- function(sctype_mat) {
  if (is.null(sctype_mat)) return(rep(NA_character_, length = ncol(sctype_mat)))
  apply(sctype_mat, 2, function(x) rownames(sctype_mat)[which.max(x)])
}

# --- łączymy wyniki w jeden tibble ---
df_top_sctype <- tibble(
  cluster                       = names(get_top_celltype(sctype_scores_clusterMeanDataRNA_scaled)),
  
  # RNA – all markers
  res0.5_RNA_scaled_allMarkers  = unname(get_top_celltype(sctype_scores_clusterMeanDataRNA_scaled)),
  res0.5_RNA_noScaled_allMarkers= unname(get_top_celltype(sctype_scores_clusterMeanDataRNA_noScaled)),
  
  # RNA – brain markers
  res0.5_RNA_scaled_brain       = unname(get_top_celltype(sctype_scores_clusterMeanDataRNA_scaled_brain)),
  res0.5_RNA_noScaled_brain     = unname(get_top_celltype(sctype_scores_clusterMeanDataRNA_noScaled_brain)),
  
  # Integrated – all markers
  res0.5_Int_scaled_allMarkers = unname(get_top_celltype(sctype_scores_clusterMeanDataIntegrated_scaled)),
  res0.5_Int_noScaled_allMarkers=unname(get_top_celltype(sctype_scores_clusterMeanDataIntegrated_noScale)),
  
  # Integrated – brain markers
  res0.5_Int_scaled_brain      = unname(get_top_celltype(sctype_scores_clusterMeanDataIntegrated_scaled_brain)),
  res0.5_Int_noScaled_brain    = unname(get_top_celltype(sctype_scores_clusterMeanDataIntegrated_noScaled_brain))
) %>%
  arrange(as.numeric(gsub("cluster_", "", cluster)))  # sortujemy wg numeru klastra

df_top_sctype
