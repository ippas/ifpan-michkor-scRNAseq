BoxPlotNoZero <- function(seurat_obj, gene, group.by = "seurat_clusters", 
                          assay = "RNA", slot = "data", min_expr = 0) {
  # Pobierz ekspresję
  expr <- GetAssayData(seurat_obj, assay = assay, slot = slot)[gene, ]
  
  # Zbuduj dataframe
  df <- data.frame(
    expr  = expr,
    group = seurat_obj@meta.data[[group.by]]
  )
  
  # Usuń zera
  df <- df[df$expr > min_expr, ]
  
  # Jeśli nic nie zostało, zwróć ostrzeżenie
  if (nrow(df) == 0) {
    warning(paste("No values >", min_expr, "for gene:", gene))
    return(NULL)
  }
  
  # Boxplot
  ggplot(df, aes(x = group, y = expr, fill = group)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
    theme_classic() +
    labs(
      title = gene,
      x = group.by,
      y = paste0("Expression (", assay, "/", slot, ")")
    ) +
    theme(legend.position = "none")
}

BoxPlotNoZero(seurat_obj = integrated_pilotDexCortCtrl, 
              gene = "ENSG00000160307-S100B",
              group.by = "seurat_clusters_res0.13") 


FeaturePlot(integrated_pilotDexCortCtrl, features = "ENSG00000160307-S100B")

integrated_pilotDexCortCtrl@meta.data %>% head

VlnPlot(
  integrated_pilotDexCortCtrl,
  features = "ENSG00000160307-S100B",       # np. "FKBP5"
  group.by = "seurat_clusters_res0.13", # grupowanie po klastrach
  pt.size = 0.1                 # wielkość kropek (0 = ukryj)
)

rownames(integrated_pilotDexCortCtrl[["RNA"]]) %>% 
  as.data.frame() %>% 
  set_colnames("name_id") %>% 
  filter(grepl("S100B", name_id))
