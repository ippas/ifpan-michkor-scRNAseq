genes_global

genes_neural

genes_neuralDown <- c(
  "CXXC5","HBEGF","PLK2","CCL2","CDKN1C","ENC1","PLAU","BIRC5","BUB1B","COL1A1",
  "DUSP10","DUSP14","EMP1","FGF2","FRMD6","FZD2","GINS2","IGFBP5","LRRC17","MCM10",
  "PER2","PER3","PIM1","PRSS23","PSAT1","SOX4","SYTL2","TCF4","TGFB3","TK1",
  "TNFSF15","VCAM1","ADAM12","ANLN","ARL4C","ATP1B1","ATP8B4","AURKB","BBC3","BCL3",
  "BHLHE40","CBR3","CCND2","CD83","CDCA5","CDCA8","CGNL1","CKAP2L","CLDN1","CNTNAP1",
  "CTNNAL1","CTSC","CXCL12","DBP","DIRAS2","DUSP6","E2F2","EGR3","EMP3","EPHB2",
  "G0S2","GAMT","GDF15","GINS4","GRAMD1B","HES1","HK2","HMGA1","HMOX1","ID2",
  "ID3","IER3","IGFBP3","IL6","INHBA","IRF1","ITGA11","LMNB1","LOXL1","LRIG1",
  "LXN","MAD2L1","MARCKS","MND1","MPP2","NFIX","NOG","NR3C1","NR4A3","PCNA",
  "PDLIM7","PHLDA1","PIF1","PLEKHG3","PLK4","PMEPA1","PTPRE","RAD51AP1","RELB","RGS4",
  "SALL2","SCARA3","SEMA3D","SEMA7A","SERTAD4","SLC6A6","SLC8A1","SMOX","SOCS3","SPON2",
  "TEF","TGFB2","TGFBI","TM4SF1","TMEM200A","TNC","TRIB3","UNC5B","WWC1"
)

genes_GlobalDown <- genes <- c(
  "CXXC5","HBEGF","PLK2","CCL2","CDKN1C","ENC1","PLAU","BIRC5","BUB1B","COL1A1",
  "DUSP10","DUSP14","EMP1","FGF2","FRMD6","FZD2","GINS2","IGFBP5","LRRC17","MCM10",
  "PER2","PER3","PIM1","PRSS23","PSAT1","SOX4","SYTL2","TCF4","TGFB3","TK1",
  "TNFSF15","VCAM1","ADAM12","ANLN","ARL4C","ATP1B1","ATP8B4","AURKB","BBC3","BCL3",
  "BHLHE40","CBR3","CCND2","CD83","CDCA5","CDCA8","CGNL1","CKAP2L","CLDN1","CNTNAP1",
  "CTNNAL1","CTSC","CXCL12","DBP","DIRAS2","DUSP6","E2F2","EGR3","EMP3","EPHB2",
  "G0S2","GAMT","GDF15","GINS4","GRAMD1B","HES1","HK2","HMGA1","HMOX1","ID2",
  "ID3","IER3","IGFBP3","IL6","INHBA","IRF1","ITGA11","LMNB1","LOXL1","LRIG1",
  "LXN","MAD2L1","MARCKS","MND1","MPP2","NFIX","NOG","NR3C1","NR4A3","PCNA",
  "PDLIM7","PHLDA1","PIF1","PLEKHG3","PLK4","PMEPA1","PTPRE","RAD51AP1","RELB","RGS4",
  "SALL2","SCARA3","SEMA3D","SEMA7A","SERTAD4","SLC6A6","SLC8A1","SMOX","SOCS3","SPON2",
  "TEF","TGFB2","TGFBI","TM4SF1","TMEM200A","TNC","TRIB3","UNC5B","WWC1"
)


# 🔹 Pomocnicza funkcja do pobrania listy genów
get_cluster_genes <- function(cluster_name, results_list) {
  dex <- results_list[[cluster_name]]$DEX_vs_CTRL %>%
    filter(FDR < 0.01, abs(logFC) > 1)
  
  cort <- results_list[[cluster_name]]$CORT_vs_CTRL %>%
    filter(FDR < 0.01, abs(logFC) > 1)
  
  unique(c(dex$gene, cort$gene))
}

# 🔹 Zbuduj listę trzech wektorów genów
gr_dependent_genes_list <- list(
  mature_neurons_genes = get_cluster_genes("cluster_Mature neurons", edgeR_anovaOneFactor_treatment_pseudobulk_res0.4_clusters16_celltypes11_results),
  astrocytes_genes = get_cluster_genes("cluster_Astrocytes", edgeR_anovaOneFactor_treatment_pseudobulk_res0.4_clusters16_celltypes11_results),
  oligodendrocytes_genes = get_cluster_genes("cluster_Oligodendrocytes", edgeR_anovaOneFactor_treatment_pseudobulk_res0.4_clusters16_celltypes11_results)
)



intersect(c(genes_neural, genes_global), gr_dependent_genes_list$oligodendrocytes_genes)



intersect(c(genes_GlobalDown, genes_neuralDown), gr_dependent_genes_list$mature_neurons_genes)
