

# seurat_featurePlot_binary_clusterSummary(
#   seurat_obj = integrated_pilotDexCortCtrl_filtered,
#   gene_pattern = "SGK1",
#   split.by = "treatment",
#   split_levels = c("CTRL", "CORT", "DEX"),
#   cluster_col = "res0.3_RNA_noScaled_brain",
#   return_data = F
# )
# 
# seurat_featurePlot_binary_clusterSummary(
#   seurat_obj = integrated_pilotDexCortCtrl_filtered,
#   gene_pattern = "SGK1",
#   split.by = "treatment",
#   split_levels = c("CTRL", "CORT", "DEX"),
#   cluster_col = "res0.3_RNA_noScaled_brain",
#   cluster_summary_figure = c("Astrocytes", "Oligodendrocytes", "Oligodendrocyte precursor cells", "Radial glial cells")
# )


# svg("results/tmp/AQP4_clusterSummary_example_binary.svg", width = 15, height = 6.5)
seurat_featurePlot_binary_clusterSummary(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  gene_pattern = "ZBTB16",
  split.by = "treatment",
  split_levels = c("CTRL", "CORT", "DEX"),
  cluster_col = "res0.3_RNA_noScaled_brain",
  label_size = 10,
  point_size = 0.6,
  detection_threshold = 0,
  cluster_summary_figure = c("Radial glial cells" = "RGCs",
                             "Astrocytes" = "Astrocytes", 
                             "Oligodendrocyte precursor cells" = "OPCs",
                             "Oligodendrocytes" = "Oligodendrocytes",
                             "Mature neurons" = "Mature neurons",
                             "Immature neurons" = "Immature neurons")
  )

# RBFOX3,MAP2,SYP,DLG4,TUBB3,MAPT,INA,GAP43,NRP1

# dev.off()

seurat_featurePlot_continuous(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  gene_pattern = "GRIK2",
  split.by = "treatment",
  split_levels = c("CTRL", "CORT", "DEX"),  # zachowane, ale kolejność i tak wymuszona
  cluster_col = "res0.3_RNA_noScaled_brain",
  # detection_threshold = 0,
  point_size = 0.4
)

# SLC6A1,GABBR1,GABBR2,GAD2,GAD1,SLC32A1

# 📁 Katalog docelowy
output_dir <- "/home/mateusz/projects/ifpan-michkor-scRNAseq/results/results_analysis_02.11.2025_v3_qcControl_res0.3_scTypeNoScaleBrain11/umap_binaryExpression_intersectGlobalNeuralGenes"

# 📜 Lista genów
genes <- c(
  "ADAMTS1", "ADARB1", "AFF1", "ARHGEF3", "ARRDC2", "BCL2L1", "BCL6",
  "CDKN1A", "CDO1", "CHKA", "CRISPLD2", "DCUN1D3", "DDIT4", "EGFR",
  "ERRFI1", "FAM107A", "FKBP5", "FOXO3", "HSD17B11", "IL6R", "KLF13",
  "KLF15", "KLF9", "LIMS2", "LPIN1", "MAP3K6", "MBNL2", "MERTK", "NEDD9",
  "NFKBIA", "PDK4", "PER1", "PIK3R1", "RASSF4", "RHOU", "RIN2", "SESN1",
  "SGK1", "SPSB1", "STIM1", "TSC22D3", "USP53", "ZBTB16", "ZFP36", "ZHX3"
)

# 🧠 Pętla po genach
for (gene in genes) {
  message(sprintf("▶️ Generating plot for gene: %s", gene))
  
  # 🔹 Ścieżka zapisu
  output_file <- file.path(output_dir, paste0(gene, "_cellSummary_binary.svg"))
  
  # 🔹 Utworzenie pliku SVG
  svg(filename = output_file, width = 15, height = 6.5)
  
  # 🔹 Wykres
  seurat_featurePlot_binary_clusterSummary(
    seurat_obj = integrated_pilotDexCortCtrl_filtered,
    gene_pattern = gene,
    split.by = "treatment",
    split_levels = c("CTRL", "CORT", "DEX"),
    cluster_col = "res0.3_RNA_noScaled_brain",
    label_size = 12,
    point_size = 0.6,
    cluster_summary_figure = c(
      "Astrocytes" = "Astrocytes",
      "Oligodendrocytes" = "Oligodendrocytes",
      "Oligodendrocyte precursor cells" = "OPCs",
      "Radial glial cells" = "RGCs"
    )
  ) %>% print()
  
  dev.off()
}

message("✅ All plots generated successfully.")



# ##############################################################################
# ---- global up 6 tissue derived ----
# ##############################################################################
genes <- c(
  "FKBP5","KLF9","NFKBIA","PDK4","SESN1","TSC22D3","TXNIP","ZBTB16","ZFP36","BCL6",
  "CD163","CEBPD","DDIT4","DUSP1","FOXO1","LPIN1","MT1M","NEDD9","RHOB","TIPARP",
  "ADAMTS9","AKAP12","CNKSR3","CRISPLD2","FOXO3","GADD45B","GLUL","IL6R","IRS2","KLF15",
  "LOX","MERTK","MT1X","MXD4","NET1","NFIL3","PER1","PIK3R1","RPS6KA2","SLC7A8",
  "TGFBR3","USP53","ABLIM3","ADAMTS1","ADARB1","AFF1","AHNAK","AKAP13","ALDH1A1","ANGPTL4",
  "ARID5B","BTG1","CDKN1A","CEBPB","COL4A1","DHCR24","DNAJB4","ERRFI1","F3","FAM107A",
  "FGD4","FSTL3","GADD45A","GPX3","HIF3A","IFNGR1","KLF6","MAOA","MBNL2","MCL1",
  "MT1E","MYC","NEXN","NFKBIZ","PTGS2","PTK2B","RGS2","S100A10","SAT1","SDC4",
  "SGK1","SLA","SLC46A3","SSH2","ST3GAL1","STK17B","TMEM100","TRIB1"
)

# 📁 Katalog docelowy
output_dir <- "/home/mateusz/projects/ifpan-michkor-scRNAseq/results/results_analysis_02.11.2025_v3_qcControl_res0.3_scTypeNoScaleBrain11/umap_binaryExpression_global6TissueDerived"


# 🧠 Pętla po genach
for (gene in genes) {
  message(sprintf("▶️ Generating plot for gene: %s", gene))
  
  # 🔹 Ścieżka zapisu
  output_file <- file.path(output_dir, paste0(gene, "_cellSummary_binary.svg"))
  
  # 🔹 Utworzenie pliku SVG
  svg(filename = output_file, width = 15, height = 6.5)
  
  # 🔹 Wykres
  seurat_featurePlot_binary_clusterSummary(
    seurat_obj = integrated_pilotDexCortCtrl_filtered,
    gene_pattern = gene,
    split.by = "treatment",
    split_levels = c("CTRL", "CORT", "DEX"),
    cluster_col = "res0.3_RNA_noScaled_brain",
    label_size = 12,
    point_size = 0.6,
    cluster_summary_figure = c(
      "Astrocytes" = "Astrocytes",
      "Oligodendrocytes" = "Oligodendrocytes",
      "Oligodendrocyte precursor cells" = "OPCs",
      "Radial glial cells" = "RGCs"
    )
  ) %>% print()
  
  dev.off()
}

message("✅ All plots generated successfully.")


# ##############################################################################
# ---- Nueral Minus Global ----
# ##############################################################################

genes <- c(
  "ACSL3","ADORA1","ADRA1D","AGTRAP","AKAP7","ARPC1A","ARRDC3","ATP10D","BCL6B","CCDC43",
  "CD63","CDC37L1","CNNM2","CNTFR","CORO2B","CRYAB","CSAD","CSNK1D","CSTF3","CYTIP",
  "DAGLB","DGAT2","DIO2","DTNBP1","DYNC1I1","FAM110B","FBN1","FOXJ3","FRMPD1","FZD1",
  "GADD45G","GHR","GOT1L1","GPD1","GPLD1","GPR155","GPR176","GPT2","HDC","HMGCS2",
  "HSPA1B","HSPB1","ID1","IFRD1","IGF2BP2","IL16","IPPK","IRF2","KCNG2","KCNK6",
  "KCTD5","KIF1B","LAMB1","LHFPL6","LPIN3","LRRC8A","LYPLA1","MANBA","MAPK4","MAST4",
  "MDH1","MGST1","MICAL2","MYO1C","NEDD4L","NKD2","NPLOC4","NPY1R","NTRK2","NUAK1",
  "NXN","PAWR","PELI1","PGM5","PLD3","PLXNA2","PNPLA2","PPP5C","PRDX6","PRKX",
  "PRODH","PRR5","PTPN14","RAB4A","RANBP9","RAPGEF3","RASL11B","RGS17","RILPL1","RIMS2",
  "SAP30","SCHIP1","SERPINA3","SH3RF1","SIPA1L2","SLC24A4","SLC25A33","SLC2A9","SLC38A2","SMPD2",
  "SNX25","SRXN1","SULF2","SVIL","SWAP70","TAPT1","TNS1","TOP1","TP63","TRAFD1",
  "TRIM2","TRIM36","TSC22D2","UXS1","ASPA","BCAT1","CABLES1","CHST1","CYP2D6","EHD3",
  "GJB6","HTRA1","IGSF11","LMOD1","MTMR2","NAPEPLD","NBL1","NMB","NT5E","NUDT9",
  "OTUD7B","PIM3","PLEKHF1","RBPMS","RDX","RNF19A","SLC25A13","SLC45A3","TBC1D9B","WIPF3",
  "WNT16","AZIN1","PEX11A","SMARCD2","THRA","TLE3","TOB2"
)

output_dir = "/home/mateusz/projects/ifpan-michkor-scRNAseq/results/results_analysis_02.11.2025_v3_qcControl_res0.3_scTypeNoScaleBrain11/umap_binaryExpression_NeuralMinusGlobal"

# 🧠 Pętla po genach
for (gene in genes) {
  message(sprintf("▶️ Generating plot for gene: %s", gene))
  
  # 🔹 Ścieżka zapisu
  output_file <- file.path(output_dir, paste0(gene, "_cellSummary_binary.svg"))
  
  # 🔹 Utworzenie pliku SVG
  svg(filename = output_file, width = 15, height = 6.5)
  
  # 🔹 Wykres
  seurat_featurePlot_binary_clusterSummary(
    seurat_obj = integrated_pilotDexCortCtrl_filtered,
    gene_pattern = gene,
    split.by = "treatment",
    split_levels = c("CTRL", "CORT", "DEX"),
    cluster_col = "res0.3_RNA_noScaled_brain",
    label_size = 12,
    point_size = 0.6,
    cluster_summary_figure = c(
      "Astrocytes" = "Astrocytes",
      "Oligodendrocytes" = "Oligodendrocytes",
      "Oligodendrocyte precursor cells" = "OPCs",
      "Radial glial cells" = "RGCs"
    )
  ) %>% print()
  
  dev.off()
}

message("✅ All plots generated successfully.")

