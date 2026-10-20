seurat_featurePlot_binary_clusterSummary_generateMultiSVG <- function(
    seurat_obj,
    genes_vector,
    output_dir,
    split.by = "treatment",
    split_levels = c("CTRL", "CORT", "DEX"),
    cluster_col = "res0.3_RNA_noScaled_brain",
    cluster_summary_figure = c(
      "Astrocytes" = "Astrocytes",
      "Oligodendrocytes" = "Oligodendrocytes",
      "Oligodendrocyte precursor cells" = "OPCs",
      "Radial glial cells" = "RGCs"
    ),
    force_create_dir = TRUE,
    verbose = TRUE
) {
  # 🧩 Require needed package
  if (!requireNamespace("progress", quietly = TRUE)) {
    install.packages("progress")
  }
  library(progress)
  
  # 📁 Create directory if needed
  if (!dir.exists(output_dir)) {
    if (force_create_dir) {
      dir.create(output_dir, recursive = TRUE)
      if (verbose) message(sprintf("📂 Created output directory: %s", output_dir))
    } else {
      stop(sprintf("❌ Directory does not exist: %s", output_dir))
    }
  }
  
  # 📊 Progress bar (only if not verbose)
  if (!verbose) {
    pb <- progress_bar$new(
      format = "Generating plots [:bar] :percent | ETA: :eta",
      total = length(genes_vector),
      width = 70
    )
  }
  
  # 🔁 Iterate over genes
  for (gene in genes_vector) {
    if (verbose) {
      message(sprintf("▶️ Generating plot for gene: %s", gene))
    } else {
      pb$tick()
    }
    
    output_file <- file.path(output_dir, paste0(gene, "_cellSummary_binary.svg"))
    svg(filename = output_file, width = 15, height = 6.5)
    
    tryCatch({
      seurat_featurePlot_binary_clusterSummary(
        seurat_obj = seurat_obj,
        gene_pattern = gene,
        split.by = split.by,
        split_levels = split_levels,
        cluster_col = cluster_col,
        label_size = 12,
        point_size = 0.6,
        cluster_summary_figure = cluster_summary_figure
      ) %>% print()
    }, error = function(e) {
      message(sprintf("⚠️ Error while plotting %s: %s", gene, e$message))
    })
    
    dev.off()
  }
  
  if (verbose) {
    message("✅ All SVG plots generated successfully.")
  } else {
    pb$terminate()
    message("\n✅ All SVG plots generated successfully.")
  }
}


neuronal_markers <- c(
  "RBFOX3","SNAP25","STMN2","SYT1","TUBB3","MAP2","NEFL","NEFM","GAP43","CPLX2","DLG4","GRIN1",
  
  # 🔹 Glutamatergic (pobudzające neurony)
  "SLC17A7","SLC17A6","SLC17A8","CAMK2A","GRIA2","GRIN2B","SATB2","TBR1",
  
  # 🔹 GABAergic (hamujące neurony)
  "GAD1","GAD2","SLC32A1","SLC6A1","PVALB","SST","VIP","RELN",
  
  # 🔹 Cortical layers and projection markers
  "CUX1","CUX2","RGS4","RORB","BCL11B","TLE4","FOXP2","TBR1",
  
  # 🔹 Dopaminergic neurons
  "TH","DDC","SLC6A3","SLC18A2",
  
  # 🔹 Cholinergic neurons
  "CHAT","SLC5A7","SLC18A3",
  
  # 🔹 Noradrenergic neurons
  "DBH","SLC6A2","TH",
  
  # 🔹 Serotonergic neurons
  "TPH2","SLC6A4","SLC18A2","DDC",
  
  # 🔹 Additional abundant neuronal genes (często używane w scRNA-seq)
  "SNCB","STX1B","SLC12A5","NCAM1","RCAN2","PCLO","SYN1"
)

seurat_featurePlot_binary_clusterSummary_generateMultiSVG(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  genes_vector = neuronal_markers,
  output_dir = "/home/mateusz/projects/ifpan-michkor-scRNAseq/results/results_analysis_02.11.2025_v3_qcControl_res0.3_scTypeNoScaleBrain11/umap_binaryExpression_NeuronalMarkersChatGPT",
  force_create_dir = TRUE,
  verbose = FALSE,
  cluster_summary_figure = c("Astrocytes" = "Astrocytes", "Oligodendrocytes" = "Oligodendrocytes",
                             "Oligodendrocyte precursor cells" = "OPCs", "Radial glial cells" = "RGCs",
                             "Immature neurons" = "Immature neurons",
                             "Mature neurons" = "Mature neurons")
)




neuronal_markers <- cluster_list$`Mature neurons` %>%
  separate(gene, into = c("ensembl_id", "gene_symbol"), sep = "-", extra = "merge") %>% .$gene_symbol


seurat_featurePlot_binary_clusterSummary_generateMultiSVG(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  genes_vector = neuronal_markers,
  output_dir = "/home/mateusz/projects/ifpan-michkor-scRNAseq/results/results_analysis_02.11.2025_v3_qcControl_res0.3_scTypeNoScaleBrain11/umap_binaryExpression_MatureNeruronsFindMarkerSeurat",
  force_create_dir = TRUE,
  verbose = FALSE,
  cluster_summary_figure = c("Astrocytes" = "Astrocytes", "Oligodendrocytes" = "Oligodendrocytes",
                             "Oligodendrocyte precursor cells" = "OPCs", "Radial glial cells" = "RGCs",
                             "Immature neurons" = "Immature neurons",
                             "Mature neurons" = "Mature neurons")
)


# ##############################################################################
# ---- immature neurons ----
# ##############################################################################

neuronal_markers <- cluster_list$`Immature neurons` %>%
  separate(gene, into = c("ensembl_id", "gene_symbol"), sep = "-", extra = "merge") %>% .$gene_symbol


seurat_featurePlot_binary_clusterSummary_generateMultiSVG(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  genes_vector = neuronal_markers,
  output_dir = "/home/mateusz/projects/ifpan-michkor-scRNAseq/results/results_analysis_02.11.2025_v3_qcControl_res0.3_scTypeNoScaleBrain11/umap_binaryExpression_ImmatureNeruronsFindMarkerSeurat",
  force_create_dir = TRUE,
  verbose = FALSE,
  cluster_summary_figure = c("Astrocytes" = "Astrocytes", "Oligodendrocytes" = "Oligodendrocytes",
                             "Oligodendrocyte precursor cells" = "OPCs", "Radial glial cells" = "RGCs",
                             "Immature neurons" = "Immature neurons",
                             "Mature neurons" = "Mature neurons")
)



# ##############################################################################
# ---- oligodendrocytes ----
# ##############################################################################
oligodendrocytes_markers <- cluster_list$Oligodendrocytes %>%
  separate(gene, into = c("ensembl_id", "gene_symbol"), sep = "-", extra = "merge") %>% .$gene_symbol


seurat_featurePlot_binary_clusterSummary_generateMultiSVG(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  genes_vector = oligodendrocytes_markers,
  output_dir = "/home/mateusz/projects/ifpan-michkor-scRNAseq/results/results_analysis_02.11.2025_v3_qcControl_res0.3_scTypeNoScaleBrain11/umap_binaryExpression_OligodendrocytesFindMarkerSeurat",
  force_create_dir = TRUE,
  verbose = FALSE,
  cluster_summary_figure = c("Astrocytes" = "Astrocytes", "Oligodendrocytes" = "Oligodendrocytes",
                             "Oligodendrocyte precursor cells" = "OPCs", "Radial glial cells" = "RGCs",
                             "Immature neurons" = "Immature neurons",
                             "Mature neurons" = "Mature neurons")
)


# ##############################################################################
# ---- OPC ----
# ##############################################################################
OPC_markers <- cluster_list$`Oligodendrocyte precursor cells` %>%
  separate(gene, into = c("ensembl_id", "gene_symbol"), sep = "-", extra = "merge") %>% .$gene_symbol


seurat_featurePlot_binary_clusterSummary_generateMultiSVG(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  genes_vector = OPC_markers,
  output_dir = "/home/mateusz/projects/ifpan-michkor-scRNAseq/results/results_analysis_02.11.2025_v3_qcControl_res0.3_scTypeNoScaleBrain11/umap_binaryExpression_OPCFindMarkerSeurat",
  force_create_dir = TRUE,
  verbose = FALSE,
  cluster_summary_figure = c("Astrocytes" = "Astrocytes", "Oligodendrocytes" = "Oligodendrocytes",
                             "Oligodendrocyte precursor cells" = "OPCs", "Radial glial cells" = "RGCs",
                             "Immature neurons" = "Immature neurons",
                             "Mature neurons" = "Mature neurons")
)


# ##############################################################################
# ---- Astrocytes ----
# ##############################################################################
astrocytes_markers <- cluster_list$Astrocytes %>%
  separate(gene, into = c("ensembl_id", "gene_symbol"), sep = "-", extra = "merge") %>% .$gene_symbol


seurat_featurePlot_binary_clusterSummary_generateMultiSVG(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  genes_vector = astrocytes_markers,
  output_dir = "/home/mateusz/projects/ifpan-michkor-scRNAseq/results/results_analysis_02.11.2025_v3_qcControl_res0.3_scTypeNoScaleBrain11/umap_binaryExpression_AstrocytesFindMarkerSeurat",
  force_create_dir = TRUE,
  verbose = FALSE,
  cluster_summary_figure = c("Astrocytes" = "Astrocytes", "Oligodendrocytes" = "Oligodendrocytes",
                             "Oligodendrocyte precursor cells" = "OPCs", "Radial glial cells" = "RGCs",
                             "Immature neurons" = "Immature neurons",
                             "Mature neurons" = "Mature neurons")
)

# ##############################################################################
# ---- Radial Glial Celss ----
# ##############################################################################
RGC_markers <- cluster_list$`Radial glial cells` %>%
  separate(gene, into = c("ensembl_id", "gene_symbol"), sep = "-", extra = "merge") %>% .$gene_symbol


seurat_featurePlot_binary_clusterSummary_generateMultiSVG(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  genes_vector = astrocytes_markers,
  output_dir = "/home/mateusz/projects/ifpan-michkor-scRNAseq/results/results_analysis_02.11.2025_v3_qcControl_res0.3_scTypeNoScaleBrain11/umap_binaryExpression_RGCFindMarkerSeurat",
  force_create_dir = TRUE,
  verbose = FALSE,
  cluster_summary_figure = c("Astrocytes" = "Astrocytes", "Oligodendrocytes" = "Oligodendrocytes",
                             "Oligodendrocyte precursor cells" = "OPCs", "Radial glial cells" = "RGCs",
                             "Immature neurons" = "Immature neurons",
                             "Mature neurons" = "Mature neurons")
)


# ##############################################################################
# ---- OPC markers GPT ----
# ##############################################################################
opc_markers <- c(
  # 🔹 Klasyczne, wysoko swoiste markery OPC
  "PDGFRA", "CSPG4", "SOX10", "OLIG1", "OLIG2",
  
  # 🔹 Wczesne i proliferujące OPC
  "PTPRZ1", "GPR17", "VCAN", "NKX2-2", "TNR", "TCF7L2", "SOX8",
  
  # 🔹 Przejście OPC → pre-oligodendrocyty
  "QKI", "ENPP6", "TMEM108", "LHFPL3", "PLPPR1", "PID1", "BCAS1",
  
  # 🔹 Regulacja i odpowiedź metaboliczna (często GR-zależne)
  "PDK4", "DDIT4", "SGK1", "SIRT2", "ERBB3", "TMEM178B",
  
  # 🔹 Dodatkowe markery często wykrywane w klastrach OPC w scRNA-seq
  "GPR37", "GPR37L1", "PDE7B", "DUSP15", "MAGI1", "CLDN11", "HOPX",
  "CXCR4", "C1QL1", "EGFR", "TTYH1", "ANO4", "SOX17"
)

seurat_featurePlot_binary_clusterSummary_generateMultiSVG(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  genes_vector = opc_markers,
  output_dir = "/home/mateusz/projects/ifpan-michkor-scRNAseq/results/results_analysis_02.11.2025_v3_qcControl_res0.3_scTypeNoScaleBrain11/umap_binaryExpression_OPCMarkersGPT",
  force_create_dir = TRUE,
  verbose = FALSE,
  cluster_summary_figure = c("Astrocytes" = "Astrocytes", "Oligodendrocytes" = "Oligodendrocytes",
                             "Oligodendrocyte precursor cells" = "OPCs", "Radial glial cells" = "RGCs",
                             "Immature neurons" = "Immature neurons",
                             "Mature neurons" = "Mature neurons")
)


# ##############################################################################
# ---- grGlobal UP ----
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
  "SGK1","SLA","SLC46A3","SSH2","ST3GAL1","STK17B","TMEM100","TRIB1","ABLIM2","ADM",
  "ADRB2","AGFG2","ALCAM","ALOX5AP","ALPP","AMPD3","ANKRD1","ANXA4","AOX1","APOD",
  "ARHGEF3","ARL4D","ARRDC2","BAIAP2","BCL2L1","BCR","BTG2","C1QTNF1","CBLB","CCND3",
  "CD59","CDO1","CHKA","CITED2","CKB","CLMN","CP","CPM","CTDSPL","CYFIP2",
  "CYP4B1","DCUN1D3","DPP4","DUSP10","DUSP6","EGFR","ELMO1","FAM13A","FBXO32","FERMT2",
  "FOSL2","FZD4","GABARAPL1","GPRC5B","HIPK2","HSD17B11","IL6","IRAK2","ITGA10","JAK1",
  "KDR","KLF13","KLF4","LCN2","LEF1","LIMS2","LRRFIP1","LRRK2","MAN1A1","MAP3K6",
  "MEGF9","MGLL","MKNK2","MLXIP","MT1A","MT2A","MTHFD2","MYLK","NDRG1","NRBP2",
  "NRCAM","P2RX5","PAG1","PAPPA","PCBP3","PER2","PKP2","PLIN2","PRKAG2","PTGDS",
  "PTPRM","RAB31","RASSF4","RBM47","RHOU","RIN2","RIN3","S100P","S1PR1","SAA1",
  "SCNN1A","SLC31A2","SLC39A14","SLC44A2","SLC7A2","SNX10","SOD2","SORBS1","SORT1","SPOCK1",
  "SPSB1","SRGN","STAB1","STIM1","SYBU","TGM2","TIMP4","TLE1","TNS3","TRIM26",
  "TSC22D1","UAP1","VCL","WIPI1","XDH","ZFAND5","ZHX3","ZNF189"
)

seurat_featurePlot_binary_clusterSummary_generateMultiSVG(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  genes_vector = genes,
  output_dir = "/home/mateusz/projects/ifpan-michkor-scRNAseq/results/results_analysis_02.11.2025_v3_qcControl_res0.3_scTypeNoScaleBrain11/umap_binaryExpression_GlobalUp5TissuesDerived",
  force_create_dir = TRUE,
  verbose = FALSE,
  cluster_summary_figure = c("Astrocytes" = "Astrocytes", "Oligodendrocytes" = "Oligodendrocytes",
                             "Oligodendrocyte precursor cells" = "OPCs", "Radial glial cells" = "RGCs",
                             "Immature neurons" = "Immature neurons",
                             "Mature neurons" = "Mature neurons")
)


# ##############################################################################
# ---- interesect Neural Global ----
# ##############################################################################
genes <- c(
  "ADAMTS1","ADARB1","AFF1","ARHGEF3","ARRDC2","BCL2L1","BCL6","CDKN1A","CDO1","CHKA",
  "CRISPLD2","DCUN1D3","DDIT4","EGFR","ERRFI1","FAM107A","FKBP5","FOXO3","HSD17B11","IL6R",
  "KLF13","KLF15","KLF9","LIMS2","LPIN1","MAP3K6","MBNL2","MERTK","NEDD9","NFKBIA",
  "PDK4","PER1","PIK3R1","RASSF4","RHOU","RIN2","SESN1","SGK1","SPSB1","STIM1",
  "TSC22D3","USP53","ZBTB16","ZFP36","ZHX3"
)


seurat_featurePlot_binary_clusterSummary_generateMultiSVG(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  genes_vector = genes,
  output_dir = "/home/mateusz/projects/ifpan-michkor-scRNAseq/results/results_analysis_02.11.2025_v3_qcControl_res0.3_scTypeNoScaleBrain11/umap_binaryExpression_IntersectGlobalNeuralTissueCell",
  force_create_dir = TRUE,
  verbose = FALSE,
  cluster_summary_figure = c("Astrocytes" = "Astrocytes", "Oligodendrocytes" = "Oligodendrocytes",
                             "Oligodendrocyte precursor cells" = "OPCs", "Radial glial cells" = "RGCs",
                             "Immature neurons" = "Immature neurons",
                             "Mature neurons" = "Mature neurons")
)


# ##############################################################################
# ---- Neural minus global ----
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

seurat_featurePlot_binary_clusterSummary_generateMultiSVG(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  genes_vector = genes,
  output_dir = "/home/mateusz/projects/ifpan-michkor-scRNAseq/results/results_analysis_02.11.2025_v3_qcControl_res0.3_scTypeNoScaleBrain11/umap_binaryExpression_NeuralMinusGlobal",
  force_create_dir = TRUE,
  verbose = FALSE,
  cluster_summary_figure = c("Astrocytes" = "Astrocytes", "Oligodendrocytes" = "Oligodendrocytes",
                             "Oligodendrocyte precursor cells" = "OPCs", "Radial glial cells" = "RGCs",
                             "Immature neurons" = "Immature neurons",
                             "Mature neurons" = "Mature neurons")
)


# ##############################################################################
# ---- astrocytes markers GPT ----
# ##############################################################################
astro_markers <- c(
  "GFAP",         # glial fibrillary acidic protein — klasyczny marker astrocytów :contentReference[oaicite:0]{index=0},
  "S100B",        # S100β — białko wapniowowe często używane jako marker astrocytów :contentReference[oaicite:1]{index=1},
  "ALDH1L1",      # dehydrogenaza/enzym glialny, szeroko używany marker astrocytów :contentReference[oaicite:2]{index=2},
  "AQP4",         # aquaporin-4, woda kanał w astrocytach, np. na wypustkach przy naczyniach :contentReference[oaicite:3]{index=3},
  "GLUL",         # glutamine synthetase, astrocyty odpowiadają za metabolizm glutaminy / glutaminianu :contentReference[oaicite:4]{index=4},
  "SLC1A2",       # transporter glutaminianu / EAAT2 (GLT-1), astrocytarne wychwytywanie glutaminianu :contentReference[oaicite:5]{index=5},
  "SLC1A3",       # EAAT1 / GLAST transporter glutaminianu w astrocytach :contentReference[oaicite:6]{index=6},
  "GJA1",         # connexin 43 (gap junction), astrocyty łączą się w sieć komunikacji przez kanały połączeń :contentReference[oaicite:7]{index=7},
  "NDRG2",        # marker astrocytów obecnie uważany za stosunkowo specyficzny :contentReference[oaicite:8]{index=8},
  "SOX9"         # czynnik transkrypcyjny związany z glialnym fenotypem astrocytów :contentReference[oaicite:9]{index=9}
)

seurat_featurePlot_binary_clusterSummary_generateMultiSVG(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  genes_vector = astro_markers,
  output_dir = "/home/mateusz/projects/ifpan-michkor-scRNAseq/results/results_analysis_02.11.2025_v3_qcControl_res0.3_scTypeNoScaleBrain11/umap_binaryExpression_AstrocytesMarkersGPT",
  force_create_dir = TRUE,
  verbose = FALSE,
  cluster_summary_figure = c("Astrocytes" = "Astrocytes", "Oligodendrocytes" = "Oligodendrocytes",
                             "Oligodendrocyte precursor cells" = "OPCs", "Radial glial cells" = "RGCs",
                             "Immature neurons" = "Immature neurons",
                             "Mature neurons" = "Mature neurons")
)

# ##############################################################################
# ---- supplementary ----
# ##############################################################################
genes <- c("FKBP5", "ZBTB16", "SESN1",
           "TOB2", "TLE3", "AZIN1",
           "AQP4", "EPHA6", "SLC1A2",
           "MAPK4", "KLF9", "GLUL", 
           "SIRT2", "MBP", "BCAS1",
           "SGK1", "CNNM2", "KLF13",
           "CACNA1B", "MYT1L", "RBFOX3")



seurat_featurePlot_binary_clusterSummary_generateMultiSVG(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  genes_vector = genes,
  output_dir = "/home/mateusz/projects/ifpan-michkor-scRNAseq/results/results_analysis_02.11.2025_v3_qcControl_res0.3_scTypeNoScaleBrain11/umap_binaryExpression_SupplementaryGenes",
  force_create_dir = TRUE,
  verbose = FALSE,
  cluster_summary_figure = c("Astrocytes" = "Astrocytes", 
                             "Radial glial cells" = "RGCs",
                             "Mature neurons" = "Mature neurons",
                             "Immature neurons" = "Immature neurons",
                             "Oligodendrocytes" = "Oligodendrocytes",
                             "Oligodendrocyte precursor cells" = "OPCs"
                             )
)



# ##############################################################################
# ---- supplementary Interesection grSingature + edgeR ----
# ##############################################################################
genes <- c("LHFPL6", "MAST4", "IGSF11",
           "FKBP5", "ZBTB16", "FOXO1",
           "ADAMTS9", "MERTK", "PIK3R1",
           "FBXO32", "TNS3")




seurat_featurePlot_binary_clusterSummary_generateMultiSVG(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  genes_vector = genes,
  output_dir = "/home/mateusz/projects/ifpan-michkor-scRNAseq/results/results_analysis_02.11.2025_v3_qcControl_res0.3_scTypeNoScaleBrain11/umap_binaryExpression_grSignautresEdgeR",
  force_create_dir = TRUE,
  verbose = FALSE,
  cluster_summary_figure = c(
    "Radial glial cells" = "RGCs",
    "Astrocytes" = "Astrocytes",
    "Immature neurons" = "Immature neurons",
    "Mature neurons" = "Mature neurons",
    "Oligodendrocyte precursor cells" = "OPCs",
    "Oligodendrocytes" = "Oligodendrocytes")
)


# ##############################################################################
# ---- supplementary supplementary v2 ----
# ##############################################################################
genes <- c("FKBP5", "ZBTB16", "SESN1",
           "TOB2", "TLE3", "AZIN1",
           "AQP4", "EPHA6", "SLC1A2",
           "MAPK4", "KLF9", "GLUL",
           "SIRT2", "MBP", "BCAS1",
           "SGK1", "CNNM2", "KLF13", 
           "CACNA1B", "MYT1L", "RBFOX3", 
           "PAPPA", "MAST4",
           "NR3C1", 
           "NR3C2")

seurat_featurePlot_binary_clusterSummary_generateMultiSVG(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  genes_vector = genes,
  output_dir = "/home/mateusz/projects/ifpan-michkor-scRNAseq/results/results_analysis_02.11.2025_v3_qcControl_res0.3_scTypeNoScaleBrain11/umap_binaryExpression_supplementaryV2_6cellTypes",
  force_create_dir = TRUE,
  verbose = FALSE,
  cluster_summary_figure = c(
    "Radial glial cells" = "RGCs",
    "Astrocytes" = "Astrocytes",
    "Immature neurons" = "Immature neurons",
    "Mature neurons" = "Mature neurons",
    "Oligodendrocyte precursor cells" = "OPCs",
    "Oligodendrocytes" = "Oligodendrocytes")
)



seurat_featurePlot_binary_clusterSummary_generateMultiSVG(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  genes_vector = genes,
  output_dir = "/home/mateusz/projects/ifpan-michkor-scRNAseq/results/results_analysis_02.11.2025_v3_qcControl_res0.3_scTypeNoScaleBrain11/umap_binaryExpression_supplementaryV2_3cellTypes",
  force_create_dir = TRUE,
  verbose = FALSE,
  cluster_summary_figure = c(
    "Astrocytes" = "Astrocytes",
    "Mature neurons" = "Mature neurons",
    "Oligodendrocytes" = "Oligodendrocytes")
)


# ##############################################################################
# ---- supplementary supplementary v2 ----
# ##############################################################################
genes <- c("IL6R", "SCHIP1", "CBLB", "NFKBIZ")



seurat_featurePlot_binary_clusterSummary_generateMultiSVG(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  genes_vector = genes,
  output_dir = "/home/mateusz/projects/ifpan-michkor-scRNAseq/results/results_analysis_02.11.2025_v3_qcControl_res0.3_scTypeNoScaleBrain11/umap_binaryExpression_supplementaryV3_3cellTypes",
  force_create_dir = TRUE,
  verbose = FALSE,
  cluster_summary_figure = c(
    "Astrocytes" = "Astrocytes",
    "Mature neurons" = "Mature neurons",
    "Oligodendrocytes" = "Oligodendrocytes")
)



seurat_featurePlot_binary_clusterSummary_generateMultiSVG(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  genes_vector = c("NR3C1", "NR3C2"),
  output_dir = "/home/mateusz/projects/ifpan-michkor-scRNAseq/results/figures/umapGeneExpression_11.12.2025/",
  force_create_dir = TRUE,
  verbose = FALSE,
  cluster_summary_figure = c(
    "Radial glial cells"              = "RGCs",
    "Astrocytes"                      = "Astrocytes", 
    "Oligodendrocytes"                = "Oligodendrocytes",
    "Oligodendrocyte precursor cells" = "OPCs", 
    "Immature neurons"                = "Immature neurons",
    "Mature neurons"                  = "Mature neurons")
)


seurat_featurePlot_binary_clusterSummary(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  gene_pattern = "NR3C1",
  split.by = "treatment",
  split_levels = c("CTRL", "CORT", "DEX"),
  cluster_col = "res0.3_RNA_noScaled_brain",
  label_size = 12,
  point_size = 0.6,
  cluster_summary_figure = c(
    "Radial glial cells" = "RGCs",
    "Astrocytes" = "Astrocytes",
    "Oligodendrocyte precursor cells" = "OPCs",
    "Oligodendrocytes" = "Oligodendrocytes",
    "Immature neurons" = "Immature neurons",
    "Mature neurons" = "Mature neurons"
  )
)
