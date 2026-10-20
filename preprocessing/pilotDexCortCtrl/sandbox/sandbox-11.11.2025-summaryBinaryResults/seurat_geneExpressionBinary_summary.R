seurat_binaryDetection_summaryTable <- function(
    seurat_obj,
    genes_vector,
    assay = "RNA",
    slot = "counts",
    split.by = "treatment",
    split_levels = c("CTRL", "CORT", "DEX"),
    detection_threshold = 0,
    verbose = TRUE
) {
  requireNamespace("dplyr")
  
  # 🔹 Przygotowanie wynikowej listy
  results_list <- list()
  
  # 🔹 Iteracja po genach
  for (gene_pattern in genes_vector) {
    if (verbose) message(sprintf("▶️ Processing gene: %s", gene_pattern))
    
    all_features <- rownames(GetAssayData(seurat_obj, assay = assay, slot = slot))
    matched <- grep(paste0("-", gene_pattern, "$"), all_features, value = TRUE)
    if (length(matched) == 0) {
      warning(sprintf("❌ Gene not found: %s", gene_pattern))
      next
    }
    if (length(matched) > 1) {
      message("⚠️ Multiple matches found for ", gene_pattern, ", using first: ", matched[1])
      matched <- matched[1]
    }
    
    # 🔹 Ekspresja i detekcja
    expr_values <- GetAssayData(seurat_obj, assay = assay, slot = slot)[matched, ]
    detected_binary <- as.numeric(expr_values > detection_threshold)
    
    # 🔹 Dodaj metadane z grupami
    meta_df <- seurat_obj@meta.data
    if (!split.by %in% colnames(meta_df))
      stop(sprintf("❌ Column '%s' not found in metadata", split.by))
    
    df <- data.frame(
      gene_symbol = gene_pattern,
      treatment = meta_df[[split.by]],
      detected = detected_binary
    )
    
    # 🔹 Agregacja po treatment
    summary_df <- df |>
      dplyr::group_by(treatment) |>
      dplyr::summarise(
        detected_count = sum(detected, na.rm = TRUE),
        total_cells = dplyr::n(),
        percent = round(100 * detected_count / total_cells, 1),
        .groups = "drop"
      )
    
    # 🔹 Uporządkuj kolejność grup
    summary_df$treatment <- factor(summary_df$treatment, levels = split_levels)
    summary_df <- summary_df |> dplyr::arrange(treatment)
    
    # 🔹 Rozszerz w jedną linię
    row_out <- data.frame(gene_symbol = gene_pattern)
    for (grp in split_levels) {
      if (grp %in% summary_df$treatment) {
        row_vals <- summary_df |> dplyr::filter(treatment == grp)
        row_out[[paste0(grp, "_detected")]] <- row_vals$detected_count
        row_out[[paste0(grp, "_total")]] <- row_vals$total_cells
        row_out[[paste0(grp, "_percent")]] <- row_vals$percent
      } else {
        row_out[[paste0(grp, "_detected")]] <- 0
        row_out[[paste0(grp, "_total")]] <- 0
        row_out[[paste0(grp, "_percent")]] <- 0
      }
    }
    
    results_list[[gene_pattern]] <- row_out
  }
  
  # 🔹 Połącz w jedną tabelę
  result_df <- dplyr::bind_rows(results_list)
  
  if (verbose) message("✅ Summary table generated successfully.")
  
  return(result_df)
}
seurat_binaryDetection_summaryByCluster_splitList <- function(
    seurat_obj,
    genes_vector,
    clusters_to_check,
    assay = "RNA",
    slot = "counts",
    split.by = "treatment",
    split_levels = c("CTRL", "CORT", "DEX"),
    cluster_col = "res0.3_RNA_noScaled_brain",
    detection_threshold = 0,
    verbose = TRUE
) {
  requireNamespace("dplyr")
  
  meta_df <- seurat_obj@meta.data
  if (!cluster_col %in% colnames(meta_df))
    stop(sprintf("❌ Cluster column '%s' not found in metadata.", cluster_col))
  if (!split.by %in% colnames(meta_df))
    stop(sprintf("❌ Split column '%s' not found in metadata.", split.by))
  
  # Upewnij się, że assay i slot są ustawione spójnie
  if (!is.null(assay)) Seurat::DefaultAssay(seurat_obj) <- assay
  
  results_by_cluster <- list()
  n_genes <- length(genes_vector)
  
  # Prelista wszystkich feature'ów, żeby szybko dopasowywać nazwę
  all_features <- rownames(Seurat::GetAssayData(seurat_obj, assay = assay, slot = slot))
  
  for (cluster_name in clusters_to_check) {
    if (verbose) message(sprintf("\n🔹 Cluster: %s", cluster_name))
    
    cluster_cells <- rownames(meta_df[meta_df[[cluster_col]] == cluster_name, ])
    if (length(cluster_cells) == 0) {
      warning(sprintf("⚠️ No cells found for cluster: %s", cluster_name))
      next
    }
    
    # treatment dla komórek w danym klastrze (w stałej kolejności)
    treatments_vec <- meta_df[cluster_cells, split.by, drop = TRUE]
    
    cluster_rows <- list()
    
    for (i in seq_along(genes_vector)) {
      gene_pattern <- genes_vector[i]
      percent_done <- round(100 * i / n_genes, 1)
      if (verbose) message(sprintf("   ▶️ Gene: %s (%d / %d, %s%%)", gene_pattern, i, n_genes, percent_done))
      
      # dopasowanie nazwy feature na końcu (twoja konwencja "-GENE$")
      matched <- grep(paste0("-", gene_pattern, "$"), all_features, value = TRUE)
      if (length(matched) == 0) {
        warning(sprintf("❌ Gene not found: %s", gene_pattern))
        next
      }
      if (length(matched) > 1) {
        if (verbose) message("   ⚠️ Multiple matches, using first: ", matched[1])
        matched <- matched[1]
      }
      
      # KLUCZOWA ZMIANA: FetchData gwarantuje gęsty wektor w kolejności komórek
      df_gene <- Seurat::FetchData(
        object = seurat_obj,
        vars = matched,
        cells = cluster_cells,
        slot = slot
      )
      # wektor ekspresji zgodny z kolejnością 'cluster_cells'
      expr_vec <- df_gene[[matched]]
      detected_binary <- as.integer(expr_vec > detection_threshold)
      
      df <- data.frame(
        gene_symbol = gene_pattern,
        treatment = treatments_vec,
        detected = detected_binary,
        stringsAsFactors = FALSE
      )
      
      summary_df <- df |>
        dplyr::group_by(treatment) |>
        dplyr::summarise(
          detected_count = sum(detected, na.rm = TRUE),
          total_cells    = dplyr::n(),
          percent        = ifelse(total_cells > 0,
                                  round(100 * detected_count / total_cells, 1),
                                  0),
          .groups = "drop"
        ) |>
        dplyr::mutate(treatment = factor(treatment, levels = split_levels))
      
      # szeroki format: CTRL/CORT/DEX x (detected/total/percent)
      wide <- data.frame(gene_symbol = gene_pattern, stringsAsFactors = FALSE)
      for (grp in split_levels) {
        if (grp %in% summary_df$treatment) {
          row_vals <- summary_df[summary_df$treatment == grp, , drop = FALSE]
          wide[[paste0(grp, "_detected")]] <- row_vals$detected_count
          wide[[paste0(grp, "_total")]]    <- row_vals$total_cells
          wide[[paste0(grp, "_percent")]]  <- row_vals$percent
        } else {
          wide[[paste0(grp, "_detected")]] <- 0L
          wide[[paste0(grp, "_total")]]    <- 0L
          wide[[paste0(grp, "_percent")]]  <- 0
        }
      }
      
      cluster_rows[[gene_pattern]] <- wide
    }
    
    results_by_cluster[[cluster_name]] <- dplyr::bind_rows(cluster_rows)
  }
  
  if (verbose) message("\n✅ List of data frames (per cluster) generated (fixed detection).")
  return(results_by_cluster)
}

genes_global <- c(
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


genes_neural <- c(
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


binary_summary_df <- seurat_binaryDetection_summaryTable(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  genes_vector = c(genes_global, genes_neural),
  split.by = "treatment",
  split_levels = c("CTRL", "CORT", "DEX")
)

print(binary_summary_df)


binary_summary_df %>% 
  mutate(diff_percent = DEX_percent - CTRL_percent) %>% 
  mutate(log2ratio = log2(DEX_percent/CTRL_percent + 0.000001)) %>% 
  filter(log2ratio > 0) %>% 
  filter(DEX_percent> 5) %>% arrange(desc(diff_percent))



# ##############################################################################
genes_to_check <- c("FKBP5", "SGK1", "ZBTB16")
clusters_to_check <- c("Immature neurons", "Mature neurons", "Astrocytes", "Oligodendrocytes")
clusters_to_check <- c("Mature neurons", "Astrocytes", "Oligodendrocytes")

# cluster_summary_list <- seurat_binaryDetection_summaryByCluster_splitList(
#   seurat_obj = integrated_pilotDexCortCtrl_filtered,
#   genes_vector = genes_to_check,
#   clusters_to_check = clusters_to_check,
#   cluster_col = "res0.3_RNA_noScaled_brain"
# )

cluster_summary_list_neural <- seurat_binaryDetection_summaryByCluster_splitList(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  genes_vector =  c(genes_neural),
  clusters_to_check = clusters_to_check,
  cluster_col = "res0.3_RNA_noScaled_brain",
  detection_threshold = 0  # np. można dać 0.1 jeśli chcesz próg
)

cluster_summary_list_neural$`Mature neurons` %>% 
  mutate(diff_percent = DEX_percent - CTRL_percent) %>% 
  mutate(log2ratio = log2(DEX_percent/CTRL_percent + 0.000001)) %>% 
  filter(log2ratio > 0) %>% 
  filter(DEX_percent> 5) %>% arrange(desc(diff_percent)) %>% 
  filter(diff_percent > 5)


cluster_summary_list_neural$Astrocytes %>% 
  mutate(diff_percent = DEX_percent - CTRL_percent) %>% 
  mutate(log2ratio = log2(DEX_percent/CTRL_percent + 0.000001)) %>% 
  filter(log2ratio > 0) %>% 
  filter(DEX_percent> 5) %>% arrange(desc(diff_percent))


cluster_summary_list_neural$Oligodendrocytes %>% 
  mutate(diff_percent = DEX_percent - CTRL_percent) %>% 
  mutate(log2ratio = log2(DEX_percent/CTRL_percent + 0.000001)) %>% 
  filter(log2ratio > 0) %>% 
  filter(DEX_percent> 5) %>% arrange(desc(diff_percent))






cluster_summary_list_GlobalNeural <- seurat_binaryDetection_summaryByCluster_splitList(
  seurat_obj = integrated_pilotDexCortCtrl_filtered,
  genes_vector =  c(genes_global, genes_neural),
  clusters_to_check = clusters_to_check,
  cluster_col = "res0.3_RNA_noScaled_brain",
  detection_threshold = 0  # np. można dać 0.1 jeśli chcesz próg
)



cluster_summary_list_GlobalNeural$`Mature neurons` %>% 
  mutate(diff_percent = DEX_percent - CTRL_percent) %>% 
  mutate(log2ratio = log2(DEX_percent/CTRL_percent + 0.000001)) %>% 
  filter(log2ratio > 0) %>% 
  filter(DEX_percent> 5) %>% arrange(desc(diff_percent)) %>% 
  filter(diff_percent > 5) %>% .$diff_percent %>% sd

cluster_summary_list_GlobalNeural$`Mature neurons` %>% 
  mutate(diff_percent = CORT_percent - CTRL_percent) %>% 
  mutate(log2ratio = log2(CORT_percent/CTRL_percent + 0.000001)) %>% 
  filter(log2ratio > 0) %>% 
  filter(CORT_percent > 5) %>% arrange(desc(diff_percent)) %>% 
  filter(diff_percent > 5) %>% head

cluster_summary_list_GlobalNeural$Astrocytes %>% 
  mutate(diff_percent = CORT_percent - CTRL_percent) %>% 
  mutate(log2ratio = log2(CORT_percent/CTRL_percent + 0.000001)) %>% 
  filter(log2ratio > 0) %>% 
  filter(CORT_percent > 5) %>% arrange(desc(diff_percent)) %>% 
  filter(diff_percent > 5) %>% head


cluster_summary_list_GlobalNeural$Astrocytes %>% 
  mutate(diff_percent = DEX_percent - CTRL_percent) %>% 
  mutate(log2ratio = log2(DEX_percent/CTRL_percent + 0.000001)) %>% 
  filter(log2ratio > 0) %>% 
  filter(DEX_percent> 5) %>% arrange(desc(diff_percent)) %>% 
  arrange(desc(diff_percent)) %>% 
  filter(diff_percent > 5) %>% .$diff_percent %>% sd


cluster_summary_list_GlobalNeural$`Oligodendrocytes` %>% 
  mutate(diff_percent = CORT_percent - CTRL_percent) %>% 
  mutate(log2ratio = log2(CORT_percent/CTRL_percent + 0.000001)) %>% 
  filter(log2ratio > 0) %>% 
  filter(CORT_percent > 5) %>% arrange(desc(diff_percent)) %>% 
  filter(diff_percent > 5) %>% .$diff_percent %>% sd


cluster_summary_list_GlobalNeural$Oligodendrocytes %>% 
  mutate(diff_percent = DEX_percent - CTRL_percent) %>% 
  mutate(log2ratio = log2(DEX_percent/CTRL_percent + 0.000001)) %>% 
  filter(log2ratio > 0) %>% 
  filter(DEX_percent> 5) %>% arrange(desc(diff_percent)) %>% 
  filter(diff_percent > 5) %>% .$diff_percent %>% sd






# cluster_summary_list_neural <- seurat_binaryDetection_summaryByCluster_splitList(
#   seurat_obj = integrated_pilotDexCortCtrl_filtered,
#   genes_vector =  c(protein_coding_genes),
#   clusters_to_check = clusters_to_check,
#   cluster_col = "res0.3_RNA_noScaled_brain",
#   detection_threshold = 0  # np. można dać 0.1 jeśli chcesz próg
# )
