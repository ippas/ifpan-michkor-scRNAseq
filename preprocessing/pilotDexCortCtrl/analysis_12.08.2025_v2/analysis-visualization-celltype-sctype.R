###############################################################################
# Load ScType database
###############################################################################
sctype_db <- readxl::read_excel(
  "data/celltype_references/ScTypeDB/ScTypeDB_full.xlsx"
)

###############################################################################
# Prepare Seurat object: join layers & check RNA layers
###############################################################################
integrated_pilotDexCortCtrl <- JoinLayers(
  integrated_pilotDexCortCtrl,
  assay = "RNA"
)

Layers(integrated_pilotDexCortCtrl[["RNA"]])

###############################################################################
# Prepare marker gene lists
###############################################################################

## --- Brain-specific markers
brain_markers <- sctype_db %>%
  dplyr::filter(tissueType == "Brain")

brain_markers_list <- lapply(
  brain_markers$geneSymbolmore1,
  function(x) stringr::str_split(x, ",")[[1]] %>% stringr::str_trim()
)
names(brain_markers_list) <- brain_markers$cellName

## --- All markers (no filter)
all_markers_list <- lapply(
  sctype_db$geneSymbolmore1,
  function(x) stringr::str_split(x, ",")[[1]] %>% stringr::str_trim()
)
names(all_markers_list) <- sctype_db$cellName

# Quick preview
brain_markers_list[1:3]

###############################################################################
# Prepare expression matrix
###############################################################################
exprMat_prepared <- prepare_proteinCodingGenesMatrix_fromSeurat(
  integrated_pilotDexCortCtrl,
  assay           = "integrated",
  layer           = "data",
  biomart_version = "v112",
  method          = "mean"
)

###############################################################################
# Run ScType scoring
###############################################################################
sctype_scores <- sctype_score(
  scRNAseqData           = exprMat_prepared,
  scaled                 = T,
  gs                     = all_markers_list,   # alternatively: brain_markers_list
  gs2                    = NULL,               # negative markers (optional)
  gene_names_to_uppercase = TRUE
)

# Preview scores
sctype_scores %>%
  head() %>%
  .[, 1:20]

###############################################################################
# Assign best cell type per cell
###############################################################################
best_type_per_cell <- apply(sctype_scores, 2, function(x) {
  rownames(sctype_scores)[which.max(x)]
})

summary_table <- table(best_type_per_cell)
summary_table %>% 
  as.data.frame %>% 
  arrange(desc(Freq))

# best_type_per_cell <- best_type_per_cell %>%
#   tolower() %>%                # zamiana na małe litery
#   gsub(" ", "_", .) 


data.frame(
  cell_id        = names(best_type_per_cell),
  celltype_raw   = best_type_per_cell,
  celltype_clean = best_type_per_cell %>%
    tolower() %>%
    gsub(" ", "_", .),
  stringsAsFactors = FALSE
) %>% 
  left_join(., {sctype_db %>% select(tissueType, cellName)}, by = c("celltype_raw" = "cellName")) %>% 
  group_by(celltype_clean) %>% 
  nest() %>% 
  mutate(n_cells = map(data, ~ .x %>% nrow)) %>% 
  unnest(n_cells) %>% 
  mutate(celltype_filter10 = ifelse(n_cells > 10, celltype_clean, "below_10_cells")) %>% 
  mutate(celltype_filter50 = ifelse(n_cells > 50, celltype_clean, "below_50_cells")) %>% 
  mutate(celltype_filter100 = ifelse(n_cells > 100, celltype_clean, "below_100_cells")) %>% 
  mutate(celltype_filter200 = ifelse(n_cells > 200, celltype_clean, "below_200_cells")) %>% 
  unnest() %>% 
  mutate(neural_cells = ifelse(tissueType == "Brain" | tissueType == "Hippocampus", celltype_clean, "no_neural")) -> best_type_per_cell_df

  
best_type_per_cell_df <- best_type_per_cell_df %>%
  dplyr::select(cell_id,
                sctypeAll_raw      = celltype_raw,
                sctypeAll_clean    = celltype_clean,
                sctypeAll_tissue   = tissueType,
                sctypeAll_n_cells  = n_cells,
                sctypeAll_min10    = celltype_filter10,
                sctypeAll_min50    = celltype_filter50,
                sctypeAll_min100   = celltype_filter100,
                sctypeAll_min200   = celltype_filter200,
                sctypeAll_neural   = neural_cells)

meta_to_add <- best_type_per_cell_df %>%
  dplyr::filter(cell_id %in% Cells(integrated_pilotDexCortCtrl)) %>%
  dplyr::arrange(match(cell_id, Cells(integrated_pilotDexCortCtrl))) %>%
  dplyr::distinct(cell_id, .keep_all = TRUE) %>%   # <-- to ważne!
  dplyr::select(-cell_id) %>%
  as.data.frame()

# ustaw rownames
rownames(meta_to_add) <- Cells(integrated_pilotDexCortCtrl)

# dodaj do Seurat
integrated_pilotDexCortCtrl <- AddMetaData(
  integrated_pilotDexCortCtrl,
  metadata = meta_to_add
)

# ##############################################################################
# ---- add column to metadata ----
# ##############################################################################
integrated_pilotDexCortCtrl <- AddMetaData(
  object   = integrated_pilotDexCortCtrl,
  metadata = best_type_per_cell[Cells(integrated_pilotDexCortCtrl)],
  col.name = "celltype_sctypeAll_v1"
)

calc_entropy_counts <- function(counts) {
  # counts = wektor liczności (np. c(50, 30, 20))
  p <- counts / sum(counts)
  entropy <- -sum(p * log2(p))
  return(entropy)
}


integrated_pilotDexCortCtrl@meta.data %>% 
  filter(celltype_trailmaker == "missing") %>% 
  .$celltype_sctypeAll_v1 %>% 
  table %>% 
  as.data.frame() %>% 
  set_colnames(c("celltype", "freq")) %>% 
  arrange(desc(freq)) %>% head(20)


plot_umap_cluster_faceted(
  integrated_pilotDexCortCtrl,
  cluster_var = "sctypeAll_neural",
  # plot_title = "UMAP Projection Across All Samplest",
  facet_var = NULL,
  # facet_order = c("CTRL", "DEX", "CORT"),
  label_clusters = TRUE,               # czy pokazywać etykiety klastrów
  text_size_axis_title = 18,           # rozmiar tytułów osi
  text_size_axis_text = 14,             # rozmiar tekstu na osiach (tick labels)
  text_size_legend_title = 16,         # rozmiar tytułu legendy
  text_size_legend_text = 16,          # rozmiar tekstu legendy
  text_size_facet = 14,                # rozmiar tekstu podtytułów paneli (facetów)
  text_size_plot_title = 18 ,           # rozmiar tytułu wykresu
  text_size_label = 6,
  n_row_legend = 2,
  filter_exclude = c("no_neural"),
  # filter_include = c("oligodendrocytes"),
  filter_var = "sctypeAll_neural"
)

