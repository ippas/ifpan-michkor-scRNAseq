# ##############################################################################
# ---- load packages and functions ----
# ##############################################################################

if (!requireNamespace("HGNChelper", quietly = TRUE)) {
  install.packages("HGNChelper")
}

require(HGNChelper)
install.packages("HGNChelper")
lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)

# źródła z GitHuba ScType
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")


# ##############################################################################
# ---- load data sctype ----
# ##############################################################################
sctype_db = readxl::read_excel("data/celltype_references/ScTypeDB/ScTypeDB_full.xlsx")

integrated_pilotDexCortCtrl <- JoinLayers(integrated_pilotDexCortCtrl, assay = "RNA")

Layers(integrated_pilotDexCortCtrl[["RNA"]])

# wybierz tylko Brain
brain_markers <- sctype_db %>%
  filter(tissueType %in% c("Brain", "Hippocampus"))

# zamień string markerów na listę wektorów
brain_markers_list <- lapply(brain_markers$geneSymbolmore1, function(x){
  str_split(x, pattern = ",")[[1]] %>% str_trim()
})
names(brain_markers_list) <- brain_markers$cellName


# use all rows of sctype_db (no filter by tissueType)
all_markers_list <- lapply(sctype_db$geneSymbolmore1, function(x) {
  stringr::str_split(x, pattern = ",")[[1]] %>% stringr::str_trim()
})

names(all_markers_list) <- sctype_db$cellName

# podgląd
brain_markers_list[1:3]

exprMat <- as.matrix(LayerData(integrated_pilotDexCortCtrl[["RNA"]], layer = "data"))

exprMat %>% .[1:10, 1:10]


rm(exprMat)

exprMat_prepared <- prepare_proteinCodingGenesMatrix_fromSeurat(
  integrated_pilotDexCortCtrl,
  assay = "RNA",
  layer = "data",
  biomart_version = "v112",
  method = "mean"
)


exprMat_prepared 


sctype_scores <- sctype_score(
  scRNAseqData = exprMat_prepared,
  scaled = TRUE,
  # gs = brain_markers_list,    # markery pozytywne
  gs = all_markers_list,
  gs2 = NULL,                 # możesz dorzucić negatywne, jeśli chcesz
  gene_names_to_uppercase = TRUE
)

sctype_scores %>% head %>% .[, 1:20]


# for each cell (column) find the row (cell type) with the highest score
best_type_per_cell <- apply(sctype_scores, 2, function(x) {
  rownames(sctype_scores)[which.max(x)]
})

summary_table <- table(best_type_per_cell)


