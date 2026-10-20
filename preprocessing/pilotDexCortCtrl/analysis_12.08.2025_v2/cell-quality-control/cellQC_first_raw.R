# ##############################################################################
# ---- cehck data ----
# ##############################################################################
integrated_pilotDexCortCtrl@assays %>% names

slotNames(integrated_pilotDexCortCtrl@assays$integrated)
slotNames(integrated_pilotDexCortCtrl@assays$RNA)


qc_filter_report <- function(
    seurat_obj,
    mt_gene_list,
    max_percent_mt = 10,
    min_counts = 500,
    max_counts = 40000,
    min_features = 200,
    max_features = 6000
){
  # pobierz dane
  counts <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
  
  # policz procent MT
  mt_genes <- intersect(rownames(counts), mt_gene_list)
  percent_mt <- (colSums(counts[mt_genes, , drop = FALSE]) / colSums(counts)) * 100
  
  meta <- seurat_obj@meta.data %>%
    rownames_to_column("cell_id") %>%
    mutate(percent_MTDNA = percent_mt)
  
  # liczba komórek na starcie
  n_start <- nrow(meta)
  
  # kolejno filtrujemy i liczymy ile odpadło
  step1 <- meta %>% filter(percent_MTDNA < max_percent_mt)
  step2 <- step1 %>% filter(nCount_RNA > min_counts)
  step3 <- step2 %>% filter(nCount_RNA < max_counts)
  step4 <- step3 %>% filter(nFeature_RNA > min_features)
  step5 <- step4 %>% filter(nFeature_RNA < max_features)
  
  # przygotowujemy raport
  report <- tibble(
    step = c("start", 
             paste0("percent.mt <", max_percent_mt),
             paste0("nCount_RNA >", min_counts),
             paste0("nCount_RNA <", max_counts),
             paste0("nFeature_RNA >", min_features),
             paste0("nFeature_RNA <", max_features),
             "final"),
    n_cells = c(n_start,
                nrow(step1),
                nrow(step2),
                nrow(step3),
                nrow(step4),
                nrow(step5),
                nrow(step5)),
    removed = c(0,
                n_start - nrow(step1),
                nrow(step1) - nrow(step2),
                nrow(step2) - nrow(step3),
                nrow(step3) - nrow(step4),
                nrow(step4) - nrow(step5),
                n_start - nrow(step5))
  )
  
  return(list(
    report = report,
    keep_cells = step5$cell_id
  ))
}


# 1. usuwamy komórki z bardzo niską liczbą odczytów (puste krople)
filtered <- subset(
  integrated_pilotDexCortCtrl,
  subset = nCount_RNA > 500
)

# 2. usuwamy komórki z ekstremalnie wysoką liczbą odczytów (doublety)
filtered <- subset(
  filtered,
  subset = nCount_RNA < 40000
)


subset(
  integrated_pilotDexCortCtrl,
  subset = nFeature_RNA > 200 & nFeature_RNA < 5000
) %>% ncol



tmp <- download_biomart_genes(version = "v112")

tmp %>% 
  .$gene_biotype %>% unique()

tmp %>% 
  filter(gene_biotype %in% c("Mt_tRNA", "Mt_rRNA"))

integrated_pilotDexCortCtrl[["RNA"]] %>% rownames %>% 
  as.data.frame() %>% 
  set_colnames("name_id") %>% 
  filter(base::grepl("-MT-", name_id )) %>% .$name_id -> mitochondrial_name_id


rna_counts <- GetAssayData(integrated_pilotDexCortCtrl, assay = "RNA", slot = "counts")

# Wybieramy tylko te geny, które mamy w naszej liście
mt_genes <- intersect(rownames(rna_counts), mitochondrial_name_id)

# sumujemy odczyty mitochondrialne per komórka
mt_counts_per_cell <- colSums(rna_counts[mt_genes, , drop = FALSE])

# sumujemy całkowite odczyty per komórka
total_counts_per_cell <- colSums(rna_counts)

# liczymy procent
percent_mt <- (mt_counts_per_cell / total_counts_per_cell) * 100

percent_mt %>% 
  as.data.frame() %>%  
  set_colnames("percent_MTDNA") %>% 
  filter(percent_MTDNA > 5) %>% dim


# ##############################################################################
rna_counts <- GetAssayData(integrated_pilotDexCortCtrl, assay = "RNA", slot = "counts")

# 1. wybieramy mitochondrialne geny
mt_genes <- intersect(rownames(rna_counts), mitochondrial_name_id)

# 2. liczymy procent MT dla każdej komórki
percent_mt <- (colSums(rna_counts[mt_genes, , drop = FALSE]) /
                 colSums(rna_counts)) * 100

# 3. tworzymy ramkę danych z metadanymi + percent.mt
qc_df <- integrated_pilotDexCortCtrl@meta.data %>%
  rownames_to_column("cell_id") %>%
  mutate(percent_MTDNA = percent_mt)

# 4. stosujemy wszystkie filtry jednocześnie
qc_filtered <- qc_df %>%
  filter(
    percent_MTDNA < 10,        # np. mniej niż 10% genów mitochondrialnych
    nCount_RNA > 1000,          # min. liczba transkryptów
    nCount_RNA < 20000,        # max. liczba transkryptów (usuwa doublety)
    nFeature_RNA > 200,        # min. liczba genów
    nFeature_RNA < 6000        # max. liczba genów
  )

# 5. sprawdzamy, ile komórek zostało
cat("Pozostało", nrow(qc_filtered), "komórek z", nrow(qc_df),
    "(", round(nrow(qc_filtered)/nrow(qc_df)*100, 2), "%)\n")


qc_summary_df <- integrated_pilotDexCortCtrl@meta.data %>%
  rownames_to_column("cell_id") %>%
  mutate(percent_MTDNA = percent_mt) %>%
  mutate(
    # --- dla progu 500 ---
    fail_reason_500 = paste(
      ifelse(percent_MTDNA > 10, "high_mt", ""),
      ifelse(nCount_RNA <= 500, "low_counts", ""),
      ifelse(nCount_RNA >= 10000, "high_counts", ""),
      ifelse(nFeature_RNA <= 200, "low_genes", ""),
      ifelse(nFeature_RNA >= 4000, "high_genes", ""),
      sep = ";"
    ),
    # usuwamy puste wpisy (same "; ")
    qc_pass_500 = ifelse(fail_reason_500 == ";;;;", T, F),
    
    # --- dla progu 1000 ---
    fail_reason_1000 = paste(
      ifelse(percent_MTDNA > 10, "high_mt", ""),
      ifelse(nCount_RNA <= 1000, "low_counts", ""),
      ifelse(nCount_RNA >= 20000, "high_counts", ""),
      ifelse(nFeature_RNA <= 200, "low_genes", ""),
      ifelse(nFeature_RNA >= 6000, "high_genes", ""),
      sep = ";"
    ),
    # fail_reason_1000 = ifelse(gsub(";", "", fail_reason_1000) == "", NA, fail_reason_1000),
    qc_pass_1000 = ifelse(fail_reason_1000 == ";;;;", T, F)
  )


# podsumowania:
cat("Dla progu 500: ", sum(qc_summary_df$qc_pass_500), "/", nrow(qc_summary_df),
    " (", round(mean(qc_summary_df$qc_pass_500)*100,2), "%)\n")

cat("Dla progu 1000:", sum(qc_summary_df$qc_pass_1000), "/", nrow(qc_summary_df),
    " (", round(mean(qc_summary_df$qc_pass_1000)*100,2), "%)\n")


# save(integrated_pilotDexCortCtrl, file = "data/pilotDexCortCtrl/r_outputs/integrated_pilotDexCortCtrl_11.09.2025_v3.RData")



qc_fail_summary_500 <- qc_summary_df %>%
  # rozdziel powody odrzucenia na osobne kolumny
  separate_rows(fail_reason_500, sep = ";") %>%
  # usuń puste wpisy
  filter(fail_reason_500 != "") %>%
  # zlicz liczbę wystąpień każdego powodu
  count(fail_reason_500, name = "n_cells_removed") %>%
  # posortuj malejąco
  arrange(desc(n_cells_removed)) %>%
  # policz procent komórek z całej populacji
  mutate(percent_removed = round(n_cells_removed / nrow(qc_summary_df) * 100, 2))

qc_fail_summary_500
