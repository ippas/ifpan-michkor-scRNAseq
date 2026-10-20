library(dplyr)
library(purrr)
library(readr)  # jeśli chcesz zapisywać pliki .tsv

# 🔹 załóżmy, że masz tabelę resolution_df jak powyżej
resolution_df <- tibble::tribble(
  ~resolution, ~n_clusters,
  "seurat_clusters_res0.01", 2,
  "seurat_clusters_res0.02", 3,
  "seurat_clusters_res0.021", 3,
  "seurat_clusters_res0.022", 3,
  "seurat_clusters_res0.025", 3,
  "seurat_clusters_res0.03", 4,
  "seurat_clusters_res0.04", 5,
  "seurat_clusters_res0.05", 6,
  "seurat_clusters_res0.06", 7,
  "seurat_clusters_res0.07", 7,
  "seurat_clusters_res0.08", 8,
  "seurat_clusters_res0.09", 8,
  "seurat_clusters_res0.1", 9,
  "seurat_clusters_res0.13", 10,
  "seurat_clusters_res0.15", 11,
  "seurat_clusters_res0.17", 13,
  "seurat_clusters_res0.2", 15,
  "seurat_clusters_res0.3", 16,
  "seurat_clusters_res0.4", 17,
  "seurat_clusters_res0.5", 20,
  "seurat_clusters_res0.8", 23,
  "seurat_clusters_res1", 32
)

# 🔹 wyciągamy tylko wartość liczbową po 'res'
resolutions <- gsub("seurat_clusters_res", "", resolution_df$resolution)

# 🔹 uruchamiamy funkcję combine_sctype_results() dla każdej rozdzielczości
sctype_results_list <- map(
  resolutions,
  ~ combine_sctype_results(
    seurat_obj         = integrated_pilotDexCortCtrl_filtered,
    res                = .x,
    all_markers_list   = all_markers_list,
    brain_markers_list = brain_markers_list,
    assay              = "RNA"
  )
)

# 🔹 nadaj nazwy elementom listy
names(sctype_results_list) <- paste0("res", resolutions)

# ✅ Wynik: lista ramek danych
sctype_results_list



sctype_summary_df <- imap_dfr(
  sctype_results_list,
  ~ {
    df <- .x
    res <- .y
    
    tibble(
      resolution = res,
      n_clusters = nrow(df),
      n_scaled_allMarkers   = n_distinct(df[[str_subset(names(df), "_scaled_allMarkers")]]),
      n_noScaled_allMarkers = n_distinct(df[[str_subset(names(df), "_noScaled_allMarkers")]]),
      n_scaled_brain        = n_distinct(df[[str_subset(names(df), "_scaled_brain")]]),
      n_noScaled_brain      = n_distinct(df[[str_subset(names(df), "_noScaled_brain")]])
    )
  }
)

# 🔹 posortuj po wartości numerycznej resolution
sctype_summary_df <- sctype_summary_df %>%
  mutate(res_numeric = as.numeric(str_replace(resolution, "res", ""))) %>%
  arrange(res_numeric) %>%
  select(-res_numeric)

sctype_summary_df %>% 
  as.data.frame()


sctype_results_list$res0.3 %>% 
  select(cluster, res0.3_RNA_noScaled_allMarkers, res0.3_RNA_noScaled_brain) %>% 
  set_colnames(c("cluster", "clustersAnnotationAllMarkers", "clustersAnnotationBrainMarkers")) %>% 
  mutate(SameCellType = ifelse(clustersAnnotationAllMarkers == clustersAnnotationBrainMarkers, TRUE, FALSE))


sctype_results_list[c("res0.3", "res0.13", "res0.5")] %>% 
  lapply(., function(x){
    x %>% select(-cluster_id)
  })

library(dplyr)
library(tibble)
library(purrr)

library(dplyr)
library(tibble)

for (res_name in c("res0.13", "res0.3", "res0.5")) {
  message("\n🔗 Przetwarzam ", res_name, " ...")
  
  # --- pobierz ramkę z wynikami ScType ---
  df_res <- sctype_results_list[[res_name]] %>%
    mutate(cluster_id = gsub("cluster_", "", cluster))
  
  join_by_df <- "cluster_id"
  join_col   <- paste0("seurat_clusters_", res_name)
  
  # --- sprawdź, czy kolumna istnieje w metadata ---
  if (!join_col %in% colnames(integrated_pilotDexCortCtrl_filtered@meta.data)) {
    warning("⚠️ Brak kolumny ", join_col, " w Seurat metadata, pomijam ", res_name)
    next
  }
  
  # --- dopasowanie typów danych ---
  df_res[[join_by_df]] <- as.character(df_res[[join_by_df]])
  integrated_pilotDexCortCtrl_filtered@meta.data[[join_col]] <-
    as.character(integrated_pilotDexCortCtrl_filtered@meta.data[[join_col]])
  
  # --- przygotuj dane do joinowania ---
  meta_tmp <- integrated_pilotDexCortCtrl_filtered@meta.data %>%
    rownames_to_column("cell_id")
  
  # dodaj tymczasową kolumnę z numerem klastra
  meta_tmp$tmp_cluster_id <- meta_tmp[[join_col]]
  
  meta_add <- meta_tmp %>%
    left_join(df_res, by = c("tmp_cluster_id" = "cluster_id")) %>%
    arrange(match(cell_id, rownames(integrated_pilotDexCortCtrl_filtered@meta.data))) %>%
    select(all_of(setdiff(colnames(df_res), c("cluster", "cluster_id")))) %>%
    as.data.frame()
  
  # --- dodaj kolumny do obiektu Seurat ---
  integrated_pilotDexCortCtrl_filtered <- AddMetaData(
    integrated_pilotDexCortCtrl_filtered,
    metadata = meta_add
  )
  
  # --- raport dopasowania ---
  n_matched <- sum(!is.na(meta_add[, 1]))
  n_total   <- nrow(meta_add)
  message("✅ ", res_name, ": dopasowano ", n_matched, "/", n_total, " komórek (", round(100 * n_matched / n_total, 1), "%)")
}



# ##############################################################################
# ---- integrated ----
# ##############################################################################


# 🔹 Zakładamy, że resolution_df jest już zdefiniowane jak wcześniej
# 🔹 i że masz obiekt integrated_pilotDexCortCtrl_filtered oraz listy markerów

# --- przygotowanie nazw rozdzielczości ---
resolutions <- gsub("seurat_clusters_res", "", resolution_df$resolution)

# --- 1️⃣ Uruchomienie ScType dla danych INTEGRATED ---
sctype_results_integrated_list <- map(
  resolutions,
  ~ combine_sctype_results(
    seurat_obj         = integrated_pilotDexCortCtrl_filtered,
    res                = .x,
    all_markers_list   = all_markers_list,
    brain_markers_list = brain_markers_list,
    assay              = "integrated"
  )
)

# --- nadanie nazw elementom listy ---
names(sctype_results_integrated_list) <- paste0("res", resolutions)

# ✅ lista ramek danych
sctype_results_integrated_list


# --- 2️⃣ Podsumowanie liczby unikalnych typów komórek ---
sctype_summary_integrated_df <- imap_dfr(
  sctype_results_integrated_list,
  ~ {
    df <- .x
    res <- .y
    
    tibble(
      resolution = res,
      n_clusters = nrow(df),
      n_scaled_allMarkers   = n_distinct(df[[str_subset(names(df), "_scaled_allMarkers")]]),
      n_noScaled_allMarkers = n_distinct(df[[str_subset(names(df), "_noScaled_allMarkers")]]),
      n_scaled_brain        = n_distinct(df[[str_subset(names(df), "_scaled_brain")]]),
      n_noScaled_brain      = n_distinct(df[[str_subset(names(df), "_noScaled_brain")]])
    )
  }
)

# --- 3️⃣ Posortuj po wartości numerycznej resolution ---
sctype_summary_integrated_df <- sctype_summary_integrated_df %>%
  mutate(res_numeric = as.numeric(str_replace(resolution, "res", ""))) %>%
  arrange(res_numeric) %>%
  select(-res_numeric)

# ✅ wynik zbiorczy
sctype_summary_integrated_df %>% as.data.frame()


# --- 4️⃣ Porównanie dwóch kolumn: allMarkers vs brain (dla wybranych resolution) ---
sctype_results_integrated_list[c("res0.3", "res0.13", "res0.5")] %>%
  lapply(., function(x) {
    x %>%
      select(cluster, contains("noScaled_allMarkers"), contains("noScaled_brain")) %>%
      set_colnames(c("cluster", "clustersAnnotationAllMarkers", "clustersAnnotationBrainMarkers")) %>%
      mutate(SameCellType = clustersAnnotationAllMarkers == clustersAnnotationBrainMarkers)
  }) -> sctype_integrated_comparison_list

# ✅ Przykład podglądu dla jednej rozdzielczości
sctype_integrated_comparison_list$res0.3



# --- przygotowanie macierzy ekspresji tak jak w combine_sctype_results ---
exprMat_cluster_meanData <- prepare_proteinCodingGenesMatrix_fromSeurat_v2(
  seurat_obj         = integrated_pilotDexCortCtrl_filtered,
  assay              = "integrated",
  layer              = "data",
  biomart_version    = "v112",
  method             = "mean",
  aggregate_cluster  = TRUE,
  group.by           = "seurat_clusters_res0.3"
)

# --- uruchom sctype_score tylko dla brain, bez skalowania ---
sctype_noscaled_brain_raw <- sctype_score(
  exprMat_cluster_meanData,
  scaled = FALSE,
  gs = brain_markers_list,
  gene_names_to_uppercase = TRUE
)


# --- uruchom sctype_score tylko dla brain, bez skalowania ---
sctype_noscaled_all_raw <- sctype_score(
  exprMat_cluster_meanData,
  scaled = FALSE,
  gs = all_markers_list,
  gene_names_to_uppercase = TRUE
)

sctype_scaled_all_raw <- sctype_score(
  exprMat_cluster_meanData,
  scaled = TRUE,
  gs = all_markers_list,
  gene_names_to_uppercase = TRUE
)


sctype_scaled_all_raw  %>% as.data.frame() %>% 
  rownames_to_column() %>% 
  arrange(desc(cluster_9))

sctype_noscaled_all_raw  %>% as.data.frame() %>% 
  rownames_to_column() %>% 
  arrange(desc(cluster_9))


for (res_name in c("res0.13", "res0.3", "res0.5")) {
  message("\n🔗 Przetwarzam ", res_name, " ...")
  
  # --- pobierz ramkę z wynikami ScType ---
  df_res <- sctype_results_integrated_list[[res_name]] %>%
    mutate(cluster_id = gsub("cluster_", "", cluster))
  
  join_by_df <- "cluster_id"
  join_col   <- paste0("seurat_clusters_", res_name)
  
  # --- sprawdź, czy kolumna istnieje w metadata ---
  if (!join_col %in% colnames(integrated_pilotDexCortCtrl_filtered@meta.data)) {
    warning("⚠️ Brak kolumny ", join_col, " w Seurat metadata, pomijam ", res_name)
    next
  }
  
  # --- dopasowanie typów danych ---
  df_res[[join_by_df]] <- as.character(df_res[[join_by_df]])
  integrated_pilotDexCortCtrl_filtered@meta.data[[join_col]] <-
    as.character(integrated_pilotDexCortCtrl_filtered@meta.data[[join_col]])
  
  # --- przygotuj dane do joinowania ---
  meta_tmp <- integrated_pilotDexCortCtrl_filtered@meta.data %>%
    rownames_to_column("cell_id")
  
  # dodaj tymczasową kolumnę z numerem klastra
  meta_tmp$tmp_cluster_id <- meta_tmp[[join_col]]
  
  meta_add <- meta_tmp %>%
    left_join(df_res, by = c("tmp_cluster_id" = "cluster_id")) %>%
    arrange(match(cell_id, rownames(integrated_pilotDexCortCtrl_filtered@meta.data))) %>%
    select(all_of(setdiff(colnames(df_res), c("cluster", "cluster_id")))) %>%
    as.data.frame()
  
  # --- dodaj kolumny do obiektu Seurat ---
  integrated_pilotDexCortCtrl_filtered <- AddMetaData(
    integrated_pilotDexCortCtrl_filtered,
    metadata = meta_add
  )
  
  # --- raport dopasowania ---
  n_matched <- sum(!is.na(meta_add[, 1]))
  n_total   <- nrow(meta_add)
  message("✅ ", res_name, ": dopasowano ", n_matched, "/", n_total, " komórek (", round(100 * n_matched / n_total, 1), "%)")
}
