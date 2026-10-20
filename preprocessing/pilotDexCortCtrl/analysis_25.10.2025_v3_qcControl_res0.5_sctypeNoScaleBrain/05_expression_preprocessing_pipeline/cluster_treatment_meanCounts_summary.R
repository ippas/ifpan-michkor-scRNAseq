rawSumCountsTreatment_res0.3_RNA_noScaledBrain <- aggregate_counts_by_cluster_sample(
  seurat_obj  = integrated_pilotDexCortCtrl_filtered,
  cluster_col = "res0.3_RNA_noScaled_brain",
  sample_col  = "treatment",
  assay       = "RNA",
  layer       = "counts",  # <- Seurat v5
  agg         = "sum",
  return      = "list",
  genes_keep  = hgnc_symbols_vector_v110,
  progress    = TRUE
)

integrated_pilotDexCortCtrl_filtered@meta.data %>%
  group_by(treatment, res0.3_RNA_noScaled_brain) %>%
  nest() %>% 
  mutate(n_cells = map_int(data, nrow)) %>% 
  select(-data) %>% 
  # mutate(cluster = paste0("cluster_", res0.3_RNA_noScaled_brain)) %>%
  mutate(cluster = paste0("cluster_", res0.3_RNA_noScaled_brain)) %>%
  ungroup %>% 
  select(treatment, cluster, n_cells) %>%
  left_join(
    rawSumCountsTreatment_res0.3_RNA_noScaledBrain %>%
      lapply(function(x) {
        colSums(x) %>%
          as.data.frame() %>%
          tibble::rownames_to_column("treatment") %>%
          dplyr::rename(count_sum = 2)
      }) %>%
      bind_rows(.id = "cluster") %>%
      select(cluster, treatment, count_sum),
    by = c("treatment", "cluster")
  ) %>% 
  mutate(mean_count = count_sum/n_cells) %>% 
  mutate(cluster = str_replace_all(cluster, "cluster_", "")) -> mean_counts_per_clusterTreatment_df


mean_counts_per_clusterTreatment_df %>% filter(cluster == "Astrocytes")

