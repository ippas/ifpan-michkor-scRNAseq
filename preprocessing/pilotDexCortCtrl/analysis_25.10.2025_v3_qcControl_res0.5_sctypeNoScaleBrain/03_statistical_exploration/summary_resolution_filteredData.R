library(dplyr)
library(tidyr)

integrated_pilotDexCortCtrl_filtered@meta.data %>%
  select(starts_with("seurat_clusters_res")) %>%
  summarise_all(~ n_distinct(.)) %>%
  pivot_longer(everything(), names_to = "resolution", values_to = "n_clusters") %>%
  arrange(resolution) %>% 
  as.data.frame()



integrated_pilotDexCortCtrl_filtered@meta.data %>%  head
  select(res0.3_RNA_noScaled_brain, treatment) %>% 
  rownames_to_column() %>% 
  group_by(res0.3_RNA_noScaled_brain, treatment) %>% 
  nest() %>% 
  mutate(n_cells = map_int(data, nrow)) %>% 
  select(-data) %>% 
  ungroup() %>% 
  pivot_wider(
    names_from = treatment,
    values_from = n_cells,
    values_fill = 0
  ) %>% 
  mutate(SUM = CTRL + DEX + CORT) %>% 
  arrange(desc(SUM))


integrated_pilotDexCortCtrl_filtered@meta.data %>%  
  select(res0.3_RNA_noScaled_brain, seurat_clusters_res0.3) %>% 
  rownames_to_column() %>% 
  select(-rowname) %>% 
  unique %>% group_by( res0.3_RNA_noScaled_brain) %>% 
  nest()
