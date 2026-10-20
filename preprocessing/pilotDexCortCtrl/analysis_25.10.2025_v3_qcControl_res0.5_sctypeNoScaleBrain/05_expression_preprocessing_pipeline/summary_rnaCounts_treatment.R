integrated_pilotDexCortCtrl_filtered@meta.data %>% 
  filter(treatment %in% c("CTRL", "CORT", "DEX")) %>% 
  mutate(
    treatment = factor(treatment, levels = c("CTRL", "CORT", "DEX")),
    log_nCount = log2(nCount_RNA)
  ) %>% 
  group_by(treatment) %>% 
  nest() %>% 
  mutate(n_cells = map(data, ~ .x %>% nrow)) %>% 
  mutate(sum_counts = map(data, ~ .x$nCount_RNA %>% sum)) %>% 
  mutate(mean = map(data, ~ .x$nCount_RNA %>% mean)) %>% 
  mutate(median = map(data, ~ .x$nCount_RNA %>% median)) %>% 
  mutate(sd = map(data, ~ .x$nCount_RNA %>% sd)) %>% 
  select(-data) %>% 
  unnest %>% 
  as.data.frame() %>% 
  mutate(se = sd/sqrt(n_cells))

  summarise(
    n_cells   = n(),
    sum_reads = sum(nCount_RNA, na.rm = TRUE),
    mean      = mean(log_nCount, na.rm = TRUE),
    sd        = sd(log_nCount, na.rm = TRUE),
    se        = sd / sqrt(n_cells),
    median    = median(log_nCount, na.rm = TRUE)
  )



integrated_pilotDexCortCtrl_filtered@meta.data %>%
  select(nCount_RNA, treatment) %>% 
  group_by(treatment) %>% 
  nest %>% 
  mutate(mean = map(data, ~ .x$nCount_RNA %>% mean)) %>% 
  mutate(median = map(data, ~ .x$nCount_RNA %>% median)) %>% 
  unnest(mean, median)

effectsize::cohens_d(nCount_RNA ~ treatment, data = integrated_pilotDexCortCtrl_filtered@meta.data %>% filter(treatment != "DEX"))
